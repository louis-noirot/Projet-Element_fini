#include "fem.h"

double geoSize(double x, double y) {
    femGeo *theGeometry = geoGetGeometry();
    double h = theGeometry->LyPlate *0.075;
    double w = theGeometry->LxPlate;
    double Ly = theGeometry->LyPlate;
    double b = Ly * (2.0 / 55.0);   // beam
    double c = w / 10.0;         // column
    double h1 = (14.0 / 55.0) * Ly; // height of the bottom rectangle
    double h2 = (18.0 / 55.0) * Ly; // height of the middle rectangle
    double h3 = (17.0 / 55.0) * Ly; // height of the top rectangle

    // Distance à chaque trou rectangulaire
    double dMin = 1e10;

    double h30 = h1 + b + h2 + b;
    // Rectangle du haut : y ∈ [h30, h30+h3], x ∈ [c, w-c]
    if (x >= c && x <= w-c) {
        if (y > h30+h3) dMin = fmin(dMin, y - (h30+h3));
        else if (y < h30) dMin = fmin(dMin, h30 - y);
        else dMin = 0;
    } else {
        double dx = (x < c) ? c - x : x - (w-c);
        if (y >= h30 && y <= h30+h3) dMin = fmin(dMin, dx);
        else {
            double dy = (y < h30) ? h30 - y : y - (h30+h3);
            dMin = fmin(dMin, sqrt(dx*dx + dy*dy));
        }
    }

    double h20 = h1 + b;
    // Rectangle du milieu : y ∈ [h20, h20+h2], x ∈ [c, w-c]
    if (x >= c && x <= w-c) {
        if (y > h20+h2) dMin = fmin(dMin, y - (h20+h2));
        else if (y < h20) dMin = fmin(dMin, h20 - y);
        else dMin = 0;
    } else {
        double dx = (x < c) ? c - x : x - (w-c);
        if (y >= h20 && y <= h20+h2) dMin = fmin(dMin, dx);
        else {
            double dy = (y < h20) ? h20 - y : y - (h20+h2);
            dMin = fmin(dMin, sqrt(dx*dx + dy*dy));
        }
    }

    // Rectangle du bas : y ∈ [0, h1], x ∈ [c, w-c]
    if (x >= c && x <= w-c) {
        if (y > h1) dMin = fmin(dMin, y - h1);
        else if (y < 0) dMin = fmin(dMin, 0 - y);
        else dMin = 0;
    } else {
        double dx = (x < c) ? c - x : x - (w-c);
        if (y >= 0 && y <= h1) dMin = fmin(dMin, dx);
        else {
            double dy = (y < 0) ? 0 - y : y - h1;
            dMin = fmin(dMin, sqrt(dx*dx + dy*dy));
        }
    }
    // Loi de raffinement douce
    double d = fmax(dMin, 1e-2);  // éviter div par 0
    double h0 = h * 0.2;
    double hfinal = h0 + (h - h0) * tanh(2.0 * d);

    return hfinal;
}

void geoMeshGenerate() {
    femGeo* theGeometry = geoGetGeometry();
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double b = h*(2.0/55.0);   // beam
    double c = w/10.0;         // column
    double h1 = (14.0/55.0)*h; // height of the bottom rectangle
    double h2 = (18.0/55.0)*h; // height of the middle rectangle
    double h3 = (17.0/55.0)*h; // height of the top rectangle
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr);
    int idBotom = gmshModelOccAddRectangle(c,0.0,0.0,w-2*c,h1,-1,0.0,&ierr);
    int idMiddle = gmshModelOccAddRectangle(c,b+h1,0.0,w-2*c,h2,-1,0.0,&ierr);
    int idTop = gmshModelOccAddRectangle(c,b*2+h1+h2,0.0,w-2*c,h3,-1,0.0,&ierr);
    int rect[] = {2,idRect};
    int bottom[] = {2,idBotom};
    int middle[] = {2,idMiddle};
    int top[] = {2,idTop};

    gmshModelOccCut(rect,2,bottom,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(rect,2,middle,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(rect,2,top,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    //geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",11,&ierr);  
        gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);
    }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);
    } return;
}

void femElasticityAssembleElements(femProblem *theProblem) {
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int map[4], mapX[4], mapY[4];
    int iElem, iInteg, i, j;
    int nLocal = theMesh->nLocalNode;
    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    double rho = theProblem->rho;
    double g = theProblem->g;
    double **A = theSystem->A;
    double *B = theSystem->B;
    femElasticCase iCase = theProblem->planarStrainStress;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {

        for (j = 0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal + j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        }

        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i]*dphidxsi[i];
                dxdeta += x[i]*dphideta[i];
                dydxsi += y[i]*dphidxsi[i];
                dydeta += y[i]*dphideta[i];
            }

            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            // AXISYM or 2D Planar ?
            if (iCase == AXISYM) {
                double r = 0.0;
                for (i = 0; i < theSpace->n; i++)
                    r += phi[i] * x[i];  // r = rayon local (x)

                double dV = 2.0 * M_PI * r * jac * weight;

                for (i = 0; i < nLocal; i++) {
                    for (j = 0; j < nLocal; j++) {
                        double sigTheta = phi[i] * phi[j] / r * r;
                        A[mapX[i]][mapX[j]] += (dphidx[i]*a*dphidx[j] + dphidy[i]*c*dphidy[j] + sigTheta*b) * dV;
                        A[mapX[i]][mapY[j]] += (dphidx[i]*b*dphidy[j] + dphidy[i]*c*dphidx[j]) * dV;
                        A[mapY[i]][mapX[j]] += (dphidy[i]*b*dphidx[j] + dphidx[i]*c*dphidy[j]) * dV;
                        A[mapY[i]][mapY[j]] += (dphidy[i]*a*dphidy[j] + dphidx[i]*c*dphidx[j]) * dV;
                    }
                    B[mapY[i]] -= phi[i] * rho * g * dV;
                }
            }
            else { // PLANAR_STRAIN ou PLANAR_STRESS
                double dV = jac * weight;
                for (i = 0; i < nLocal; i++) {
                    for (j = 0; j < nLocal; j++) {
                        A[mapX[i]][mapX[j]] += (dphidx[i]*a*dphidx[j] + dphidy[i]*c*dphidy[j]) * dV;
                        A[mapX[i]][mapY[j]] += (dphidx[i]*b*dphidy[j] + dphidy[i]*c*dphidx[j]) * dV;
                        A[mapY[i]][mapX[j]] += (dphidy[i]*b*dphidx[j] + dphidx[i]*c*dphidy[j]) * dV;
                        A[mapY[i]][mapY[j]] += (dphidy[i]*a*dphidy[j] + dphidx[i]*c*dphidx[j]) * dV;
                    }
                    B[mapY[i]] -= phi[i] * rho * g * dV;
                }
            }
        }
    }
}

void femElasticityAssembleNeumann(femProblem *theProblem) {
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule   = theProblem->ruleEdge;
    femDiscrete    *theSpace  = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    femElasticCase iCase     = theProblem->planarStrainStress;

    double x[2], y[2], phi[2];
    int iBnd, iElem, iInteg, iEdge, i, j;
    int map[2], mapU[2];
    int nLocal = 2;
    double *B = theSystem->B;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        double value = theCondition->value;

        int shift = -1;
        if (type == NEUMANN_X) shift = 0;
        if (type == NEUMANN_Y) shift = 1;
        if (shift == -1) continue;

        for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
            iElem = theCondition->domain->elem[iEdge];
            for (j = 0; j < nLocal; j++) {
                map[j]  = theEdges->elem[iElem * nLocal + j];
                mapU[j] = 2 * map[j] + shift;
                x[j]    = theNodes->X[map[j]];
                y[j]    = theNodes->Y[map[j]];
            }

            double length = sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]));
            double jac = length / 2.0;

            for (iInteg = 0; iInteg < theRule->n; iInteg++) {
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];
                femDiscretePhi(theSpace, xsi, phi);

                // Calcul du rayon local moyen (x correspond à r)
                double r = phi[0] * x[0] + phi[1] * x[1];
                double dS;

                if (iCase == AXISYM) {
                    dS = 2.0 * M_PI * r * jac * weight;
                } else {
                    dS = jac * weight;
                }

                for (i = 0; i < nLocal; i++) {
                    B[mapU[i]] += dS * phi[i] * value;
                }
            }
        }
    }
}


double *femElasticitySolve(femProblem *theProblem) {
    femFullSystem  *theSystem = theProblem->system;
    femFullSystemInit(theSystem);
    femElasticityAssembleElements(theProblem);    
    femElasticityAssembleNeumann(theProblem);
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value);
        }
    }
    femFullSystemEliminate(theSystem);
    memcpy(theProblem->soluce, theSystem->B, theSystem->size * sizeof(double));
    return theProblem->soluce;
}

double * femElasticityForces(femProblem *theProblem) {
    femFullSystem  *theSystem = theProblem->system;
    femFullSystemInit(theSystem);
    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);         
    double *theResidual = theProblem->residuals;
    for(int i=0; i < theSystem->size; i++) {
        theResidual[i] = 0.0;}
    for(int i=0; i < theSystem->size; i++){
        for(int j=0; j < theSystem->size; j++) {
            theResidual[i] += theSystem->A[i][j] * theProblem->soluce[j];
        }
        theResidual[i] -= theSystem->B[i];
    } return theResidual;
}