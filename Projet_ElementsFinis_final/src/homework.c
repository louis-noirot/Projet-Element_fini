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
    //femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int map[4], mapX[4], mapY[4];
    int iElem, iInteg, i, j;
    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    double rho = theProblem->rho;
    double g = theProblem->g;
    int nLocal = theProblem->space->n;
    femElasticCase iCase = theProblem->planarStrainStress;

    int *number = theNodes->number;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j = 0; j < nLocal; j++) {
            int node = theMesh->elem[iElem*nLocal + j];
            int renumNode = (theProblem->typeSolver == FEM_BAND) ? number[node] : node;
    
            map[j]  = renumNode;
            mapX[j] = 2 * renumNode;
            mapY[j] = 2 * renumNode + 1;
    
            x[j] = theNodes->X[node];
            y[j] = theNodes->Y[node];
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
                        double sigTheta = phi[i] * phi[j] / (r * r);
                        femSystemAdd(theProblem, mapX[i], mapX[j], (a*dphidx[i]*dphidx[j] + c*dphidy[i]*dphidy[j] + sigTheta*b) * dV);
                        femSystemAdd(theProblem, mapX[i], mapY[j], ((a-c)*dphidx[i]*a*dphidx[j] + b*dphidy[i]*dphidy[j]) * dV);
                        femSystemAdd(theProblem, mapY[i], mapX[j], ((a-c)*dphidx[i]*a*dphidx[j] + b*dphidy[i]*dphidy[j]) * dV);
                        femSystemAdd(theProblem, mapY[i], mapY[j], (a*dphidx[i]*dphidx[j] + c*dphidy[i]*dphidy[j]) * dV);
                    }
                    //B[mapY[i]] -= phi[i] * rho * g * dV;
                    femSystemAddRightHandSide(theProblem, mapY[i], -(phi[i] * rho * g * dV));

                }
            }
            else { // PLANAR_STRAIN ou PLANAR_STRESS
                double dV = jac * weight;
                for (i = 0; i < nLocal; i++) {
                    for (j = 0; j < nLocal; j++) {
                        // xx
                        femSystemAdd(theProblem, mapX[i], mapX[j], (a*dphidx[i]*dphidx[j] + c*dphidy[i]*dphidy[j]) * dV);
                        // xy
                        femSystemAdd(theProblem, mapX[i], mapY[j], ((a-c)*dphidx[i]*dphidy[j] + b*dphidy[i]*dphidx[j]) * dV);
                        // yx
                        femSystemAdd(theProblem, mapY[i], mapX[j], ((a-c)*dphidy[i]*dphidx[j] + b*dphidx[i]*dphidy[j]) * dV);
                        // yy
                        femSystemAdd(theProblem, mapY[i], mapY[j], (a*dphidy[i]*dphidy[j] + c*dphidx[i]*dphidx[j]) * dV);
                    } femSystemAddRightHandSide(theProblem, mapY[i], -(phi[i] * rho * g * dV));
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
                    femSystemAddRightHandSide(theProblem, mapU[i], dS * phi[i] * value);
                }
            }
        }
    }
}

double *femElasticitySolve(femProblem *theProblem) {
    femFullSystem  *theSystem = theProblem->system;
    femFullSystemInit(theSystem);                   // Remise à zéro
    femElasticityAssembleElements(theProblem);      // Assemblage des éléments
    femElasticityAssembleNeumann(theProblem);       // Conditions de Neumann

    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femSystemConstraint(theProblem, i, value);
        }
    }
    if (theProblem->typeSolver == FEM_NONE) {
        femNoneSystem *s = (femNoneSystem *)theProblem->systemNone;
        femNoneSystemEliminate(s);
        memcpy(theProblem->soluce, s->X, s->size * sizeof(double));
    } else if (theProblem->typeSolver == FEM_BAND) {
        femBandSystem *s = theProblem->systemBand;
        double *solution = femBandSystemEliminate(s);
        memcpy(theProblem->soluce, solution, s->size * sizeof(double));
    } else {
        femFullSystemEliminate(theSystem);
        memcpy(theProblem->soluce, theSystem->B, theSystem->size * sizeof(double));
    }
    

    return theProblem->soluce;
}

double *femElasticityForces(femProblem *theProblem) {
    int size = theProblem->geometry->theNodes->nNodes * 2;
    double *soluce = theProblem->soluce;
    double *residual = theProblem->residuals;

    // Initialisation du système et assemblage
    if (theProblem->typeSolver == FEM_FULL) {
        femFullSystem *theSystem = theProblem->system;
        femFullSystemInit(theSystem);
        femElasticityAssembleElements(theProblem);
        femElasticityAssembleNeumann(theProblem);

        for (int i = 0; i < size; i++) {
            residual[i] = -theSystem->B[i];
            for (int j = 0; j < size; j++) {
                residual[i] += theSystem->A[i][j] * soluce[j];
            }
        }

    } else if (theProblem->typeSolver == FEM_BAND) {
        femBandSystem *theSystem = theProblem->systemBand;
        femBandSystemInit(theSystem);
        femElasticityAssembleElements(theProblem);
        femElasticityAssembleNeumann(theProblem);

        int band = theSystem->band;
        for (int i = 0; i < size; i++) {
            residual[i] = -theSystem->B[i];
            for (int j = fmax(0, i - band); j <= fmin(size - 1, i + band); j++) {
                int colBandIndex = j - i + band;
                if (colBandIndex >= 0 && colBandIndex < 2 * band + 1)
                    residual[i] += theSystem->A[i][colBandIndex] * soluce[j];
            }
        }

    } else if (theProblem->typeSolver == FEM_NONE) {
        femNoneSystem *theSystem = theProblem->systemNone;
        for (int i = 0; i < size; i++) {
            residual[i] = -theSystem->B[i];
            for (int j = 0; j < size; j++) {
                residual[i] += theSystem->A[i][j] * soluce[j];
            }
        }
    } else {
        Error("Type de solveur inconnu dans femElasticityForces");
    } return residual;
}

void femElasticityComputeStresses(femProblem *theProblem) {
    femGeo         *theGeometry = theProblem->geometry;
    femMesh        *theMesh     = theGeometry->theElements;
    femNodes       *theNodes    = theGeometry->theNodes;
    femIntegration *theRule     = theProblem->rule;
    femDiscrete    *theSpace    = theProblem->space;
    double         *sol         = theProblem->soluce;

    int nNodes = theNodes->nNodes;
    int nLocal = theMesh->nLocalNode;
    int nElem  = theMesh->nElem;

    // Allocation si nécessaire
    if (theProblem->sigmaXX == NULL) {
        theProblem->sigmaXX = malloc(nNodes * sizeof(double));
        theProblem->sigmaYY = malloc(nNodes * sizeof(double));
        theProblem->sigmaXY = malloc(nNodes * sizeof(double));
        for (int i = 0; i < nNodes; ++i) {
            theProblem->sigmaXX[i] = 0.0;
            theProblem->sigmaYY[i] = 0.0;
            theProblem->sigmaXY[i] = 0.0;
        }
    }

    // Pour moyenne pondérée
    int *count = calloc(nNodes, sizeof(int));

    // Constantes du matériau
    double A = theProblem->A;
    double B = theProblem->B;
    double C = theProblem->C;

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];

    for (int iElem = 0; iElem < nElem; ++iElem) {

        int map[4], mapX[4], mapY[4];
        for (int i = 0; i < nLocal; i++) {
            map[i]  = theMesh->elem[iElem * nLocal + i];
            x[i]    = theNodes->X[map[i]];
            y[i]    = theNodes->Y[map[i]];
            mapX[i] = 2 * map[i];
            mapY[i] = 2 * map[i] + 1;
        }

        for (int iInteg = 0; iInteg < theRule->n; ++iInteg) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
            for (int i = 0; i < nLocal; ++i) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }

            double jac = dxdxsi * dydeta - dxdeta * dydxsi;
            for (int i = 0; i < nLocal; ++i) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            double u[4], v[4];
            for (int i = 0; i < nLocal; ++i) {
                u[i] = sol[mapX[i]];
                v[i] = sol[mapY[i]];
            }

            double eps_xx = 0.0, eps_yy = 0.0, eps_xy = 0.0;
            for (int i = 0; i < nLocal; ++i) {
                eps_xx += u[i] * dphidx[i];
                eps_yy += v[i] * dphidy[i];
                eps_xy += 0.5 * (u[i] * dphidy[i] + v[i] * dphidx[i]);
            }

            double sigma_xx = A * eps_xx + B * eps_yy;
            double sigma_yy = B * eps_xx + A * eps_yy;
            double sigma_xy = 2.0 * C * eps_xy;

            // Contribution aux nœuds
            for (int i = 0; i < nLocal; ++i) {
                int node = map[i];
                theProblem->sigmaXX[node] += sigma_xx;
                theProblem->sigmaYY[node] += sigma_yy;
                theProblem->sigmaXY[node] += sigma_xy;
                count[node]++;
            }
        }
    }

    // Moyenne sur les contributions
    for (int i = 0; i < nNodes; ++i) {
        if (count[i] > 0) {
            theProblem->sigmaXX[i] /= count[i];
            theProblem->sigmaYY[i] /= count[i];
            theProblem->sigmaXY[i] /= count[i];
        }
    }

    free(count);
}

