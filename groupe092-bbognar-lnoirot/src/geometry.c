#include "elast_fem.h"
#include <math.h>

double geoSize(double x, double y) {
    femGeo *theGeometry = geoGetGeometry();
    double h = theGeometry->h;
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
    printf("Lx = %f, Ly = %f\n",w,h);
    
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
    geoSetSizeCallback(geoSize);
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