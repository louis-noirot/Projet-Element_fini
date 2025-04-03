#include "fem.h"
#include <math.h>

double geoSize(double x, double y) {
    femGeo *theGeometry = geoGetGeometry();
    double h = theGeometry->h;

    // Distance à chaque trou rectangulaire
    double dMin = 1e10;

    // Rectangle du haut : y ∈ [36, 53], x ∈ [2, 18]
    if (x >= 2 && x <= 18) {
        if (y > 53) dMin = fmin(dMin, y - 53);
        else if (y < 36) dMin = fmin(dMin, 36 - y);
        else dMin = 0;
    } else {
        double dx = (x < 2) ? 2 - x : x - 18;
        if (y >= 36 && y <= 53) dMin = fmin(dMin, dx);
        else {
            double dy = (y < 36) ? 36 - y : y - 53;
            dMin = fmin(dMin, sqrt(dx*dx + dy*dy));
        }
    }

    // Rectangle du milieu : y ∈ [16, 34], x ∈ [2, 18]
    if (x >= 2 && x <= 18) {
        if (y > 34) dMin = fmin(dMin, y - 34);
        else if (y < 16) dMin = fmin(dMin, 16 - y);
        else dMin = 0;
    } else {
        double dx = (x < 2) ? 2 - x : x - 18;
        if (y >= 16 && y <= 34) dMin = fmin(dMin, dx);
        else {
            double dy = (y < 16) ? 16 - y : y - 34;
            dMin = fmin(dMin, sqrt(dx*dx + dy*dy));
        }
    }

    // Rectangle du bas : y ∈ [0, 14], x ∈ [2, 18]
    if (x >= 2 && x <= 18) {
        if (y > 14) dMin = fmin(dMin, y - 14);
        else if (y < 0) dMin = fmin(dMin, 0 - y);
        else dMin = 0;
    } else {
        double dx = (x < 2) ? 2 - x : x - 18;
        if (y >= 0 && y <= 14) dMin = fmin(dMin, dx);
        else {
            double dy = (y < 0) ? 0 - y : y - 14;
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
    femGeo *theGeometry = geoGetGeometry();
    int ierr;

    double w = 20.0;
    double h = 55.0;

    // Boîte externe
    int idMain = gmshModelOccAddRectangle(0, 0, 0, w, h, -1, 0.0, &ierr);

    // Trois rectangles internes (découpes)
    int idTop    = gmshModelOccAddRectangle(2, 36, 0, 16, 17, -1, 0.0, &ierr);
    int idMiddle = gmshModelOccAddRectangle(2, 16, 0, 16, 18, -1, 0.0, &ierr);
    int idBottom = gmshModelOccAddRectangle(2, 0, 0, 16, 14, -1, 0.0, &ierr);

    int mainBox[] = {2, idMain};
    int cuts[]    = {2, idTop, 2, idMiddle, 2, idBottom};

    gmshModelOccCut(mainBox, 2, cuts, 6, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccSynchronize(&ierr);

    // Appliquer la taille locale via callback
    geoSetSizeCallback(geoSize);

    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);

    /*gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 11, &ierr);     // Auparavant = 8 (=11 dixit Alexandre)
    gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);    // Nouvelle option magique
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  
    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
    gmshModelMeshGenerate(2, &ierr);*/

    /*gmshFltkInitialize(&ierr);
    gmshFltkRun(&ierr);*/
}
