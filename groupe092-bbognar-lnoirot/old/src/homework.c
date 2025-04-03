#include "fem.h"

double geoSize(double x, double y){
    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;
      
    double hfinal = h;
    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - r0;

    if (d < d0) {
        double a = (-2*h + 2*h0)/(d0*d0*d0);
        double b = (3*h  - 3*h0)/(d0*d0);
        double c = 0;
        hfinal = a*d*d*d + b*d*d + c*d + h0;
    }
    d = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1)) - r1;
    if (d < d1) {
        double a = (-2*h + 2*h1)/(d1*d1*d1);
        double b = (3*h  - 3*h1)/(d1*d1);
        double c = 0;
        hfinal = fmin(hfinal,a*d*d*d + b*d*d + c*d + h1);
    } return hfinal;
}

void geoMeshGenerate() {
    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(-w/2.0,-h/2.0,0.0,w,h,-1,0.0,&ierr); 
    int idNotch = gmshModelOccAddDisk(x0,y0,0.0,r0,r0,-1,NULL,0,NULL,0,&ierr); 
    int idHole  = gmshModelOccAddDisk(x1,y1,0.0,r1,r1,-1,NULL,0,NULL,0,&ierr); 
    
    int plate[] = {2,idPlate};
    int notch[] = {2,idNotch};
    int hole[]  = {2,idHole};
    gmshModelOccCut(plate,2,notch,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(plate,2,hole ,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
                  
    geoSetSizeCallback(geoSize);   
    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads (avec quelques triangles...) :-)
//
    /*gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 11, &ierr);     // Auparavant = 8 (=11 dixit Alexandre)
    gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);    // Nouvelle option magique
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  
    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
    gmshModelMeshGenerate(2, &ierr);*/
   
    //Plot avec Fltk et vous avez accès au graphique de gmsh :-)
    /*gmshFltkInitialize(&ierr);
    gmshFltkRun(&ierr);
    chk(ierr);*/
}