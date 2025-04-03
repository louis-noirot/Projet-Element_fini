/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Utilisation de l'API de GMSH pour cr�er un maillage
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#ifdef TEST1
#include "elast_glfem.h"
#include "elast_fem.h"
#endif

#ifdef TEST2
#include "slv_glfem.h"
#include "slv_fem.h"
#endif

#ifdef TEST3
#include "elast_glfem.h"
#include "elast_fem.h"
#endif

#include <time.h>


int main(void){ 

    #ifdef TEST1
    printf("\n\n    V : Mesh and size mesh field \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n");

    double Lx = 20.0;
    double Ly = 55.0;
    int ierr;
       
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
     
    theGeometry->LxPlate     =  Lx;
    theGeometry->LyPlate     =  Ly;     
    theGeometry->h           =  Lx * 0.075;    
    theGeometry->elementType = FEM_TRIANGLE;
    //theGeometry->elementType = FEM_QUAD;
   
    geoMeshGenerate();
    geoMeshImport();

    geoSetDomainName(0,"2_beam_T");
    geoSetDomainName(1,"R_column_3");
    geoSetDomainName(2,"3_beam_B");
    geoSetDomainName(3,"L_column_3");
    geoSetDomainName(4,"L_column");
    geoSetDomainName(5,"3_beam_T");
    geoSetDomainName(6,"Bottom1");
    geoSetDomainName(7,"R_column");
    geoSetDomainName(8,"L_column_1");
    geoSetDomainName(9,"Bottom2");
    geoSetDomainName(10,"1_beam_B");
    geoSetDomainName(11,"R_column_1");
    geoSetDomainName(12,"1_beam_T");
    geoSetDomainName(13,"L_column_2");
    geoSetDomainName(14,"R_column_2");
    geoSetDomainName(15,"2_beam_B");
    geoMeshWrite("../data/mesh.txt");
    

//
//  -2- Creation du fichier du maillage
//
    char filename[] = "../data/mesh.txt";
    geoMeshWrite(filename);

//
//  -3- Champ de la taille de r�f�rence du maillage
//
    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);  
    double hMax = femMax(meshSizeField,theNodes->nNodes);  
    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);
 
//
//  -4- Visualisation du maillage
//  
    int mode = 1; // Change mode by pressing "j", "k", "l"
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[256];
    double pos[2] = {20,460};
 
 
    GLFWwindow* window = glfemInit("EPL1110 : Mesh generation ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
    //    glfemChangeState(&mode, theMeshes->nMesh);
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}

        if (t-told > 0.5) {freezingButton = FALSE; }

        if (mode == 1) {
            glfemPlotField(theGeometry->theElements, meshSizeField);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
            }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 

            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);

            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
            }
            
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);  
    geoFinalize();
    glfwTerminate(); 
    exit(EXIT_SUCCESS);
    #endif

    #ifdef TEST2
    printf("\n\n    V : Plot results \n");
    printf("    S : Spy matrix \n");
    printf("    F-B-I : Full solver - Band solver - Iterative solver \n");
    printf("    X-Y-N : Renumbering along x - along y - No renumbering \n");
    
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_YNUM;
    char meshFileName[] = "../data/mesh.txt";

    femDiffusionProblem* theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
    clock_t tic = clock();
    femDiffusionCompute(theProblem);  
    femSolverPrintInfos(theProblem->solver);
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);
    
    int option = 1;    
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;
    femMesh *theMesh = theProblem->geo->theElements;

    GLFWwindow* window = glfemInit("LEPL1110 : Band Solver ");
    glfwMakeContextCurrent(window);
     
    do {int testConvergence,w,h;
        char theMessage[256];
        sprintf(theMessage, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        glfwGetFramebufferSize(window,&w,&h);
      
        if (option == 1) {
            glfemReshapeWindows(theMesh,w,h);
            glfemPlotField(theMesh,theProblem->soluce);   }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); }
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
      
        if (solverType != newSolverType || renumType != newRenumType) { 
            solverType = newSolverType;
            renumType = newRenumType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            theMesh = theProblem->geo->theElements;
            clock_t tic = clock();
            do {
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; }
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); }
        if (glfwGetKey(window,'V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey(window,'F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey(window,'B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        if (glfwGetKey(window,'I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
       
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed
    glfwTerminate(); 
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
    #endif
    
    #ifdef TEST3
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n\n\n");
 
    double Lx = 20.0;
    double Ly = 55.0;
       
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
     
    theGeometry->LxPlate     =  Lx;
    theGeometry->LyPlate     =  Ly;     
    theGeometry->h           =  Lx * 0.075;    
    theGeometry->elementType = FEM_TRIANGLE;
    //theGeometry->elementType = FEM_QUAD;
   
    geoMeshGenerate();
    geoMeshImport();
    geoSetDomainName(1,"Symmetry");

    geoSetDomainName(0,"2_beam_T");
    geoSetDomainName(1,"R_column_3");
    geoSetDomainName(2,"3_beam_B");
    geoSetDomainName(3,"L_column_3");
    geoSetDomainName(4,"L_column");
    geoSetDomainName(5,"3_beam_T");
    geoSetDomainName(6,"Bottom1");
    geoSetDomainName(7,"R_column");
    geoSetDomainName(8,"L_column_1");
    geoSetDomainName(9,"Bottom2");
    geoSetDomainName(10,"1_beam_B");
    geoSetDomainName(11,"R_column_1");
    geoSetDomainName(12,"1_beam_T");
    geoSetDomainName(13,"L_column_2");
    geoSetDomainName(14,"R_column_2");
    geoSetDomainName(15,"2_beam_B");
    geoMeshWrite("../data/mesh.txt");
         
//
//  -2- Creation probleme 
//
    double E   = 211.e9;
    double nu  = 0.3;
    double rho = 7.85e3; 
    double g   = 9.81;

    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN);
    femElasticityAddBoundaryCondition(theProblem,"R_column_1",DIRICHLET_X,-1.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom1",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom1",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom2",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom2",DIRICHLET_X,0.0);
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem); 
    
//
//  -3- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e0;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
     
    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + theSoluce[2*i+1]*theSoluce[2*i+1]); }
   
    double hMin = femMin(normDisplacement,theNodes->nNodes);  
    double hMax = femMax(normDisplacement,theNodes->nNodes);  
    printf(" ==== Minimum displacement          : %14.7e \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e \n",hMax);
  
//
//  -4- Visualisation du maillage
//  
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
 
    GLFWwindow* window = glfemInit("EPL1110 : Linear elasticity ");
    glfwMakeContextCurrent(window);

     do {int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        
        if (t-told > 0.5) {freezingButton = FALSE; }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);  }

            glfwSwapBuffers(window);
            glfwPollEvents();
    } while(glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
            glfwWindowShouldClose(window) != 1 );
             
    // Check if the ESC key was pressed or the window was closed
    free(normDisplacement);
    femElasticityFree(theProblem) ; 
    geoFinalize();
    glfwTerminate(); 
    exit(EXIT_SUCCESS);
    return 0;
    #endif

    return 0;
}