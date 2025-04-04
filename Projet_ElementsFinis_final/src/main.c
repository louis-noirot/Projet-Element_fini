/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *  Calcul des densités de force aux noeuds contraints
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"
# include <time.h>

femElementType elementType  = FEM_TRIANGLE;      // FEM_QUAD or FEM_TRIANGLE
femElasticCase theCase      = PLANAR_STRESS; // PLANAR_STRESS or PLANAR_STRAIN or AXISYM (pas utile dans notre cas mais implémenté)
femSolverType theSolver    = FEM_BAND;     // FEM_FULL or FEM_BAND or FEM_NONE (ne fait rien)
femRenumType renumType      = FEM_NO;  // FEM_NO or FEM_XNUM or FEM_YNUM

double zoomFactor = 1.0;
double panX = 0.0, panY = 0.0;
double lastX = 0.0, lastY = 0.0;
int isDragging = 0;


void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    if (yoffset > 0) zoomFactor *= 1.1;
    if (yoffset < 0) zoomFactor /= 1.1;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            isDragging = 1;
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            lastX = xpos;
            lastY = ypos;
        } else if (action == GLFW_RELEASE) {
            isDragging = 0;
        }
    }
}

double fun(double x, double y) {
    return 1;
}

int main(void) {  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    X : Horizontal residuals for unconstrained equations \n");
    printf("    Y : Horizontal residuals for unconstrained equations \n");
    printf("    N : Next domain highlighted\n");
    printf("    1 : Champs sigmaXX \n");
    printf("    2 : Champs sigmaYY \n");
    printf("    3 : Champs sigmaXY \n");
    printf("    U : Champs de la norme du déplacement \n");

    double Lx = 20.0;
    double Ly = 55.0;
      
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    
    theGeometry->LxPlate     =  Lx;
    theGeometry->LyPlate     =  Ly;     
    theGeometry->h           =  Lx * 0.075;
    theGeometry->elementType = elementType;
  
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
//  -2- Creation probleme 
//
    
    double E   = 200.e9;
    double nu  = 0.3;
    double rho = 7.85e3; 
    double g   = 9.81;
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,theCase, theSolver, renumType);
    femElasticityAddBoundaryCondition(theProblem,"R_column_1",NEUMANN_X,-12000000);
    femElasticityAddBoundaryCondition(theProblem,"3_beam_T",NEUMANN_Y,-40875000/2);
    femElasticityAddBoundaryCondition(theProblem,"Bottom1",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom1",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom2",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom2",DIRICHLET_Y,0.0);
    femElasticityPrint(theProblem);

//
//  -3- Resolution du probleme et calcul des forces
//

    clock_t tStart = clock();
    double *theSoluce = femElasticitySolve(theProblem);
    clock_t tSolve = clock();

    femElasticityComputeStresses(theProblem);
    clock_t tStresses = clock();

    double *theForces = femElasticityForces(theProblem);
    clock_t tForces = clock();

    double area = femElasticityIntegrate(theProblem, fun);
    clock_t tIntegrate = clock();

    double clk = (double)CLOCKS_PER_SEC;
    printf("\n⏱️  Benchmark (%s solver) :\n", 
        theProblem->typeSolver == FEM_FULL ? "FULL" :
        theProblem->typeSolver == FEM_BAND ? "BAND" :
        theProblem->typeSolver == FEM_NONE ? "None" : "UNKNOWN");
    printf("   - Solve time              : %10.4f s\n", (tSolve - tStart)/clk);
    printf("   - Stresses compute time   : %10.4f s\n", (tSolve - tStart)/clk);
    printf("   - Force recovery time     : %10.4f s\n", (tForces - tSolve)/clk);
    printf("   - Integration time        : %10.4f s\n", (tIntegrate - tForces)/clk);
    printf("   - Total FEM pipeline      : %10.4f s\n", (tIntegrate - tStart)/clk);


   
//
//  -4- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 100;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    double *displaySoluce = theSoluce;
    double *unmappedResiduals = theProblem->residuals;
    double *vecU = malloc(theNodes->nNodes * sizeof(double));
    double *vecV = malloc(theNodes->nNodes * sizeof(double));

    if (theProblem->typeSolver == FEM_BAND && theProblem->renumType != FEM_NO) {
        displaySoluce = femElasticityUnmapSoluce(theProblem);
        unmappedResiduals = femElasticityUnmapResiduals(theProblem);
    }
    

    for (int i = 0; i < theNodes->nNodes; i++) {
        int idx = i;
        theNodes->X[i] += displaySoluce[2*idx+0] * deformationFactor;
        theNodes->Y[i] += displaySoluce[2*idx+1] * deformationFactor;
        normDisplacement[i] = sqrt(displaySoluce[2*idx+0]*displaySoluce[2*idx+0] +
                                   displaySoluce[2*idx+1]*displaySoluce[2*idx+1]);
        forcesX[i] = unmappedResiduals[2*i+0];
        forcesY[i] = unmappedResiduals[2*i+1];
        vecU[i] = displaySoluce[2*i + 0];
        vecV[i] = displaySoluce[2*i + 1];
    }
    
    
    

  
    double hMin = femMin(normDisplacement,theNodes->nNodes);  
    double hMax = femMax(normDisplacement,theNodes->nNodes);  
    printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);

//
//  -5- Calcul de la force globaleresultante
//

    double theGlobalForce[2] = {0, 0};
    for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
        theGlobalForce[0] += theForces[2*i+0];
        theGlobalForce[1] += theForces[2*i+1]; }
    printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);

    femElasticityWrite(theProblem,"../data/problem.txt");


//
//  -6- Visualisation du maillage
//  
    
    int mode = 1;
    int fieldSelected = 0;  // 0 = sigmaXX, 1 = sigmaYY, 2 = sigmaXY
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
    glfwMakeContextCurrent(window);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);


    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindowsZoom(theGeometry->theNodes, w, h, zoomFactor);


        t = glfwGetTime();
        if (isDragging) {
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
        
            // Inversion du signe pour que la vue suive la souris
            double dx = xpos - lastX;
            double dy = ypos - lastY;
        
            double panSpeed = 0.1 / zoomFactor;  // Ajuste cette valeur pour la vitesse
        
            panX -= dx * panSpeed;
            panY += dy * panSpeed;
        
            lastX = xpos;
            lastY = ypos;
        }
           
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 2;}
        if (glfwGetKey(window,'Y') == GLFW_PRESS) { mode = 3;}
        if (glfwGetKey(window, '1') == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_KP_1) == GLFW_PRESS) {
            mode = 4; fieldSelected = 0;
        } if (glfwGetKey(window, '2') == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_KP_2) == GLFW_PRESS) {
            mode = 4; fieldSelected = 1;
        } if (glfwGetKey(window, '3') == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_KP_3) == GLFW_PRESS) {
            mode = 4; fieldSelected = 2;
        } if (glfwGetKey(window,'U') == GLFW_PRESS) { mode = 5; }
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        if (t-told > 0.5) {freezingButton = FALSE; }

        //glClearColor(0.0, 0.0, 0.0, 0.0);  // fond noir
        glClearColor(0.96f, 0.94f, 0.88f, 1.0f);  // beige / crème doux
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements,forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 3) {
            glfemPlotField(theGeometry->theElements,forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
            if (mode == 4) {
                double *fieldToPlot = NULL;
                if (fieldSelected == 0) {
                    fieldToPlot = theProblem->sigmaXX;
                    sprintf(theMessage, "Champ : sigmaXX");
                }
                if (fieldSelected == 1) {
                    fieldToPlot = theProblem->sigmaYY;
                    sprintf(theMessage, "Champ : sigmaYY");
                }
                if (fieldSelected == 2) {
                    fieldToPlot = theProblem->sigmaXY;
                    sprintf(theMessage, "Champ : sigmaXY");
                }
                glfemPlotField(theGeometry->theElements, fieldToPlot);
                glfemPlotMesh(theGeometry->theElements);
                glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage);
            }
            if (mode == 5) {
                glfemPlotMesh(theGeometry->theElements); 
                glfemDrawVectors(theGeometry->theElements, vecU, vecV, 50.0); // 50.0 = facteur d'affichage
                sprintf(theMessage, "Vecteurs déplacement (U,V)");
                glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); 
            }
            
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    if (unmappedResiduals != theProblem->residuals) free(unmappedResiduals);
    free(forcesX);
    free(forcesY);
    free(vecU);
    free(vecV);
    femElasticityFree(theProblem) ; 
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}