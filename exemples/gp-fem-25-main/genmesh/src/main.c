/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"

int main(int argc, char *argv[]) {

  // ================ Modify to your needs ================

  // Global mesh size factor (i.e. lower = smaller triangles)
  // NOTE : don't put it too low, 0.01 will already generate more than 10^5 tri
  double clscale = 1.0;

  // Do you want the mesh shown graphically after generation ?
  // 0 = no
  // 1 = yes
  int plot = 1;

  // =================== DON'T TOUCH (you can read though) ===================
  femBoundaryType arctype = UNDEFINED;
  double arcvalue0 = 0.0;
  double arcvalue1 = 0.0;
  femBoundaryType toptype = UNDEFINED;
  double topvalue0 = 0.0;
  double topvalue1 = 0.0;
  int axisym = 0;

  for (int i = 0; i < argc; ++i) {
    if (strncasecmp(argv[i], "-axisym", 8) == 0)
      axisym = 1;
    if (strncasecmp(argv[i], "-arc", 5) == 0) {
      arctype = atoi(argv[i + 1]);
      arcvalue0 = atof(argv[i + 2]);
      arcvalue1 = atof(argv[i + 3]);
    }
    if (strncasecmp(argv[i], "-top", 5) == 0) {
      toptype = atoi(argv[i + 1]);
      topvalue0 = atof(argv[i + 2]);
      topvalue1 = atof(argv[i + 3]);
    }
    if (strncasecmp(argv[i], "-clscale", 9) == 0) {
      clscale = atof(argv[i + 1]);
    }
  }

  //
  //  -1- Construction de la geometrie
  //
  geoInitialize();
  femGeo *theGeometry = geoGetGeometry();

  // OPTION 1 : Utilisation de GMSH avec OpenCascade
  // theGeometry->h = 0.5;
  geoMeshGenerate(0.1 * clscale);

  // OPTION 2 : Utilisation de GMSH directement
  // theGeometry->h = 0.05;
  // geoMeshGenerateGeo();

  // OPTION 3 : Lecture d'un fichier .geo
  // theGeometry->h = 0.05;
  // geoMeshGenerateGeoFile("../../data/mesh.geo");

  // OPTION 4 : Lecture d'un fichier .msh
  // geoMeshGenerateMshFile("../../data/mesh.msh");

  geoMeshImport();
  geoSetDomainName(0, "symmetry");
  geoSetDomainName(1, "top");
  geoSetDomainName(4, "arc");
  geoSetDomainName(7, "base");
  geoMeshWrite("../../data/mesh.txt");

  //
  //  -2- Champ de la taille de référence du maillage (uniquement pour la visualisation)
  //

  double *meshSizeField = malloc(theGeometry->theNodes->nNodes * sizeof(double));
  femNodes *theNodes = theGeometry->theNodes;
  for (int i = 0; i < theNodes->nNodes; ++i)
    meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
  double hMin = femMin(meshSizeField, theNodes->nNodes);
  double hMax = femMax(meshSizeField, theNodes->nNodes);
  printf(" ==== Global requested h : %14.7e \n", theGeometry->h);
  printf(" ==== Minimum h          : %14.7e \n", hMin);
  printf(" ==== Maximum h          : %14.7e \n", hMax);
  printf(" ==== Mesh has %d elements and %d nodes\n", theGeometry->theElements->nElem, theGeometry->theNodes->nNodes);

  //
  //  -3- Visualisation
  //

  if (plot) {
    int mode = 1;
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];

    GLFWwindow *window = glfemInit("EPL1110 : Project 2022-23 ");
    glfwMakeContextCurrent(window);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    do {
      int w, h;
      glfwGetFramebufferSize(window, &w, &h);
      glfemReshapeWindows(window, theGeometry->theNodes, w, h);

      t = glfwGetTime();
      if (glfwGetKey(window, 'D') == GLFW_PRESS) {
        mode = 0;
      }
      if (glfwGetKey(window, 'V') == GLFW_PRESS) {
        mode = 1;
      }
      if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) {
        domain++;
        freezingButton = TRUE;
        told = t;
      }

      if (t - told > 0.5) {
        freezingButton = FALSE;
      }
      if (mode == 1) {
        glfemPlotField(theGeometry->theElements, meshSizeField);
        glfemPlotMesh(theGeometry->theElements);
        sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
        glColor3f(1.0, 0.0, 0.0);
        glfemMessage(theMessage);
      }
      if (mode == 0) {
        domain = domain % theGeometry->nDomains;
        glfemPlotDomain(theGeometry->theDomains[domain]);
        sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
        glColor3f(1.0, 0.0, 0.0);
        glfemMessage(theMessage);
      }

      glfwSwapBuffers(window);
      glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);

    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
  }
  return 0;
}
