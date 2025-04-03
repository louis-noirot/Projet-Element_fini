
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ErrorScan(a) femErrorScan(a, __LINE__, __FILE__)
#define Error(a) femError(a, __LINE__, __FILE__)
#define Warning(a) femWarning(a, __LINE__, __FILE__)
#define FALSE 0
#define TRUE 1
#define MAXNAME 256

typedef enum { FEM_TRIANGLE, FEM_QUAD, FEM_EDGE } femElementType;
typedef enum { DIRICHLET_X, DIRICHLET_Y } femBoundaryType;
typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } femElasticCase;

/**
 * @brief Structure containing the nodes of the global mesh
 * @param nNodes Number of nodes
 * @param X x-coordinates of the nodes
 * @param Y y-coordinates of the nodes
 */
typedef struct {
  int nNodes;
  double *X;
  double *Y;
} femNodes;

/**
 * @brief Structure containing the elements of a mesh. Typically there are 2 meshes in a project:
 * one for the edges (segments) and one for the elements inside (e.g. triangles or quads)
 * @param nLocalNode Number of nodes per element, e.g. 2 for a segment (e.g. elements are edges), 3 for a triangle, 4 for a quad
 * @param nElem Number of elements
 * @param elem List of elements' nodes indices (size: `nElem * nLocalNode`).
 * It is arranged as follows: `elem[i*nLocalNode + j]` is the index of the j-th node of the i-th element.
 * For example, for a triangle, `elem[3*i]`, `elem[3*i+1]` and `elem[3*i+2]` are the indices of the 3 nodes of the i-th triangle
 * @param nodes Pointer to the nodes' structure
 */
typedef struct {
  int nLocalNode;
  int nElem;
  int *elem;
  femNodes *nodes;
} femMesh;

/**
 * @brief Structure containing one domain (a boundary) of the general geometry (e.g. the left side of a rectangle)
 * @param mesh Pointer to the mesh structure: contains elements of 2 nodes (edges)
 * @param nElem Number of elements in the domain
 * @param elem List of domain's elements' indices
 * @param name Name of the domain
 */
typedef struct {
  femMesh *mesh;
  int nElem;
  int *elem;
  char name[MAXNAME];
} femDomain;

/**
 * @brief Structure containing the general geometry of the problem
 * @param elementType Type of elements used in the mesh (e.g. triangles or quads)
 * @param theNodes Pointer to the nodes structure containing the nodes of the global mesh
 * @param theElements Pointer to the elements structure
 * @param theEdges Pointer to the edges structure
 * @param nDomains Number of domains in the geometry
 * @param theDomains List of domains
 */
typedef struct {
  femElementType elementType;
  femNodes *theNodes;
  femMesh *theElements;
  femMesh *theEdges;
  int nDomains;
  femDomain **theDomains;
} femGeo;

/**
 * @brief Structure the parent element for a given element type (e.g. triangles or quads or edges)
 * @param n Number of nodes of an element (e.g. 3 for a triangle)
 * @param x2 Function that computes the parent coordinates of the element nodes (if 2D)
 * @param phi2 Function that computes the shape functions at a given integration point (if 2D)
 * @param dphi2dx Function that computes the derivatives of the shape functions at a given integration point (if 2D)
 * @param x Function that computes the parent coordinates of the element nodes (if 1D)
 * @param phi Function that computes the shape functions at a given integration point (if 1D)
 * @param dphidx Function that computes the derivatives of the shape functions at a given integration point (if 1D)
 */
typedef struct {
  int n;
  void (*x2)(double *xsi, double *eta);
  void (*phi2)(double xsi, double eta, double *phi);
  void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
  void (*x)(double *xsi);
  void (*phi)(double xsi, double *phi);
  void (*dphidx)(double xsi, double *dphidxsi);
} femDiscrete;

/**
 * @brief Structure containing the integration rule for a given element type (e.g. triangles or quads)
 * @param n Number of integration points
 * @param xsi Array of size n containing the xsi coordinates of the integration points
 * @param eta Array of size n containing the eta coordinates of the integration points
 * @param weight Array of size n containing the weights of the integration points
 */
typedef struct {
  int n;
  const double *xsi;
  const double *eta;
  const double *weight;
} femIntegration;

/**
 * @brief Structure containing the full system of equations
 * @param B Right-hand side of the system
 * @param A Stiffness matrix of the system
 * @param size Size of the system ( = 2*Nnodes)
 */
typedef struct {
  double *B;
  double **A;
  int size;
} femFullSystem;

typedef struct {
  femDomain *domain;
  femBoundaryType type;
  double value;
} femBoundaryCondition;

typedef struct {
  double E, nu, rho, g;
  double A, B, C;
  int planarStrainStress;
  int nBoundaryConditions;
  femBoundaryCondition **conditions;
  int *constrainedNodes;
  femGeo *geometry;
  femDiscrete *space;
  femIntegration *rule;
  femFullSystem *system;
} femProblem;

void geoInitialize(femGeo *theGeometry);
void geoMeshPrint(femGeo *theGeometry);
femGeo *geoMeshRead(const char *filename);
void femGeoFree(femGeo *theGeometry);
void geoSetDomainName(femGeo *theGeometry, int iDomain, char *name);
int geoGetDomain(femGeo *theGeometry, char *name);
void geoFinalize();

femProblem *femElasticityCreate(femGeo *theGeometry, double E, double nu, double rho, double g, femElasticCase iCase);
void femElasticityFree(femProblem *theProblem);
void femElasticityPrint(femProblem *theProblem);
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
double *femElasticitySolve(femProblem *theProblem);

femIntegration *femIntegrationCreate(int n, femElementType type);
void femIntegrationFree(femIntegration *theRule);

femDiscrete *femDiscreteCreate(int n, femElementType type);
void femDiscreteFree(femDiscrete *mySpace);
void femDiscretePrint(femDiscrete *mySpace);
void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta);
void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi);
void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femFullSystem *femFullSystemCreate(int size);
void femFullSystemFree(femFullSystem *mySystem);
void femFullSystemPrint(femFullSystem *mySystem);
void femFullSystemInit(femFullSystem *mySystem);
void femFullSystemAlloc(femFullSystem *mySystem, int size);
double *femFullSystemEliminate(femFullSystem *mySystem);
void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double value);

double femMin(double *x, int n);
double femMax(double *x, int n);
void femError(char *text, int line, char *file);
void femErrorScan(int test, int line, char *file);
void femWarning(char *text, int line, char *file);

void femSolutionWrite(int nNodes, int nfields, double *data, const char *filename);

#endif
