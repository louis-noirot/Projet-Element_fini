
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gmshc.h"


#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define ErrorGmsh(a)   femErrorGmsh(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define MAXNAME 256

typedef enum {FEM_TRIANGLE,FEM_QUAD,FEM_EDGE} femElementType;
typedef enum {DIRICHLET_X,DIRICHLET_Y,NEUMANN_X,NEUMANN_Y} femBoundaryType;
typedef enum {PLANAR_STRESS,PLANAR_STRAIN,AXISYM} femElasticCase;
typedef enum {FEM_FULL, FEM_BAND, FEM_NONE} femSolverType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

typedef struct {
    int nNodes;
    double *X;
    double *Y;
    int *number;
} femNodes;

typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
} femMesh;

typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[MAXNAME];
} femDomain;

typedef struct {
    double LxPlate, LyPlate;
    double h;
    femElementType elementType;
    double (*geoSize)(double x, double y);
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeo;

typedef struct {
    int n;
    femElementType type;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x)(double *xsi);
    void (*phi)(double xsi, double *phi);
    void (*dphidx)(double xsi, double *dphidxsi);
} femDiscrete;
    
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

typedef struct {
    int size;
    double *B;
    double *X;
    double **A;
} femNoneSystem;

typedef struct {
    double *B;
    double *X;
    double **A;        
    int size;
    int band;        
} femBandSystem;


typedef struct {
    femDomain* domain;
    femBoundaryType type; 
    double value;
} femBoundaryCondition;

typedef struct {
    double E,nu,rho,g;
    double A,B,C;
    int planarStrainStress;
    int nBoundaryConditions;
    femBoundaryCondition **conditions;  
    int *constrainedNodes; 
    double *soluce;
    double *residuals;
    femGeo *geometry;
    femDiscrete *space;
    femIntegration *rule;
    femDiscrete *spaceEdge;
    femIntegration *ruleEdge;
    femFullSystem *system;       // Solveur complet
    femSolverType typeSolver;    // Choix du solveur
    void *systemFull;            // Solveur full
    void *systemNone;           // Solveur None
    void *systemBand;            // Solveur bande
    femRenumType renumType;      // Type de renumerotation
    double *sigmaXX;
    double *sigmaYY;
    double *sigmaXY;
} femProblem;



void                geoInitialize();
femGeo*             geoGetGeometry();
double              geoSize(double x, double y);
double              geoSizeDefault(double x, double y);
void                geoSetSizeCallback(double (*geoSize)(double x, double y));
void                geoMeshGenerate();
void                geoMeshImport();
void                geoMeshPrint();
void                geoMeshWrite(const char *filename);
void                geoMeshRead(const char *filename);
void                geoSetDomainName(int iDomain, char *name);
int                 geoGetDomain(char *name);
void                geoFinalize();

femProblem*         femElasticityCreate(femGeo* theGeometry, double E, double nu, double rho, double g, femElasticCase iCase, femSolverType typeSolver, femRenumType renumType);
void                femElasticityFree(femProblem *theProblem);
void                femElasticityPrint(femProblem *theProblem);
void                femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
void                femElasticityAssembleElements(femProblem *theProblem);
void                femElasticityAssembleNeumann(femProblem *theProblem);
double*             femElasticitySolve(femProblem *theProblem);
double*             femElasticityForces(femProblem *theProblem);
double              femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y));


femIntegration*     femIntegrationCreate(int n, femElementType type);
void                femIntegrationFree(femIntegration *theRule);

femDiscrete*        femDiscreteCreate(int n, femElementType type);
void                femDiscreteFree(femDiscrete* mySpace);
void                femDiscretePrint(femDiscrete* mySpace);
void                femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void                femDiscreteXsi(femDiscrete* mySpace, double *xsi);
void                femDiscretePhi(femDiscrete* mySpace, double xsi, double *phi);
void                femDiscreteDphi(femDiscrete* mySpace, double xsi, double *dphidxsi);

femFullSystem*      femFullSystemCreate(int size);
void                femFullSystemFree(femFullSystem* mySystem);
void                femFullSystemPrint(femFullSystem* mySystem);
void                femFullSystemInit(femFullSystem* mySystem);
void                femFullSystemAlloc(femFullSystem* mySystem, int size);
double*             femFullSystemEliminate(femFullSystem* mySystem);
void                femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);

double              femMin(double *x, int n);
double              femMax(double *x, int n);
void                femError(char *text, int line, char *file);
void                femErrorScan(int test, int line, char *file);
void                femErrorGmsh(int test, int line, char *file);
void                femWarning(char *text, int line, char *file);
void                femElasticityWrite(femProblem *theProblem, const char *filename);


// None solver
femNoneSystem*     femNoneSystemCreate(int size);
void                femNoneSystemFree(femNoneSystem* mySystem);
double*             femNoneSystemEliminate(femNoneSystem* mySystem);
void                femSystemAddRightHandSide(femProblem *theProblem, int i, double value);
void                femSystemAdd(femProblem *theProblem, int i, int j, double value);
void                femSystemConstraint(femProblem *theProblem, int i, double value);


// Band solver
int                 compare(const void *nodeOne, const void *nodeTwo);
void                femMeshRenumber(femMesh *theMesh, femRenumType renumType);
int                 femMeshComputeBand(femMesh *theMesh);
void                femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
double*             femBandSystemEliminate(femBandSystem *myBand);
femBandSystem*      femBandSystemCreate(int size, int band);
void                femBandSystemInit(femBandSystem *myBandSystem);
void                femBandSystemConstrain(femBandSystem *myBand, int myNode, double value);
double*             femElasticityUnmapSoluce(femProblem *theProblem);
double*             femElasticityUnmapResiduals(femProblem *theProblem);
void                femElasticityComputeStresses(femProblem *theProblem);



#endif
