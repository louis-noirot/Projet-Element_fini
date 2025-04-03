#include "fem.h"

double *femElasticitySolve(femProblem *theProblem) {

  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;

  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int map[4], mapX[4], mapY[4];

  int nLocal = theMesh->nLocalNode;

  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double g = theProblem->g;
  double **A = theSystem->A;
  double *B = theSystem->B;

  //
  //  A faire :-)
  //
  int *theConstrainedNodes = theProblem->constrainedNodes;
  for (int i = 0; i < theSystem->size; i++) {
    if (theConstrainedNodes[i] != -1) {
      double value = theProblem->conditions[theConstrainedNodes[i]]->value;
      femFullSystemConstrain(theSystem, i, value);
    }
  }

  return femFullSystemEliminate(theSystem);
}


double *femFullSystemEliminate(femFullSystem *mySystem) {
  double **A, *B, factor;
  int i, j, k, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  /* Gauss elimination */

  for (k = 0; k < size; k++) {
    if (fabs(A[k][k]) <= 1e-16) {
      printf("Pivot index %d  ", k);
      printf("Pivot value %e  ", A[k][k]);
      Error("Cannot eliminate with such a pivot");
    }
    for (i = k + 1; i < size; i++) {
      factor = A[i][k] / A[k][k];
      for (j = k + 1; j < size; j++)
        A[i][j] = A[i][j] - A[k][j] * factor;
      B[i] = B[i] - B[k] * factor;
    }
  }

  /* Back-substitution */

  for (i = size - 1; i >= 0; i--) {
    factor = 0;
    for (j = i + 1; j < size; j++)
      factor += A[i][j] * B[j];
    B[i] = (B[i] - factor) / A[i][i];
  }

  return (mySystem->B);
}

void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue) {
  double **A, *B;
  int i, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  for (i = 0; i < size; i++) {
    B[i] -= myValue * A[i][myNode];
    A[i][myNode] = 0;
  }

  for (i = 0; i < size; i++)
    A[myNode][i] = 0;

  A[myNode][myNode] = 1;
  B[myNode] = myValue;
}
