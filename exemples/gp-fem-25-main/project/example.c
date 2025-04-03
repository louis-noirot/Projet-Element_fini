#include "../benchmark.h"
#include <stdio.h>
#include "fem.h"

void elasticity_solve(const char *meshfile, const char *outfile, double E, double nu, double rho, double g) {
  // Read the mesh and the problem
  femGeo *theGeometry = geoMeshRead(meshfile);
  femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN);

  // Boundary conditions are Dirichlet X and Y
  femElasticityAddBoundaryCondition(theProblem, "Base", DIRICHLET_X, 0.0);
  femElasticityAddBoundaryCondition(theProblem, "Base", DIRICHLET_Y, 0.0);    
  femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_X, 0.0);
  femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_Y, 0.0);    

  // Assemble and solve
  // femElasticityPrint(theProblem);
  double *theSoluce = femElasticitySolve(theProblem);
  
  // Write out the solution
  int nNodes = theGeometry->theNodes->nNodes;
  femSolutionWrite(nNodes, 2, theSoluce, outfile);

  // free the allocated ressources
  femElasticityFree(theProblem);
  femGeoFree(theGeometry);
}
