#include "../../Project/src/fem.h"
#include "gmsh.h"

// GMSH Functions:

int createRectangle(double x, double y, double width, double height){
    int ierr;
    return gmshModelOccAddRectangle(x, y, 0.0, width, height, -1, 0, &ierr);
}

void cutElement(int *mainElement, int *cutElement) {
    int ierr;
    gmshModelOccCut(mainElement, 2, cutElement, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
}

int createDisk(double xc, double yc, double rx, double ry) {
    int ierr;
    return gmshModelOccAddDisk(xc, yc, 0.0, rx, ry, -1, NULL, 0, NULL, 0, &ierr);
}

// Position Geometry: /

// Materials:

double *getMaterialProperties(char *material){
    double *properties = malloc(3 * sizeof(double));
    if (properties == NULL) { 
        Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return NULL;
    }
    if (strcmp(material, "steel") == 0)
    {
        properties[0] = 200.0e9; // E [Pa]
        properties[1] = 0.3;     // nu [/]
        properties[2] = 7850.0;  // rho [kg/m^3]
    }
    else if (strcmp(material, "reinforced_concrete") == 0)
    {
        properties[0] = 30.0e9; // E [Pa]
        properties[1] = 0.2;    // nu [/]
        properties[2] = 2400.0; // rho [kg/m^3]
    }
    else {
        Error("Material Unknown.");
    } return properties;
}

// Geosize:
double geoSize(double x, double y) {
    femGeometry *theGeometry = geoGetGeometry();
    double h = theGeometry->LxPlate;
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

// Generate Mesh:
void geoMeshGenerate(femDiscreteType discreteType, int bridgeSimplified) {
    femGeometry *theGeometry = geoGetGeometry();
    
    // Define the geometry parameters
    theGeometry->LxPlate = 20.0;
    theGeometry->LyPlate = 55.0;

    theGeometry->geoSize = geoSize;
    theGeometry->getMaterialProperties = getMaterialProperties;

    int ierr;

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    double b = h*(2.0/55.0);   // beam
    double c = w/10.0;         // column
    double h1 = (14.0/55.0)*h; // height of the bottom rectangle
    double h2 = (18.0/55.0)*h; // height of the middle rectangle
    double h3 = (17.0/55.0)*h; // height of the top rectangle
    int idMain = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr);
    int idBotom = gmshModelOccAddRectangle(c,0.0,0.0,w-2*c,h1,-1,0.0,&ierr);
    int idMiddle = gmshModelOccAddRectangle(c,b+h1,0.0,w-2*c,h2,-1,0.0,&ierr);
    int idTop = gmshModelOccAddRectangle(c,b*2+h1+h2,0.0,w-2*c,h3,-1,0.0,&ierr);
    int main[] = {2,idMain};
    int bottom[] = {2,idBotom};
    int middle[] = {2,idMiddle};
    int top[] = {2,idTop};

    gmshModelOccCut(main,2,bottom,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(main,2,middle,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(main,2,top,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);
    
    // Generate triangles meshing
    if (theGeometry->elementType == FEM_TRIANGLE) {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { 
            gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr);
        }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    // Generate quads meshing
    else if (theGeometry->elementType == FEM_QUAD) {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC){
            gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr);
        }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    gmshFltkInitialize(&ierr);
}

// Domain Name:

void setDomainsName(int bridgeSimplified) {
    typedef struct domainMapping {
        int id;
        char *name;
    } domain_Mapping_t;

    domain_Mapping_t domain_mapping[] = {
        {0,"2_beam_T"},
        {1,"R_column_3"},
        {2,"3_beam_B"},
        {3,"L_column_3"},
        {4,"L_column"},
        {5,"3_beam_T"},
        {6,"Bottom1"},
        {7,"R_column"},
        {8,"L_column_1"},
        {9,"Bottom2"},
        {10,"1_beam_B"},
        {11,"R_column_1"},
        {12,"1_beam_T"},
        {13,"L_column_2"},
        {14,"R_column_2"},
        {15,"2_beam_B"},
    };

    const int NB_DOMAINS = sizeof(domain_mapping) / sizeof(domain_mapping[0]);
    for (int i = 0; i < NB_DOMAINS; ++i) {
        domain_Mapping_t domain = domain_mapping[i];
        geoSetDomainName(domain.id, domain.name);
    }
}

void createBoundaryConditions(femProblem *theProblem, int bridgeSimplified) {
    typedef struct domainBoundaryMapping {
        char *name;
        femBoundaryType type;
        double value1;
        double value2;
    } domainBoundaryMapping_t;

    // Define constants
    double force_bateau = 1e15; // TODOOOOOO
    const double width_bateau = 2.0;

    domainBoundaryMapping_t mapping[] = {

            {"Bottom1", NEUMANN_X, 0.0, NAN},
            {"Bottom1", NEUMANN_Y, 0.0, NAN},
            {"Bottom2", NEUMANN_X, 0.0, NAN},
            {"Bottom2", NEUMANN_Y, 0.0, NAN},
            {"R_column_1", NEUMANN_Y, 1e10000, NAN},
    };

    const int NB_DOMAINS = sizeof(mapping) / sizeof(mapping[0]);

    for (int i = 0; i < NB_DOMAINS; ++i) {
        domainBoundaryMapping_t domain = mapping[i];
        femElasticityAddBoundaryCondition(theProblem, domain.name, domain.type, domain.value1, domain.value2);
    }
}

/***************************************************/
/********* MESH + GEO SIZE (ON AN EXAMPLE) *********/
/***************************************************/

double geoSizeExample(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    return theGeometry->defaultSize;
}

void geoMeshGenerateExample(femDiscreteType discreteType, int beam_example)
{
    if (beam_example == TRUE) { geoMeshGenerate_Beam(discreteType); }
    else                      { geoMeshGenerate_UForm(discreteType); }
}

void geoMeshGenerate_UForm(femDiscreteType discreteType)
{
    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
                \        |
                )        |
                /        |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    femGeometry *theGeometry = geoGetGeometry();
    double Lx = 1.0;
    double Ly = 1.0;
    theGeometry->LxPlate = Lx;
    theGeometry->LyPlate = Ly;
    theGeometry->defaultSize = Lx * 0.02;
    theGeometry->geoSize = geoSizeExample;

    geoSetSizeCallback(theGeometry->geoSize);

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;

    int ierr;
    double r = w / 4;
    int idRect = createRectangle(0.0, 0.0, w, h);
    int idDisk = createDisk(w / 2.0, h / 2.0, r, r);
    int idSlit = createRectangle(w / 2.0, h / 2.0 - r, w, 2.0 * r);
    int rect[] = {2, idRect};
    int disk[] = {2, idDisk};
    int slit[] = {2, idSlit};

    cutElement(rect, disk);
    cutElement(rect, slit);

    gmshModelOccSynchronize(&ierr);

    if (theGeometry->elementType == FEM_QUAD)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }

    if (theGeometry->elementType == FEM_TRIANGLE)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    
    return;
}

void geoMeshGenerate_Beam(femDiscreteType discreteType)
{
    femGeometry *theGeometry = geoGetGeometry();
    double Lx = 8.0;
    double Ly = 1.0;
    theGeometry->LxPlate = Lx;
    theGeometry->LyPlate = Ly;
    theGeometry->defaultSize = Lx * 0.01;
    theGeometry->geoSize = geoSizeExample;

    geoSetSizeCallback(theGeometry->geoSize);

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;

    int ierr;
    double r = w / 4;
    int idRect = createRectangle(0.0, 0.0, w, h);

    gmshModelOccSynchronize(&ierr);

    if (theGeometry->elementType == FEM_QUAD)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }

    if (theGeometry->elementType == FEM_TRIANGLE)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    
    return;
}