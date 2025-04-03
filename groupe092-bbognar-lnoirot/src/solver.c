#include"slv_fem.h"
#ifndef NORENUMBER 
double *theGlobalCoord;

int compare(const void *nodeOne, const void *nodeTwo) {
    int *iOne = (int *)nodeOne;
    int *iTwo = (int *)nodeTwo;
    double diff = theGlobalCoord[*iOne] - theGlobalCoord[*iTwo];
    if (diff < 0)    return  1;
    if (diff > 0)    return -1;
    return  0;  
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType) {
    int i, *inverse;
    femNodes *theNodes = theMesh->nodes;
    int nNodes = theNodes->nNodes;
    int *number = theNodes->number;
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < nNodes; i++) 
                number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*nNodes);
            for (i = 0; i < nNodes; i++) 
                inverse[i] = i; 
            theGlobalCoord = theNodes->X;
            qsort(inverse, nNodes, sizeof(int), compare);
            for (i = 0; i < nNodes; i++)
                number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*nNodes);
            for (i = 0; i < nNodes; i++) 
                inverse[i] = i; 
            theGlobalCoord = theNodes->Y;
            qsort(inverse, nNodes, sizeof(int), compare);
            for (i = 0; i < nNodes; i++)
                number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option");
    }
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh) {
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    femNodes *theNodes = theMesh->nodes;
    int *number = theNodes->number;

    myBand = 0;
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = number[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin);
        }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin;
    } return(++myBand);
}

#endif
#ifndef NOBANDASSEMBLE

void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc) {
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            if (myCol >= myRow)  myBandSystem->A[myRow][myCol] += Aloc[i*nLoc+j];
        } myBandSystem->B[myRow] += Bloc[i];
    }
}

#endif
#ifndef NOBANDELIMINATE

double  *femBandSystemEliminate(femBandSystem *myBand) {
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Incomplete Cholesky factorization */ 
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot");
        }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor;
        }
    }
        
    /* Back-substitution */
    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i];
    } return(myBand->B);
}

#endif