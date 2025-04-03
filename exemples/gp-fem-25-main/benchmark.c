#include "benchmark.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef WARMUP
#define WARMUP 1
#endif

#ifndef NRUNS
#define NRUNS 5
#endif

int main(int argc, char** argv){

    // NOTE: The default input mesh and output file are 
    // "../data/mesh.txt" and "../data/UV.txt" but you can 
    // change them either by modifying the lines below 
    // or by passing the argument in command line :
    //    $ ./benchmark <mesh-file> <output-file>
    char* meshfile = "../data/mesh.txt";
    char* outfile  = "../data/UV.txt";
    if (argc > 1) meshfile = argv[1];
    if (argc > 2) outfile  = argv[2];

    struct timespec t0, t1;
    double* const restrict times = malloc(sizeof(*times)*NRUNS);
    double E   = 211.e9;
    double nu  = 0.3;
    double rho = 7.85e3; 
    double g   = 9.81;

    // ====== READY, SET... ========
    for (int i = 0; i < WARMUP; i++){
        elasticity_solve(meshfile, outfile, E, nu, rho, g);
    }
    printf("Warmup runs done, now benchmarking\n");

    // =========  GO ! =============
    for (int i = 0; i < NRUNS; i++){
        timespec_get(&t0, TIME_UTC);
        elasticity_solve(meshfile, outfile, E, nu, rho, g);
        timespec_get(&t1, TIME_UTC);

        // Do not print the value here to not incur IO overhead
        times[i] = (t1.tv_sec - t0.tv_sec)*1.0 + (t1.tv_nsec - t0.tv_nsec)*1e-9;
    }

    double tsum = 0;
    double tsum2 = 0;
    for (int i = 0; i < NRUNS; i++){
        tsum += times[i];
        tsum2 += times[i] * times[i];
    }
    double tmean = tsum / NRUNS;
    double tvar = (tsum2 - NRUNS*tmean*tmean) / (NRUNS-1);
    double tstd = sqrt(tvar);

    char prefix[2] = {0};
    if (tmean < 1e-3){
        prefix[0] = 'u';
        tmean *= 1e6;
        tstd *= 1e6;
    } else if (tmean < 1.0){
        prefix[0] = 'm';
        tmean *= 1e3;
        tstd *= 1e3;
    }
    printf("Your code runs in %.4f ± %.4f %ss for mesh file '%s'\n", tmean, tstd, prefix, meshfile);

    free(times);
    return 0;
}
