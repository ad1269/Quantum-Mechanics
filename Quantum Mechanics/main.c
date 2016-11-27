//
//  main.c
//  Quantum Mechanics
//
//  Created by Mohan Lakshmanan on 7/22/15.
//  Copyright (c) 2015 AD. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

//Wavefunction Struct
typedef struct Wavefunction {
    double f[100000];
    double g[100000];
} Wavefunction;

//Constants
const int Nx = 16000;
const int Nt = 400000;
const double dx = 0.01;
const double dt = 0.0001;
const double k0 = 4.0;
const double A = 1.0;
const double hbar = 1.0;
const double m = 2.0;
const double pi = M_PI;

const char* path = "/Users/mohanraj/Dropbox/Github Projects/On Resume/Quantum Mechanics";

//Global Count Variables
int fileCount = 1;
int threadCount = 1;

//Method Declarations
void initialize_WF(double f[], double g[]);
void initialize_PE(double v[]);
void PE(double f[], double g[], double v[]);
void KE(double f[], double g[], int start);
void verifyNorm(double f[], double g[]);
void evolve_WF(double f[], double g[], double v[]);

//Thread creating methods
void printThread(double f[], double g[]);

//Subthread methods
void *print_WF(void *args);

int main() {
    //Declare all variables
    double f[100000], g[100000];
    double v[100000];
    
    //Initialize Wavefunction/Potential and print it out
    initialize_WF(f, g);
    printThread(f, g);
    initialize_PE(v);
    
    //Run main time loop
    evolve_WF(f, g, v);
    
    //Print calculated wavefunction
    printThread(f, g);
    
    return 0;
}

void initialize_WF(double f[], double g[]) {
    int ix;
    double norm, x, x0;
    x0 = Nx * dx/2.;
    norm = 0.;
    
    //Initialize the Schrodinger Wavefunction
    for(ix = 0; ix < Nx; ix++) {
        x = dx * ix;
        f[ix] = exp(-(x-x0)*(x-x0)/(4.0*A)) * cos(k0*(x-x0));
        g[ix] = exp(-(x-x0)*(x-x0)/(4.0*A)) * sin(k0*(x-x0));
        
        norm += f[ix] * f[ix] + g[ix] * g[ix];
    }
    if(fabs(norm - 1.) <= 0.000001) {
        printf("Failure!\n%lf\n", norm);
    }
    
    //Normalize the wavefunction
    norm = sqrt(norm);
    for(ix = 0; ix < Nx; ix++) {
        f[ix] = f[ix] / norm;
        g[ix] = g[ix] / norm;
    }
}

void initialize_PE(double v[]) {
    int ix;
    //Initialize potential energy to 0
    for(ix = 0; ix < Nx; ix++) {
        v[ix] = 0;
    }
    
    //Add barriers here
    for(ix = 12000; ix < 16000; ix++) {
        v[ix] = 100000000;
    }
    for(ix = 0; ix < 4000; ix++) {
        v[ix] = 100000000;
    }
}

void PE(double f[], double g[], double v[]) {
    int ix;
    double phi, fp, gp;
    double B = hbar * hbar / (2.*m*dx*dx);
    
    //Potential Energy
    for(ix = 0; ix < Nx; ix++) {
        phi = (2.0*B-v[ix]) * dt;
        fp = cos(phi) * f[ix] + sin(phi)*g[ix];
        gp = -sin(phi) * f[ix] + cos(phi)*g[ix];
        
        f[ix] = fp;
        g[ix] = gp;
    }

}

void KE(double f[], double g[], int start) {
    int ix;
    double phi, fp, gp, fpp, gpp;
    double B = hbar * hbar / (2.*m*dx*dx);

    //Kinetic Energy
    for(ix = start; ix < Nx - 1; ix += 2) {
        phi = B * dt;
        fp = cos(phi) * f[ix] - sin(phi) * g[ix+1];
        gp = sin(phi) * f[ix+1] + cos(phi) * g[ix];
        fpp = cos(phi) * f[ix + 1] - sin(phi) * g[ix];
        gpp = sin(phi) * f[ix] + cos(phi) * g[ix+1];
        
        f[ix] = fp;
        g[ix] = gp;
        f[ix+1] = fpp;
        g[ix+1] = gpp;
    }
}

void printThread(double f[], double g[]) {
    //Increase thread count
    threadCount++;
    
    //Declare threads for multitasking
    pthread_t printThread;
    Wavefunction *args = malloc(sizeof *args);
    int ix;

    for (ix = 0; ix < Nx; ix++) {
        args->f[ix] = f[ix];
        args->g[ix] = g[ix];
    }
    
    //Creating thread
    pthread_create(&printThread, NULL, print_WF, args);
}

void *print_WF(void *args) {
    Wavefunction *wf = args;
    int ix;
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    
    char name[1024];
    sprintf(name, "%s/Quantum Mechanics/Well/schrodinger%i", path, fileCount);
    FILE *fileout = fopen(name, "w");
    
    //Print out square of wavefunction
    for(ix = 0; ix < Nx; ix++) {
        double mod = wf->f[ix] * wf->f[ix] + wf->g[ix] * wf->g[ix];
        if(ix % 10 == 0) {
            fprintf(fileout, "%8i %12.6lf\n", ix, mod);
        }
    }
    
    fclose(fileout);
    fileCount++;
    threadCount--;
    return NULL;
}

void evolve_WF(double f[], double g[], double v[]) {
    //Main time loop
    int it;
    for(it = 1; it < Nt; it++) {
        
        //Print the wavefunction every 5000 runs
        if(it % 5000 == 0) {
            printThread(f, g);
        }

        PE(f, g, v);
        KE(f, g, 0);
        KE(f, g, 1);
        
        //Verify no calculation errors have been made
        //verifyNorm(f, g);
    }
}

void verifyNorm(double f[], double g[]) {
    double norm = 0;
    int ix;
    //Calculate norm to confirm the calculation has been done properly
    for(ix = 0; ix < Nx; ix++) {
        norm += f[ix] * f[ix] + g[ix] * g[ix];
    }
    if(fabs(norm - 1.) <= 0.000001) {
        printf("Failure!\n%lf\n", norm);
    }
}
