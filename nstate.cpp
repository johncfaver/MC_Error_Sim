#include <cmath>
#include <iostream>
#include <stdlib.h>

/*
    This simulation computes the free energy and estimated uncertainty
    in free energy at increasing ensemble sizes. Each microstate sampled
    is assumed to have uniform systematic and random errors in energy,
    defined by errsys and errrand. The 1D potential energy surface is modeled 
    by eFunction, which currently is a Lennard-Jones surface.
    
    For each size ensemble, the free energy is calculated, along with the overall
    propagated systematic and random errors in free energy. Systematic error
    magnitudes tend to increase and converge to some value, while random error
    magnitudes tend to decrease toward zero.

*/

/* ADJUSTABLE PARAMETERS */
int   nstates   =   500;                 // Maximum ensemble size
int   nsteps    =   10000;               // Number of MC steps
float errsys    =   0.0;                 // Systematic error per microstate
float errrand   =   1.0;                 // Random error per microstate
float beta      =   1.6785192002458318;  // (RT)^-1 in (kcal/mol)^-1 at 300K
/* ADJUSTABLE PARAMETERS */

float *energies = new float[nstates]; // Array for holding microstate energies

//Box-Muller algorithm for generating normally distributed random numbers.
float gauNoise(float u,float stddev){
    float s=0.0;
    float u1=0;
    float u2=0;
    while(s<=0.0 or s>=1.0){
        u1 =2.0*(float)rand()/(float)RAND_MAX-1.0;
        u2 =2.0*(float)rand()/(float)RAND_MAX-1.0;
        s = u1*u1 + u2*u2;
    }
    return sqrt(-2.0*log(s)/s)*u1*stddev+u;
}

// LJ-type potential energy function
float eFunction(float x){ 
       return 5.0*(pow(x,-12.0)-pow(x,-6.0));
}

int main(int nargs,char* args[]){
    
    srand((unsigned)time(0));

    double freeenergy=0.0;                
    double avgfreeenergy=0.0;
    double stdfreeenergy=0.0;
    double realfreeenergy=0.0;
    double M=0.0;
    double davg=0.0;
    double eguess=0.0;
    
    std::cout << "#N    RealA    AvgA    TotSysErr    TotRandErr\n";

/* 
    We now loop through different ensemble sizes t, from 1 to nstates.
    For each ensemble, we generate t microstate energies, which are selected
        randomly from the PES. We restrict our sampling to the region between
        0.9 and 4.9 (attractive region of the potential), and only allow 
        bound states with negative energies.
    We then do an MC error propagation for this ensemble of t microstates.
        For each MC step s, we estimate free energy by perturbing the 
        microstate energies by randomly drawing from a normal distribution
        with mean errsys and standard deviation errrand. We finally report 
        the mean and standard deviation of free energy estimations for this
        ensemble of size t.
*/

    for(int t=1;t<=nstates;t++){    
        realfreeenergy=0.0;
        for(int j=0;j<t;j++){
            while(eguess>=0.0){
                eguess=eFunction((float)rand()/(float)RAND_MAX*4.0+0.9);    
            }
            energies[j]=eguess;
            realfreeenergy+=exp(-beta*energies[j]);
        }
        realfreeenergy=log(realfreeenergy)/-beta;
        
        avgfreeenergy=0.0;
        M=0.0;
        for(int s=0;s<nsteps;s++){
            freeenergy=0.0;
            for(int n=0;n<t;n++){
                freeenergy+=exp(-beta*(energies[n]+gauNoise(errsys,errrand)));
            }
            freeenergy=log(freeenergy)/-beta;
            davg=freeenergy-avgfreeenergy;
            avgfreeenergy+=davg/(s+1);
            M+=davg*(freeenergy-avgfreeenergy);
        }
        stdfreeenergy=sqrt(M/(nsteps-1));
    
        std::cout << t << " " <<  realfreeenergy << " " << avgfreeenergy << " ";
        std::cout << (avgfreeenergy-realfreeenergy) << " " << stdfreeenergy << "\n";

    }
    return 0;
}
