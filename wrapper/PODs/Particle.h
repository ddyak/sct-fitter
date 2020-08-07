#pragma once

#include <iostream>

namespace POD {

typedef struct Particle {
    int pdg;
    int daughters_size;
    struct Particle** daughters;
    
    // dummy particle
    double momentum[3];
    double covariance[3][3];

    // Track*;
    // Cluster*;

    Particle(int pdg, double momentum[3], double covariance[3][3]):
        pdg(pdg), daughters_size(0) 
    {
        for (int i = 0; i < 3; ++i) {
            this->momentum[i] = momentum[i];
            for (int j = 0; j < 3; ++j) {
                this->covariance[i][j] = covariance[i][j];
            }
        }
    }

    Particle(int pdg, int daughters_size, Particle** daughters):  
        pdg(pdg), daughters_size(daughters_size) 
    {
        this->daughters = new Particle*[daughters_size];
        
        for (int i = 0; i < daughters_size; ++i) {
            this->daughters[i] = daughters[i];
        } 
        
        for (int daughter_idx = 0; daughter_idx < daughters_size; ++daughter_idx) {
            for (int i = 0; i < 3; ++i) {
                this->momentum[i] += daughters[daughter_idx]->momentum[i]; 
            }
        }
    }

} Particle;

}