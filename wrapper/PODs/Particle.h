#pragma once

#include <iostream>

typedef struct Particle {
    double momentum[3];
    int pdg;

    int daughters_size;
    struct Particle** daughters;
    // Track*;
    // Cluster*;

    Particle(int pdg, double momentum[3]): pdg(pdg), daughters_size(0) {
        for (int i = 0; i < 3; ++i) {
            this->momentum[i] = momentum[i]; 
        }
            std::cerr << pdg << std::endl;
    }

    Particle(int pdg, int daughters_size, Particle** daughters) : 
        pdg(pdg), daughters_size(daughters_size), daughters(daughters) {}

} Particle;
