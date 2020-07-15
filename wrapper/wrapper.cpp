#include <iostream>
#include <stdio.h>
#include <map>

#include "dataobjects/Particle.h"
#include "PODs/Particle.h"

std::map<Particle*, sct::ana::ParticlePtr> pod2particle;
std::map<sct::ana::ParticlePtr, Particle*> particle2pod;

sct::ana::ParticlePtr createParticle(Particle* part) {
    part->pdg = 7; // temp for testing
    if (part->daughters_size) {
        std::vector<sct::ana::ParticlePtr> daughters;
        for (int i = 0; i < part->daughters_size; ++i) {
            daughters.push_back(createParticle(part->daughters[i]));
        }
        return std::make_shared<sct::ana::Particle>(daughters, part->pdg);
    } else {
        auto momentum = sct::kine::ThreeVector<double>{part->momentum[0], part->momentum[1], part->momentum[2]};
        return std::make_shared<sct::ana::Particle>(momentum, 0.139, part->pdg);
    }
}

extern "C" {

Particle* Particle_from_daughters(int pdg, int daughters_size, Particle** daughters) { 
    return new Particle(pdg, daughters_size, daughters);
}

Particle* Particle_from_momentum(int pdg, double momentum[3]) { 
    return new Particle(pdg, momentum); 
}

void fit(Particle* part) {
    createParticle(part);
    part->pdg = 7;
}

}
