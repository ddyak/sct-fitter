#include <iostream>
#include <stdio.h>
#include <map>

#include "dataobjects/Particle.h"
#include "PODs/Particle.h"

#include "treefitter/include/ConstraintConfiguration.h"
#include "treefitter/include/FitManager.h"


sct::ana::ParticlePtr createParticle(Particle* part, std::map<sct::ana::ParticlePtr, Particle*>& particle2pod) {
    if (part->daughters_size > 0) {
        std::vector<sct::ana::ParticlePtr> daughters;
        for (int i = 0; i < part->daughters_size; ++i) {
            auto daughter = createParticle(part->daughters[i], particle2pod);
            daughters.push_back(daughter);
        }
        auto particle = std::make_shared<sct::ana::Particle>(daughters, part->pdg);
       particle2pod[particle] = part; // for pod
        return particle;
    } else {
        auto momentum = sct::kine::ThreeVector<double>{part->momentum[0], part->momentum[1], part->momentum[2]};
        auto particle = std::make_shared<sct::ana::Particle>(momentum, 0.139, part->pdg);
       particle2pod[particle] = part; // for pod
        return particle;
    }
}

extern "C" {

Particle* Particle_from_daughters(int pdg, int daughters_size, Particle** daughters) { 
    return new Particle(pdg, daughters_size, daughters);
}

Particle* Particle_from_momentum(int pdg, double momentum[3]) { 
    return new Particle(pdg, momentum); 
}

double fit(Particle* part) {
    std::map<sct::ana::ParticlePtr, Particle*> particle2pod;
    sct::ana::ParticlePtr particle = createParticle(part, particle2pod);
    sct::ana::FitManager fitmanager(particle, {});
    fitmanager.fit();
    for (auto& [particle, pod]: particle2pod) {
        pod->momentum[0] = particle->momentum()[0];
        pod->momentum[1] = particle->momentum()[1];
        pod->momentum[2] = particle->momentum()[2];
    }
    return fitmanager.chiSquare();
}

}
