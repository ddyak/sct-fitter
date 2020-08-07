#include <iostream>
#include <stdio.h>
#include <map>

#include "dataobjects/Particle.h"
#include "PODs/Particle.h"

#include "treefitter/include/ConstraintConfiguration.h"
#include "treefitter/include/FitManager.h"


sct::ana::ParticlePtr createParticle(POD::Particle* part, std::map<sct::ana::ParticlePtr, POD::Particle*>& particle2pod) {
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
        auto covariance = Eigen::Matrix<double, 3, 3>{};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                covariance(i, j) = part->covariance[i][j]; 
            }
        }
        auto particle = std::make_shared<sct::ana::Particle>(part->pdg, momentum, covariance); // FIX IT RIGHT NOW!!!
       particle2pod[particle] = part; // for pod
        return particle;
    }
}

extern "C" {

POD::Particle* Particle_from_daughters(int pdg, int daughters_size, POD::Particle** daughters) { 
    return new POD::Particle(pdg, daughters_size, daughters);
}

POD::Particle* Particle_from_momentum(int pdg, double momentum[3], double covariance[3][3]) { 
    return new POD::Particle(pdg, momentum, covariance); 
}

double fit(POD::Particle* part) {
    std::map<sct::ana::ParticlePtr, POD::Particle*> particle2pod;
    sct::ana::ParticlePtr particle = createParticle(part, particle2pod);
    sct::ana::FitManager fitmanager(particle, {});
    fitmanager.fit();
    for (auto& [particle, pod]: particle2pod) {
        // update Python pod particles
        pod->momentum[0] = particle->momentum()[0];
        pod->momentum[1] = particle->momentum()[1];
        pod->momentum[2] = particle->momentum()[2];
    }
    return fitmanager.chiSquare();
}

}
