#include <iostream>
#include <stdio.h>
#include <map>

#include "dataobjects/Particle.h"

#include "treefitter/include/ConstraintConfiguration.h"
#include "treefitter/include/FitManager.h"
#include <iostream>


extern "C" {

sct::ana::FitManager* create_manager(sct::ana::Particle* part) {
    sct::ana::ParticlePtr p = sct::ana::ParticlePtr(part, [](sct::ana::Particle*){});
    return new sct::ana::FitManager(p, {});
}

void manager_fit(sct::ana::FitManager* fitmanager) {
    fitmanager->fit();
}

double* manager_momentum(sct::ana::FitManager* fitmanager, sct::ana::Particle* part) {
    sct::ana::ParticlePtr p = sct::ana::ParticlePtr(part, [](sct::ana::Particle*){});
    sct::FourVector fv = fitmanager->getMomentum(p);
    double* mom = new double[4];
    for (int i = 0; i < 4; ++i) mom[i] = fv(i);
    return mom;
}

double manager_chi(sct::ana::FitManager* fitmanager) {
    return fitmanager->chiSquare();
}

sct::ana::Particle* particle_from_momentum(int pdgCode, double momentum[3], double covariance[3][3]) {
    sct::ThreeVector mom; mom << momentum[0], momentum[1], momentum[2];
    Eigen::Matrix<double, 3, 3> cov; 
    for (int i = 0; i < 3; ++i) 
        for (int j = 0; j < 3; ++j)
            cov(i, j) = covariance[i][j];

    return new sct::ana::Particle(pdgCode, mom, cov);
}

sct::ana::Particle* particle_from_daughters(int pdgCode, int daughters_size, sct::ana::Particle** daughters) {
    std::vector<sct::ana::ParticlePtr> daughters_vector;
    for (int i = 0; i < daughters_size; ++i)
        daughters_vector.push_back(sct::ana::ParticlePtr(daughters[i], [](sct::ana::Particle* part){}));

    return new sct::ana::Particle(daughters_vector, pdgCode);
}

const double* get_momentum(sct::ana::Particle* particle) {
    return particle->fourMomentum().data();
}

}
