#include "v3.h"
#include <cmath>
#include <iostream>

// Source - https://stackoverflow.com/a/1727896
// Posted by Ferenc Deak, modified by community. See post 'Timeline' for change
// history Retrieved 2026-05-02, License - CC BY-SA 3.0
#define _USE_MATH_DEFINES
#include <math.h>

struct SimulationParameters {
  public:
    // Simulation Relevant Parameters (once I am done at least)
    double volumeParticle;
    double dragRadius; // can kill
    double dragCoeffTransl;
    double dragCoeffRot;
    double chiEffLongAxesAB;
    double chiEffShortAxisC;
    double chiEffShapeAnisotropyFactor; // can kill
    double totalMagDipoleMomentParticle;
    double lengthSimulationCube;

    // Experiment Parameters
    int numberOfParticles = 100;
    double simulationTime = 1000; // ms
    double volumeFraction = 0.01; // 1.0 = 100%
    double viscosity = 3.5;       // Pa*s
    double magnitudeMagFieldExternal = 0;
    double magnitudeElFieldExternal = 0;

    // Particle/Matrix Parameters
    double magMomentDensityParticle = 0; // can kill
    double aspectRatioParticle = 7 / 2;  // Width/Height
    double longSemiaxesAB = 2.5;         // um
    // double shortSemiaxisC;
    //   v3 orientationParticle;
    double relPermittivityParticle = 10; // can kill?
    double relPermittivityMatrix = 2;    // can kill?

    // Correction Factors
    double corrFactorRepulsiveForce = 40;
    double corrFactorVelocity = 1 / 3;

    /*
    void set_numberOfParticles(int val_numberOfParticles) {
        numberOfParticles = val_numberOfParticles;
    };
    void set_simulationTime(double val_simulationTime) {
        simulationTime = val_simulationTime;
    };
    void set_volumeFraction(double val_volumeFraction) {
        volumeFraction = val_volumeFraction;
    };
    void set_viscosity(double val_viscosity) { viscosity = val_viscosity; };
    void set_magnitudeMagFieldExternal(double val_magnitudeMagFieldExternal) {
        magnitudeMagFieldExternal = val_magnitudeMagFieldExternal;
    };
    void set_magnitudeElFieldExternal(double val_magnitudeElFieldExternal) {
        magnitudeElFieldExternal = val_magnitudeElFieldExternal;
    };

    void set_magMomentDensityParticle(double val_magMomentDensityParticle) {
        magMomentDensityParticle = val_magMomentDensityParticle;
    };
    void set_aspectRatioParticle(double val_aspectRatioParticle) {
        aspectRatioParticle = val_aspectRatioParticle;
    };
    void set_longSemiaxesAB(double val_longSemiaxesAB) {
        longSemiaxesAB = val_longSemiaxesAB;
    };
    void set_relPermittivityParticle(double val_relPermittivityParticle) {
        relPermittivityParticle = val_relPermittivityParticle;
    };
    void set_relPermittivityMatrix(double val_relPermittivityMatrix) {
        relPermittivityMatrix = val_relPermittivityMatrix;
    };

    void set_corrFactorRepulsiveForce(double val_corrFactorRepulsiveForce) {
        corrFactorRepulsiveForce = val_corrFactorRepulsiveForce;
    };
    void set_corrFactorVelocity(double val_corrFactorVelocity) {
        corrFactorVelocity = val_corrFactorVelocity;
    };
    */
};

class SimBuilder {
  private:
    SimulationParameters parameters;

  public:
    // SimBuilder () {parameters = SimulationParameters();};

    SimBuilder &set_numParticles(int numberOfParticles) {
        parameters.numberOfParticles = numberOfParticles;
        return *this;
    };
    SimBuilder &set_duration(double simulationTime) {
        parameters.simulationTime = simulationTime;
        return *this;
    };
    SimBuilder &set_volFrac(double volumeFraction) {
        parameters.volumeFraction = volumeFraction;
        return *this;
    };
    SimBuilder &set_viscosity(double viscosity) {
        parameters.viscosity = viscosity;
        return *this;
    };
    SimBuilder &set_magField(double magnitudeMagFieldExternal) {
        parameters.magnitudeMagFieldExternal = magnitudeMagFieldExternal;
        return *this;
    };
    SimBuilder &set_elField(double magnitudeElFieldExternal) {
        parameters.magnitudeElFieldExternal = magnitudeElFieldExternal;
        return *this;
    };

    SimBuilder &set_magMomentDensityParticle(double magMomentDensityParticle) {
        parameters.magMomentDensityParticle = magMomentDensityParticle;
        return *this;
    };
    SimBuilder &set_aspectRatio(double aspectRatioParticle) {
        parameters.aspectRatioParticle = aspectRatioParticle;
        return *this;
    };
    SimBuilder &set_longSemiaxesAB(double longSemiaxesAB) {
        parameters.longSemiaxesAB = longSemiaxesAB;
        return *this;
    };
    SimBuilder &set_relPermittivityParticle(double relPermittivityParticle) {
        parameters.relPermittivityParticle = relPermittivityParticle;
        return *this;
    };
    SimBuilder &set_relPermittivityMatrix(double relPermittivityMatrix) {
        parameters.relPermittivityMatrix = relPermittivityMatrix;
        return *this;
    };

    SimBuilder &set_corrFactorRepulsiveForce(double corrFactorRepulsiveForce) {
        parameters.corrFactorRepulsiveForce = corrFactorRepulsiveForce;
        return *this;
    };
    SimBuilder &set_corrFactorVelocity(double corrFactorVelocity) {
        parameters.corrFactorVelocity = corrFactorVelocity;
        return *this;
    };

    SimulationParameters build() {
        double shortSemiaxisC =
            parameters.longSemiaxesAB / parameters.aspectRatioParticle;
        parameters.volumeParticle = 4.0 / 3.0 * M_PI *
                                    std::pow(parameters.longSemiaxesAB, 2) *
                                    shortSemiaxisC;

        parameters.dragRadius =
            std::cbrt(std::pow(parameters.longSemiaxesAB, 2) * shortSemiaxisC);

        parameters.dragCoeffTransl =
            6 * M_PI * parameters.viscosity * parameters.dragRadius;

        parameters.dragCoeffRot = 8 * M_PI * parameters.viscosity *
                                  std::pow(parameters.dragRadius, 3);

        parameters.chiEffShapeAnisotropyFactor =
            (std::pow(parameters.longSemiaxesAB, 2) * shortSemiaxisC) /
            (2 * (std::pow(parameters.longSemiaxesAB, 2) -
                  std::pow(shortSemiaxisC, 2))) *
            ((M_PI / 2 *
              std::sqrt(std::pow(parameters.longSemiaxesAB, 2) -
                        std::pow(shortSemiaxisC, 2))) -
             shortSemiaxisC / std::pow(parameters.longSemiaxesAB, 2));

        parameters.chiEffLongAxesAB =
            (parameters.relPermittivityParticle -
             parameters.relPermittivityMatrix) /
            (1 + parameters.chiEffShapeAnisotropyFactor *
                     (parameters.relPermittivityParticle -
                      parameters.relPermittivityMatrix) /
                     parameters.relPermittivityMatrix);

        parameters.chiEffShortAxisC =
            (parameters.relPermittivityParticle -
             parameters.relPermittivityMatrix) /
            (1 + (1 - 2 * parameters.chiEffShapeAnisotropyFactor) *
                     (parameters.relPermittivityParticle -
                      parameters.relPermittivityMatrix) /
                     parameters.relPermittivityMatrix);

        parameters.totalMagDipoleMomentParticle =
            parameters.volumeParticle * parameters.magMomentDensityParticle;

        parameters.lengthSimulationCube =
            std::cbrt(parameters.numberOfParticles * parameters.volumeParticle /
                      parameters.volumeFraction);

        return parameters;
    }
};

/*
int main() {
    SimBuilder builder1;
    builder1.set_aspectRatio(40);
    SimulationParameters test1 = builder1.build();

    SimBuilder builder2;
    SimulationParameters test2 = builder2.build();

    std::cout << test2.volumeParticle << std::endl;
    std::cout << test1.volumeParticle << std::endl;

    builder1.set_longSemiaxesAB(5);
    SimulationParameters test3 = builder1.build();

    std::cout << test3.volumeParticle << std::endl;
}
*/