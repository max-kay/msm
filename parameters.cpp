#include "v3.h"
#include <cmath>
#include <iostream>

// Source - https://stackoverflow.com/a/1727896
// Posted by Ferenc Deak, modified by community. See post 'Timeline' for change
// history Retrieved 2026-05-02, License - CC BY-SA 3.0
#define _USE_MATH_DEFINES
#include <math.h>

struct SimBuilder {

  private:
    SimulationParameters parameters;

  public:
    void set_numberOfParticles(int numberOfParticles) {
        parameters.numberOfParticles = numberOfParticles;
    };
    void set_simulationTime(double simulationTime) {
        parameters.simulationTime = simulationTime;
    };
    void set_volumeFraction(double volumeFraction) {
        parameters.volumeFraction = volumeFraction;
    };
    void set_viscosity(double viscosity) { parameters.viscosity = viscosity; };
    void set_magnitudeMagFieldExternal(double magnitudeMagFieldExternal) {
        parameters.magnitudeMagFieldExternal = magnitudeMagFieldExternal;
    };
    void set_magnitudeMagFieldExternal(double magnitudeMagFieldExternal) {
        parameters.magnitudeMagFieldExternal = magnitudeMagFieldExternal;
    };

    void set_magMomentDensityParticle(double magMomentDensityParticle) {
        parameters.magMomentDensityParticle = magMomentDensityParticle;
    };
    void set_aspectRatioParticle(double aspectRatioParticle) {
        parameters.aspectRatioParticle = aspectRatioParticle;
    };
    void set_longSemiaxesAB(double longSemiaxesAB) {
        parameters.longSemiaxesAB = longSemiaxesAB;
    };
    void set_relPermittivityParticle(double relPermittivityParticle) {
        parameters.relPermittivityParticle = relPermittivityParticle;
    };
    void set_relPermittivityMatrix(double relPermittivityMatrix) {
        parameters.relPermittivityMatrix = relPermittivityMatrix;
    };

    void set_corrFactorRepulsiveForce(double corrFactorRepulsiveForce) {
        parameters.corrFactorRepulsiveForce = corrFactorRepulsiveForce;
    };
    void set_corrFactorVelocity(double corrFactorVelocity) {
        parameters.corrFactorVelocity = corrFactorVelocity;
    };

    void build() {
        double shortSemiaxisC =
            parameters.longSemiaxesAB / parameters.aspectRatioParticle;
        parameters.volumeParticle = 4 / 3 * M_PI *
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

        return /*???!!!!!!???*/;
    }
};

struct SimulationParameters {
  public:
    // Simulation Relevant Parameters
    double volumeParticle;
    double dragRadius;
    double dragCoeffTransl;
    double dragCoeffRot;
    double chiEffLongAxesAB;
    double chiEffShortAxisC;
    double chiEffShapeAnisotropyFactor;
    double totalMagDipoleMomentParticle;
    double lengthSimulationCube;

    // Experimental Parameters
    int numberOfParticles = 100;
    double simulationTime = 1000; // ms
    double volumeFraction = 0.01; // 1.0 = 100%
    double viscosity = 3.5;       // Pa*s
    double magnitudeMagFieldExternal = 0;
    double magnitudeElFieldExternal = 0;

    // Particle/Matrix Parameters
    double magMomentDensityParticle = 0;
    double aspectRatioParticle = 7 / 2; // Width/Height
    double longSemiaxesAB = 2.5;        // um
    // double shortSemiaxisC;
    //   v3 orientationParticle;
    double relPermittivityParticle = 10;
    double relPermittivityMatrix = 2;

    // Correction Factors
    double corrFactorRepulsiveForce = 40;
    double corrFactorVelocity = 1 / 3;
};