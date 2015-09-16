#include "vibes.h"
#include "simulator.h"
#include "ParticleFilter.h"
#include "../Eigen/Dense"
#include <unistd.h>
#include <cmath>

using namespace Eigen;
using namespace std;

int main(int argc, char** argv) {

    // Set Vibes Params
    vibes::beginDrawing();
    vibes::newFigure("test");
    vibes::setFigureProperties("test", vibesParams("x", 0, "y", 0, "width", 1000, "height", 1000));
    vibes::axisLimits(-10, 140, -10, 140);

    Simulator simu("map.txt", "traj.txt");
    simu.setSonarFreq(5); // 5hz
    simu.setNbOutliers(0); // 5 % outliers

    ParticleFilter pf;
    pf.setWallsFile("map.txt");
    Pose &p0 = simu.currentPose();

    Vector2d p(p0.x, p0.y);
    Matrix2d pCov = Matrix2d::Identity();

    pf.init(p, pCov);

    Eigen::Vector2d u;
    Eigen::Matrix2d uCov = Matrix2d::Identity();
    uCov(0, 0) = 1 * 1;
    uCov(1, 1) = 7 * 7;

    // Resample every 10 scanlines
    pf.setResampleEvery(20);
    pf.setResampleMethod(RESIDUAL);
    pf.setAlphas(0.1,10);
    simu.drawMap();
    vibes::newGroup("greenlines", "green");
    vibes::newGroup("redlines", "red");

    vibes::newGroup("particles", vibesParams("color", "red"));

    vibes::newGroup("AUV", vibesParams("color", "black"));

    while (simu.nextStep() != -1) {

        // Current point of the trajectory
        Pose &p = simu.currentPose();

        // Draw the current situation
        u << p.speed, p.theta * 180 / M_PI;
        pf.predict(p.t, u, uCov);

        double d = simu.genSonarDist(p.x, p.y, p.theta, p.t);
        if (d >= 0) { // the measurment is valid
            vibes::clearGroup("particles");

            Eigen::Matrix<double, 2, PARTICLE_NUMBER> particles = pf.getParticles();
            for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
                vibes::drawCircle(particles.col(i)[0], particles.col(i)[1], 0.01, vibesParams("group", "particles"));
            }
            double alpha = simu.getSonarAngle(p.t) + p.theta; // sonar beam angle
            pf.update_walls(d, 25, alpha, 1);

            Vector2d mean = pf.computeWeightedMean();
            Matrix2d cov = pf.computeWeightedCovariance();
            vibes::drawConfidenceEllipse(mean[0], mean[1], cov(0, 0), cov(0, 1), cov(1, 1), 3, vibesParams("group", "particles", "color", "red"));
            //vibes::drawCircle(mean[0], mean[1], 5, vibesParams("group", "particles", "color", "red"));

            //Vector2d wMean = pf.computeWeightedMean();
            //vibes::drawCircle(wMean[0], wMean[1], 5, vibesParams("group", "particles", "color", "red"));
        }
        vibes::clearGroup("AUV");
        vibes::drawAUV(p.x, p.y, p.theta * 180 / M_PI, 1, vibesParams("group", "AUV"));
        usleep(1000);
    }
    return 0;
}
