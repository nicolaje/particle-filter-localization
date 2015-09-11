#include "ParticleFilter.h"

using namespace Eigen;
using namespace std;

ParticleFilter::ParticleFilter() {
    this->distribution = std::normal_distribution<double>(0, 1);
}

void ParticleFilter::init(const Eigen::Vector2d &initPos, const Eigen::Matrix2d &initPosCov) {
    vector<VectorXd> initPositions = this->drawSamples(PARTICLE_NUMBER, initPos, initPosCov);

    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        this->particles.col(i) = initPositions[i];
    }
}

VectorXd ParticleFilter::drawSample(const Eigen::VectorXd& mean, const Eigen::MatrixXd& cov) {
    // eigen.tuxfamily.org/bz/show_bug.cgi?id=720
    auto normal = [&](double) {
        return distribution(generator);
    };

    // Cf: Robotics lesson
    this->solver = SelfAdjointEigenSolver<MatrixXd>(cov);
    return this->solver.operatorSqrt() * VectorXd::NullaryExpr(mean.rows(), normal);
}

vector<VectorXd> ParticleFilter::drawSamples(const int& nbParticle, const VectorXd& mean, const MatrixXd& cov) {
    vector<VectorXd> res;

    for (unsigned int i = 0; i < nbParticle; i++) {
        VectorXd toPush = this->drawSample(mean, cov);
        res.push_back(toPush);
    }

    return res;
}

void ParticleFilter::predict(const double &t, const Vector2d &u, const Matrix2d &uCov) {
    if (!this->lastTimeinit) {
        this->lastTimeinit = true;
    } else {
        double dt = t - this->lastTime;

        vector<VectorXd> inputs = this->drawSamples(PARTICLE_NUMBER, u, uCov);

        // TODO: parallelize this loop
        for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
            Vector2d dX(dt * inputs[i](0) * cos(inputs[i](1) * M_PI / 180.), dt * inputs[i](0) * sin(inputs[i](1) * M_PI / 180.));
            this->particles.col(i) += dX;
        }
    }

    this->lastTime = t;
}

void ParticleFilter::update_range(const Vector2d &emitterPos, const double& range, const double& rangeVar) {
    double log_inv_sqrt = -0.5 * 2.5066282746310002 * sqrt(rangeVar);
    double inv_var = 1. / (2 * rangeVar);
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double err_sqr = pow((emitterPos - particles.col(i)).norm() - range, 2);
        // Update the log-odds
        weights[i] += log_inv_sqrt - err_sqr*inv_var;
    }
}