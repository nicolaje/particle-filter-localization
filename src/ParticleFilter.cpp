#include "ParticleFilter.h"

using namespace Eigen;
using namespace std;
using namespace octomap;

ParticleFilter::ParticleFilter() {
    this->distribution = std::normal_distribution<double>(0, 1);
    this->uniformDistribution = std::uniform_real_distribution<double>(0, 1);
    this->isInit = false;
    this->lastTime = false;
}

//*****************************************
// Private methods
//*****************************************

void ParticleFilter::normalize() {
    // Get the highest likelyhood
    double maxLW = logWeights[0];
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        if (logWeights[1] > maxLW) {
            maxLW = logWeights[i];
        }
    }

    // Using the log-sum of exponentials transform to avoid overflow
    // cf: http://lingpipe-blog.com/2009/06/25/log-sum-of-exponentials/
    double sumExp = 0;
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        sumExp += exp(logWeights[i]);
    }

    double logSumExp = maxLW + log(sumExp);

    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        logWeights[i] -= logSumExp;
    }
}

//*****************************************
// Public
//*****************************************

void ParticleFilter::init(const Eigen::Vector2d &initPos, const Eigen::Matrix2d &initPosCov) {
    vector<VectorXd> initPositions = this->drawSamples(PARTICLE_NUMBER, initPos, initPosCov);

    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        this->particles.col(i) = initPositions[i];
    }
    this->isInit = true;
}

bool ParticleFilter::isInitialized() {
    return this->isInit;
}

void ParticleFilter::setMaxRange(const double& range) {
    this->maxRange = range;
}

void ParticleFilter::setBathymetricFile(const std::string& filePath) {
    AbstractOcTree* read = OcTree::read(filePath);
    this->bathyMap = dynamic_cast<OcTree*> (read);
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
        
        double vRand=sqrt(uCov(0,0))*distribution(generator);
        double thetaRand=sqrt(uCov(1,1))*distribution(generator);
        
        Vector2d dX;
        // TODO: parallelize this loop
        for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
            dX<<
                    dt * vRand * cos(thetaRand * M_PI / 180.),
                    dt * vRand * sin(thetaRand * M_PI / 180.);
            this->particles.col(i) += dX;
        }
    }

    this->lastTime = t;
}

void ParticleFilter::update_walls(const double &rho, const double &rhoVar, const double &alpha, const double &alphaVar, const double &theta, const double &thetaVar) {
    double rhoLocal, alphaLocal,dist,err_sqr;
    double log_inv_sqrt = -0.5*2.5066282746310002 * sqrt(rhoVar);
    double inv_var = 1./(2*rhoVar);
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        rhoLocal = rho + sqrt(rhoVar) * distribution(generator);
        alphaLocal = alpha + sqrt(alphaVar) * distribution(generator);
        dist=getWallsRangeAt(particles.col(i),theta,alpha);
        err_sqr = pow(dist-rho,2);
        
        logWeights[i]+=log_inv_sqrt-err_sqr*inv_var;
    }
}

void ParticleFilter::update_range(const Vector2d &emitterPos, const double& range, const double& rangeVar) {
    double log_inv_sqrt = -0.5 * 2.5066282746310002 * sqrt(rangeVar);
    double inv_var = 1. / (2 * rangeVar);
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double err_sqr = pow((emitterPos - particles.col(i)).norm() - range, 2);

        // Update the log-odds
        logWeights[i] += log_inv_sqrt - err_sqr*inv_var;
    }
}

void ParticleFilter::update_sonar_vertical(const double& beamRange, const double& beamRangeVar, const double& beamAngle, const double& robotYaw, const double& z) {
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double xEnd = particles.col(i)[0] + beamRange * cos(robotYaw * M_PI / 180) * cos(beamAngle * M_PI / 180);
        double yEnd = particles.col(i)[1] + beamRange * sin(robotYaw * M_PI / 180) * cos(beamAngle * M_PI / 180);
        double zEnd = z + beamRange * sin(beamAngle * M_PI / 180);
        OcTreeNode* search = bathyMap->search(xEnd, yEnd, z);
        if (search != NULL) {
            this->logWeights[i] += search->getLogOdds();
        } else {
            // Decide what to do when sonar hits free space
            // Maybe cast a ray and then compute the error !
        }
    }
}

void ParticleFilter::update_sonar_horizontal(const double& beamRange, const double &beamRangeVar, const double &beamAngle, const double& z) {
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double xEnd = particles.col(i)[0] + beamRange * cos(beamAngle * M_PI / 180);
        double yEnd = particles.col(i)[1] + beamRange * sin(beamAngle * M_PI / 180);
        double zEnd = z;
        OcTreeNode* search = bathyMap->search(xEnd, yEnd, z);
        if (search != NULL) {
            this->logWeights[i] += search->getLogOdds();
        } else {
            // Decide what to do when sonar hits free space
            // Maybe cast a ray and then compute the error !
        }
    }
}

void ParticleFilter::update_GPS(const double& n, const double& nVar, const double& e, const double& eVar) {

}

void ParticleFilter::update_echosounder(const double& beam, const double& beamVar, const double &z) {
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        OcTreeNode* search = bathyMap->search(particles.col(i)[0], particles.col(i)[1], z);
        if (search != NULL) {
            this->logWeights[i] += search->getLogOdds();
        } else {
            // Decide what to do when sonar hits free space
            // Maybe cast a ray and then compute the error !
        }
    }
}

void ParticleFilter::resample() {
    // Draw PARTICLE_NUMBER uniformly distributed number between 0,1 and compute their log
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        uniforms[i] = uniformDistribution(generator);
    }

    // For each logUniform, compute the cumulative sum of the weights
    // and stop when the logUniform value is reached: keep this particle
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double cumSumTarget = uniforms[i];
        double cumSum = 0;
        unsigned int j = 0;
        for (j = 0; j < PARTICLE_NUMBER + 1 && cumSum < cumSumTarget; j++) {
            cumSum += exp(logWeights[j]);
        }
        if (j == PARTICLE_NUMBER) {
            cout << "Problem!" << endl;
            exit(EXIT_FAILURE);
            // Find a smarter way to admit things went wrong
        } else {
            particlesSwap.col(i) = particles.col(j);
            logWeightsSwap[i] = logWeights[j];
        }
    }

    particles = particlesSwap;
    memcpy(&logWeights[0], &logWeights[0], sizeof (double)*PARTICLE_NUMBER);
    //logWeights=logWeightsSwap;
}

Vector2d ParticleFilter::computeMean() {
    return particles.rowwise().mean();
}

Matrix2d ParticleFilter::computeCovariance() {
    Matrix<double, 2, PARTICLE_NUMBER> centered = particles.colwise() - particles.rowwise().mean();
    Matrix2d prod = centered*centered.transpose();
    return (prod) / double(particles.cols() - 1);
}