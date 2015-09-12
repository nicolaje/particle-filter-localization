#include "ParticleFilter.h"

using namespace Eigen;
using namespace std;
using namespace octomap;

ParticleFilter::ParticleFilter() {
    this->distribution = std::normal_distribution<double>(0, 1);
}

//*****************************************
// Private methods
//*****************************************

void ParticleFilter::normalize()
{
    // Get the highest likelyhood
    double max=*max_element(logWeights.begin(),logWeights.end());
    
    // Using the log-sum of exponentials transform to avoid overflow
    // cf: http://lingpipe-blog.com/2009/06/25/log-sum-of-exponentials/
    double sumExp=0;
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        sumExp+=exp(logWeights[i]);
    }
    
    double logSumExp=max+log(sumExp);
    
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        logWeights-=logSumExp;
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
        logWeights[i] += log_inv_sqrt - err_sqr*inv_var;
    }
}

void ParticleFilter::update_sonar_vertical(const double& beamRange, const double& beamRangeVar, const double& beamAngle, const double& robotYaw, const double& z)
{
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        double xEnd=particles.col(i)[0]+beamRange*cos(robotYaw*M_PI/180)*cos(beamAngle*M_PI/180);
        double yEnd=particles.col(i)[1]+beamRange*sin(robotYaw*M_PI/180)*cos(beamAngle*M_PI/180);
        double zEnd=z+beamRange*sin(beamAngle*M_PI/180);
        OcTreeNode* search = bathyMap->search(xEnd, yEnd, z);
        if (search != NULL) {
            this->logWeights[i]+=search->getLogOdds();
        }
        else{
            // Decide what to do when sonar hit free space
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
            this->logWeights[i]+=search->getLogOdds();
        }
        else{
            // Decide what to do when sonar hit free space
            // Maybe cast a ray and then compute the error !
        }
    }
}

void ParticleFilter::update_echosounder(const double& beam, const double& beamVar, const double &z) {
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        OcTreeNode* search=bathyMap->search(particles.col(i)[0],particles.col(i)[1],z);
        if(search!=NULL){
            this->logWeights[i]+=search->getLogOdds();
        }
        else{
            // Decide what to do when sonar hit free space
            // Maybe cast a ray and then compute the error !
        }
    }
}