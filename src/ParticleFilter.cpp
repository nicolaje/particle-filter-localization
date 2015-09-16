#include "ParticleFilter.h"
#include <float.h>

using namespace Eigen;
using namespace std;
using namespace octomap;
int blah = 0;
int blah2 = 0;

ParticleFilter::ParticleFilter() {
    generator.seed(100);
    this->distribution = std::normal_distribution<double>(0, 1);
    this->uniformDistribution = std::uniform_real_distribution<double>(0, 1);

    this->isInit = false;
    this->lastTime = false;
    nbScanline = 0;
    scanLineUpdates = 0;
    double norm = 1. / PARTICLE_NUMBER;
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        weights[i] = norm;
    }
    wSlow = 0;
    wFast = 0;
}

void ParticleFilter::setResampleEvery(const int &nbScanLine) {
    this->nbScanline = nbScanLine;
}

void ParticleFilter::setResampleMethod(const int& method) {
    this->sampling_method = method;
}

void ParticleFilter::setAlphas(const double &alphaSlow, const double &alphaFast) {
    this->aSlow = alphaSlow;
    this->aFast = alphaFast;
}

void ParticleFilter::setWallsFile(const std::string &filePath) {
    std::ifstream in_file;
    in_file.open(filePath, ios::in);
    if (in_file.fail()) {
        std::stringstream s;
        s << "Simulator [load]: cannot open file " << filePath << "for reading the map";
        std::cerr << s.str() << std::endl;
        exit(-1);
    }
    Wall wall;
    while (in_file >> wall) {
        walls.push_back(wall);
    }
    in_file.close();
}

Eigen::Matrix<double, 2, PARTICLE_NUMBER> &ParticleFilter::getParticles() {
    return particles;
}

//*****************************************
// Private methods
//*****************************************

void ParticleFilter::normalize() {
    double sumW = weights[0];
    for (unsigned int i = 1; i < PARTICLE_NUMBER; i++) {
        sumW += weights[i];
    }

    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        weights[i] /= sumW;
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
    this->particlesInertial = this->particles;
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
    return mean + this->solver.operatorSqrt() * VectorXd::NullaryExpr(mean.rows(), normal);
}

vector<VectorXd> ParticleFilter::drawSamples(const int& nbParticle, const VectorXd& mean, const MatrixXd& cov) {
    vector<VectorXd> res;

    for (unsigned int i = 0; i < nbParticle; i++) {
        VectorXd toPush = this->drawSample(mean, cov);
        res.push_back(toPush);
    }

    return res;
}

double ParticleFilter::getWallsRangeAt(const Eigen::Vector2d &point, const double &theta) {
    /*double dist = 1000;

    for (unsigned int j = 0; j < walls.size(); j++) {
        Wall &w = walls[j];
        double dj, phij;
        CalcDistanceDirSegment(dj, phij, point[0], point[1], theta, w[0], w[1], w[2], w[3]);
        dist = (dj < dist) ? dj : dist;
    }
    return dist;*/

    double d = DBL_MAX;
    for (uint j = 0; j < walls.size(); j++) {
        Wall &w = walls[j];
        double dj, phij;
        DistanceDirSegment2(dj, phij, point[0], point[1], theta, w[0], w[1], w[2], w[3]);
        d = (dj < d) ? dj : d;
    }
    return d;
}

void ParticleFilter::predict(const double &t, const Vector2d &u, const Matrix2d &uCov) {
    if (!this->lastTimeinit) {
        this->lastTimeinit = true;
    } else {
        double dt = t - this->lastTime;
        double vRand;
        double thetaRand;
        Vector2d dX;
        for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
            vRand = sqrt(uCov(0, 0)) * distribution(generator) + u(0);
            thetaRand = (sqrt(uCov(1, 1)) * distribution(generator) + u(1)) * M_PI / 180.;
            dX <<
                    dt * vRand * cos(thetaRand),
                    dt * vRand * sin(thetaRand);
            this->particles.col(i) += dX;
            this->particlesInertial.col(i) += dX;
        }
    }

    this->lastTime = t;
}

void ParticleFilter::update_walls(const double &rho, const double &rhoVar, const double &alpha, const double &alphaVar) {

    scanLineUpdates += 1;
    blah2++;
    double rhoLocal, dist, err_sqr;
    double inv_sqrt = 1 / sqrt(rhoVar * 2 * M_PI);
    double inv_var = 1. / (2 * rhoVar);

    normalize();
    double wAvg = 0;
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        rhoLocal = rho + sqrt(rhoVar) * distribution(generator);
        dist = getWallsRangeAt(particles.col(i), alpha);
        err_sqr = pow((dist - rhoLocal), 2);
        double newW = inv_sqrt * exp(-err_sqr * inv_var);
        weights[i] *= newW;
        wAvg += weights[i] / PARTICLE_NUMBER;
    }
    wSlow += aSlow * (wAvg - wSlow);
    wFast += aFast * (wAvg - wFast);
    if (nbScanline == scanLineUpdates) {
        switch (sampling_method) {
            case LOW_VARIANCE:
                resampleLowVariance();
                break;
            case RESIDUAL:
                resampleResidual();
                break;
            case MULTINOMIAL:
                resample();
                break;
        }
        scanLineUpdates = 0;
    }
}

void ParticleFilter::update_range(const Vector2d &emitterPos, const double& range, const double& rangeVar) {
    double log_inv_sqrt = -0.5 * 2.5066282746310002 * sqrt(rangeVar);
    double inv_var = 1. / (2 * rangeVar);
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double err_sqr = pow((emitterPos - particles.col(i)).norm() - range, 2);

        // Update the log-odds
        //logWeights[i] += log_inv_sqrt - err_sqr*inv_var;
    }
}

void ParticleFilter::update_sonar_vertical(const double& beamRange, const double& beamRangeVar, const double& beamAngle, const double& robotYaw, const double& z) {
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double xEnd = particles.col(i)[0] + beamRange * cos(robotYaw * M_PI / 180) * cos(beamAngle * M_PI / 180);
        double yEnd = particles.col(i)[1] + beamRange * sin(robotYaw * M_PI / 180) * cos(beamAngle * M_PI / 180);
        double zEnd = z + beamRange * sin(beamAngle * M_PI / 180);
        OcTreeNode* search = bathyMap->search(xEnd, yEnd, z);
        if (search != NULL) {
            //this->logWeights[i] += search->getLogOdds();
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
            //this->logWeights[i] += search->getLogOdds();
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
            //this->logWeights[i] += search->getLogOdds();
        } else {
            // Decide what to do when sonar hits free space
            // Maybe cast a ray and then compute the error !
        }
    }
}

void ParticleFilter::resample() {
    normalize();
    bool initSqrtm = false;

    Vector2d mn;
    Matrix2d sqrtCov;
    
    this->bernoulliDistribution = bernoulli_distribution(1 - wFast / wSlow > 0 ? 1 - wFast / wSlow : 0);

    double cumulativeSum[PARTICLE_NUMBER];

    cumulativeSum[0] = weights[0];

    for (unsigned int i = 1; i < PARTICLE_NUMBER; i++) {
        cumulativeSum[i] = cumulativeSum[i - 1] + weights[i];
    }

    // For each logUniform, compute the cumulative sum of the weights
    // and stop when the logUniform value is reached: keep this particle
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double cumSumTarget = uniformDistribution(generator);
        unsigned int idx = 0;

        while (cumulativeSum[idx] < cumSumTarget) {
            idx++;
        }
        if (idx >= PARTICLE_NUMBER) {
            cout << "Problem!" << endl;
            cout << "CumSum[" << idx << "] : " << cumulativeSum[idx] << endl;
            cout << "CumSumTarget: " << cumSumTarget << endl;
            for (unsigned int k = 0; k < PARTICLE_NUMBER; k++) {
                cout << "weights[" << k << "]: " << weights[k] << ", w[" << k << "]= " << weights[k] << endl;
            }
            exit(EXIT_FAILURE);
            // Find a smarter way to admit things went wrong
        } else {
            if (bernoulliDistribution(generator) || weights[i] == 0) {
                if (!initSqrtm) {
                    initSqrtm = true;
                    SelfAdjointEigenSolver<Matrix2d> s2d(computeCovariance());
                    sqrtCov = s2d.operatorSqrt();
                }
                particlesSwap.col(i) = computeWeightedMean() + sqrtCov * (Vector2d(distribution(generator), distribution(generator)));
            } else {
                particlesSwap.col(i) = particles.col(idx);
            }
            weightsSwap[i] = 1;
        }
    }
    particles = particlesSwap;
    memcpy(&weights[0], &weightsSwap[0], sizeof (double)*PARTICLE_NUMBER);
}

void ParticleFilter::resampleResidual() {
    this->normalize();
    bool initSqrtm = false;

    Vector2d mn;
    Matrix2d sqrtCov;

    // Nb of copies of each particle
    int nbCopies[PARTICLE_NUMBER];
    int indexes[PARTICLE_NUMBER];

    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        nbCopies[i] = abs((int) (weights[i] * PARTICLE_NUMBER));
    }
    int k = 0;
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        for (unsigned int j = 0; j < nbCopies[i]; j++) {
            indexes[k] = i;
            k++;
        }
    }

    double residual[PARTICLE_NUMBER];
    double residualSum = 0;

    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        residual[i] = weights[i] * PARTICLE_NUMBER - nbCopies[i];
        residualSum += residual[i];
    }

    // Normalize the residuals
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        residual[i] /= residualSum;
    }

    // multinomial resampling on the residual
    double cumulativeSum[PARTICLE_NUMBER];
    cumulativeSum[0] = residual[0];
    for (unsigned int i = 1; i < PARTICLE_NUMBER; i++) {
        cumulativeSum[i] = cumulativeSum[i - 1] + residual[i];
    }

    for (unsigned int i = 0; i < PARTICLE_NUMBER - k; i++) {
        double cumSumTarget = min(uniformDistribution(generator), cumulativeSum[PARTICLE_NUMBER - 1]);
        unsigned int idx = 0;

        while (cumulativeSum[idx] < cumSumTarget) {
            idx++;
        }
        indexes[k + i] = idx;
    }

    this->bernoulliDistribution = bernoulli_distribution(1 - wFast / wSlow > 0 ? 1 - wFast / wSlow : 0);

    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        if (bernoulliDistribution(generator) || weights[i] == 0) {
            if (!initSqrtm) {
                initSqrtm = true;
                SelfAdjointEigenSolver<Matrix2d> s2d(computeCovariance());
                sqrtCov = s2d.operatorSqrt();
            }
            particlesSwap.col(i) = computeWeightedMean() + sqrtCov * (Vector2d(distribution(generator), distribution(generator)));
        } else {
            particlesSwap.col(i) = particles.col(indexes[i]);
        }
        weightsSwap[i] = 1;
    }
    particles = particlesSwap;
    memcpy(&weights[0], &weightsSwap[0], sizeof (double)*PARTICLE_NUMBER);
}

void ParticleFilter::resampleLowVariance() {
    normalize();
    
    bool initSqrtm = false;

    Vector2d mn;
    Matrix2d sqrtCov;

    double r = uniformDistribution(generator) / PARTICLE_NUMBER;
    double c = weights[0];
    int idx = 0;
    for (unsigned int i = 1; i < PARTICLE_NUMBER; i++) {
        double U = r + ((double) (i - 1)) / PARTICLE_NUMBER;
        while (U > c) {
            idx++;
            c += weights[idx];
        }
        if (bernoulliDistribution(generator) || weights[i] == 0) {
            if (!initSqrtm) {
                initSqrtm = true;
                SelfAdjointEigenSolver<Matrix2d> s2d(computeCovariance());
                sqrtCov = s2d.operatorSqrt();
            }
            particlesSwap.col(i) = computeWeightedMean() + sqrtCov * (Vector2d(distribution(generator), distribution(generator)));
        } else {
        particlesSwap.col(i) = particles.col(idx);
        }
        weightsSwap[i] = 1;
    }
    particles = particlesSwap;
    memcpy(&weights[0], &weightsSwap[0], sizeof (double)*PARTICLE_NUMBER);
}

Vector2d ParticleFilter::computeMean() {
    this->normalize();
    return particles.rowwise().mean();
}

Matrix2d ParticleFilter::computeCovariance() {
    this->normalize();
    Matrix<double, 2, PARTICLE_NUMBER> centered = particles.colwise() - particles.rowwise().mean();
    Matrix2d prod = centered * centered.transpose();
    return (prod) / double(particles.cols() - 1);
}

Vector2d ParticleFilter::computeWeightedMean() {
    normalize();
    Map<Matrix<double, 1, PARTICLE_NUMBER> > weightVec(weights);
    Vector2d res = (weightVec * particles.transpose()).transpose();
    return res;
}

Matrix2d ParticleFilter::computeWeightedCovariance() {
    normalize();
    
    Vector2d wMean=computeWeightedMean();
    Map<Matrix<double, 1, PARTICLE_NUMBER> > weightVec(weights);
    
    double sumSquare = weightVec.squaredNorm();
    
    //Matrix<double, 2, PARTICLE_NUMBER> centered = particles.colwise() - particles*(weightVec.transpose());//.rowwise().mean();
    
    Matrix2d res=Matrix2d::Zero();;//(centered*(weightVec.cwiseSqrt()).transpose())*(weightVec.cwiseSqrt()*centered.transpose());//
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        res(0,0)+=weights[i]*(particles(0,i)-wMean[0])*(particles(0,i)-wMean[0]);
    }
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        res(1,1)+=weights[i]*(particles(1,i)-wMean[1])*(particles(1,i)-wMean[1]);
    }
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        res(0,1)+=weights[i]*(particles(0,i)-wMean[0])*(particles(1,i)-wMean[1]);
    }
    res(1,0)=res(0,1);
    res*=(1/(1-sumSquare));
    return res;
}