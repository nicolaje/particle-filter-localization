#include "ParticleFilter.h"

using namespace Eigen;
using namespace std;
using namespace octomap;
int blah=0;

ParticleFilter::ParticleFilter() {
    generator.seed(100);
    this->distribution = std::normal_distribution<double>(0, 1);
    this->uniformDistribution = std::uniform_real_distribution<double>(0, 1);
    this->isInit = false;
    this->lastTime = false;
    double norm = 1. / PARTICLE_NUMBER;
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        weights[i] = norm;
    }
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
    cout << "debug"<<endl;
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
    double dist = 1000;

    for (unsigned int j = 0; j < walls.size(); j++) {
        Wall &w = walls[j];
        double dj, phij;
        CalcDistanceDirSegment(dj, phij, point[0], point[1], theta, w[0], w[1], w[2], w[3]);
        dist = (dj < dist) ? dj : dist;
    }
    return dist;
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
        }
    }

    this->lastTime = t;
}

void ParticleFilter::update_walls(const double &rho, const double &rhoVar, const double &alpha, const double &alphaVar) {
    double rhoLocal, dist, err_sqr;
    double inv_sqrt = 1 / sqrt(rhoVar * 2 * M_PI);
    double inv_var = 1. / (2 * rhoVar);
    if(blah>=84)
    {
        cout << "debug me"<<endl;
    }
    this->normalize();
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        rhoLocal = rho + sqrt(rhoVar) * distribution(generator);
        dist = getWallsRangeAt(particles.col(i), alpha);
        err_sqr = pow((dist - rhoLocal), 2);
        double newW = inv_sqrt * exp(err_sqr * inv_var);
        weights[i] *= newW;
        /*std::cout << "dist: "<<dist<<std::endl;
        std::cout << "rho: "<<rho<<std::endl;
        std::cout << "alpha: "<<alpha*180./M_PI<<std::endl;
        std::cout << "alpha%360: "<<fmod(((alpha*180./M_PI)),360)<<std::endl;
        std::cout << "err: "<<sqrt(err_sqr)<<std::endl;
        std::cout << "particle: "<<std::endl<<particles.col(i)<<std::endl;
        cout << "inv_sqrt= "<<inv_sqrt<<endl;
        cout << "inv_var= "<<inv_var<<endl;
        cout << "exp(err_sqr*inv_var): "<<exp(err_sqr*inv_var)<<endl;
        std::cout << "loGweight["<<i<<"]: "<<weights[i]<<std::endl;*/
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
            particlesSwap.col(i) = particles.col(idx);
            weightsSwap[i] = weights[idx];
        }
    }
    particles = particlesSwap;
    memcpy(&weights[0], &weightsSwap[0], sizeof (double)*PARTICLE_NUMBER);
}
void ParticleFilter::resampleResidual() {
    blah++;
    if(blah>=84)
    {
        cout << "debug me"<<endl;
    }
    this->normalize();
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
        double cumSumTarget = min(uniformDistribution(generator),cumulativeSum[PARTICLE_NUMBER-1]);
        unsigned int idx = 0;

        while (cumulativeSum[idx] < cumSumTarget) {
            idx++;
        }
        indexes[k + i] = idx;
    }

    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        particlesSwap.col(i) = particles.col(indexes[i]);
        weightsSwap[i] = weights[indexes[i]];
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
    Map<Matrix<double, 1, PARTICLE_NUMBER> > weightVec(weights);
    double sum = weightVec.sum();
    double sumSquare = weightVec.squaredNorm();
    double squareSum = sum*sum;
    // TODO
    //Matrix2d res = 
}