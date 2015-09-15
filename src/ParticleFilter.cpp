#include "ParticleFilter.h"

using namespace Eigen;
using namespace std;
using namespace octomap;

ParticleFilter::ParticleFilter() {
    generator.seed(100);
    this->distribution = std::normal_distribution<double>(0, 1);
    this->uniformDistribution = std::uniform_real_distribution<double>(0, 1);
    this->isInit = false;
    this->lastTime = false;
    double norm=1./PARTICLE_NUMBER;
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        weights[i]=norm;
    }
}

void ParticleFilter::setWallsFile(const std::string &filePath)
{
    std::ifstream in_file;
    in_file.open(filePath, ios::in);
    if(in_file.fail()) {
        std::stringstream s;
        s << "Simulator [load]: cannot open file " << filePath << "for reading the map";
        std::cerr << s.str() << std::endl;
        exit(-1);
    }
    Wall wall;
    while(in_file >> wall){
        walls.push_back(wall);
    }
    in_file.close();
}

Eigen::Matrix<double, 2, PARTICLE_NUMBER> &ParticleFilter::getParticles()
{
    return particles;
}

//*****************************************
// Private methods
//*****************************************

void ParticleFilter::normalize() {
    // Get the highest likelyhood
    //double maxLW = logWeights[0];
    //for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
    //    if (logWeights[1] > maxLW) {
    //        maxLW = logWeights[i];
    //    }
    //}

    // Using the log-sum of exponentials transform to avoid overflow
    // cf: http://lingpipe-blog.com/2009/06/25/log-sum-of-exponentials/
    //double sumExp = exp(logWeights[0]-maxLW);
    //for (unsigned int i = 1; i < PARTICLE_NUMBER; i++) {
    //    sumExp += exp(logWeights[i]-maxLW);
    //}

    //double logSumExp = maxLW + log(sumExp);

    //for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
    //    logWeights[i] -= logSumExp;
    //}
    double sumW=weights[0];
    for(unsigned int i=1;i<PARTICLE_NUMBER;i++)
    {
        sumW+=weights[i];
    }
    for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
    {
        weights[i]/=sumW;
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
    return mean+ this->solver.operatorSqrt() * VectorXd::NullaryExpr(mean.rows(), normal);
}

vector<VectorXd> ParticleFilter::drawSamples(const int& nbParticle, const VectorXd& mean, const MatrixXd& cov) {
    vector<VectorXd> res;

    for (unsigned int i = 0; i < nbParticle; i++) {
        VectorXd toPush = this->drawSample(mean, cov);
        res.push_back(toPush);
    }

    return res;
}

double ParticleFilter::getWallsRangeAt(const Eigen::Vector2d &point, const double &theta)
{
    double dist=1000;
    
    for(unsigned int j=0;j<walls.size();j++)
        {
            Wall &w=walls[j];
            double dj,phij;
            CalcDistanceDirSegment(dj,phij,point[0],point[1],theta,w[0],w[1],w[2],w[3]);
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
            vRand=sqrt(uCov(0,0))*distribution(generator)+u(0);
            thetaRand=(sqrt(uCov(1,1))*distribution(generator)+u(1))*M_PI/180.;
            dX<<
                    dt*vRand * cos(thetaRand),
                    dt*vRand * sin(thetaRand);
            this->particles.col(i) += dX;
        }
    }

    this->lastTime = t;
}

void ParticleFilter::update_walls(const double &rho, const double &rhoVar, const double &alpha, const double &alphaVar) {
    double rhoLocal, dist,err_sqr;
    double inv_sqrt = 1/sqrt(rhoVar*2*M_PI);
    double inv_var = 1./(2*rhoVar);
    std::cout<<"update_walls"<<std::endl;
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        rhoLocal = rho + sqrt(50*50) * distribution(generator);
        dist=getWallsRangeAt(particles.col(i),alpha);
        err_sqr = pow((dist-rhoLocal),2);
        weights[i]*=inv_sqrt*exp(err_sqr*inv_var);
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
        this->normalize();
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
    // Draw PARTICLE_NUMBER uniformly distributed number between 0,1 and compute their log
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        uniforms[i] = uniformDistribution(generator);
    }

    // For each logUniform, compute the cumulative sum of the weights
    // and stop when the logUniform value is reached: keep this particle
    for (unsigned int i = 0; i < PARTICLE_NUMBER; i++) {
        double cumSumTarget = uniforms[i];
        double cumSum = weights[0];
        unsigned int idx = 0;
        
        while(cumSum<cumSumTarget)
        {
            idx++;
            cumSum+=weights[idx];
        }
        if (idx >= PARTICLE_NUMBER) {
            cout << "Problem!" << endl;
            cout << "CumSum: "<<cumSum<<endl;
            cout << "CumSumTarget: "<<cumSumTarget<<endl;
            for(unsigned int k=0;k<PARTICLE_NUMBER;k++)
            {
                cout << "weights["<<k<<"]: "<<weights[k]<<", w["<<k<<"]= "<<weights[k]<<endl;
            }
            exit(EXIT_FAILURE);
            // Find a smarter way to admit things went wrong
        } else {
            particlesSwap.col(i) = particles.col(idx);
            weightsSwap[i] = weights[idx];
        }
    }
    particles = particlesSwap;
    memcpy(&weights[0], &weights[0], sizeof (double)*PARTICLE_NUMBER);
}

Vector2d ParticleFilter::computeMean() {
    return particles.rowwise().mean();
}

Matrix2d ParticleFilter::computeCovariance() {
    Matrix<double, 2, PARTICLE_NUMBER> centered = particles.colwise() - particles.rowwise().mean();
    Matrix2d prod = centered*centered.transpose();
    return (prod) / double(particles.cols() - 1);
}