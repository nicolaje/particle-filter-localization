#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#define PARTICLE_NUMBER 500 // Number of particles

#include <random>
#include "../Eigen/Dense"

class ParticleFilter {
private:
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    
    bool lastTimeinit=false;
    double lastTime;
    
    FILE *fBathy;
    FILE *fWalls;
    
    Eigen::Matrix<double, 2, PARTICLE_NUMBER> particles;
    Eigen::Matrix<double, 1, PARTICLE_NUMBER> weights;

    Eigen::Vector2d computeMean();
    Eigen::Matrix2d computeCovariance();

    Eigen::Vector2d computeWeightedMean();
    Eigen::Matrix2d computeWeightedCovariance();
    
    bool pointIsBetweenWalls(const Eigen::Vector2d &point);
    double getWallsRangeAt(const Eigen::Vector2d &point,const double &theta);
    double getBathyAt(const Eigen::Vector2d &point);

    Eigen::VectorXd drawSample(const Eigen::VectorXd &mean, const Eigen::MatrixXd &cov);
    std::vector<Eigen::VectorXd> drawSamples(const int &nbParticle, const Eigen::VectorXd &mean, const Eigen::MatrixXd &cov);
    
    double mahalanobisDistance(const Eigen::Vector2d &point, const Eigen::Matrix2d &shapeMatrix);
public:

    ParticleFilter();
    
    void init(const Eigen::Vector2d &initPos, const Eigen::Matrix2d &initPosCov);
    
    /**
     * 
     * @param filePath ABSOLUTE PATH!!
     */
    void setWallsFile(const std::string &filePath);
    
    /**
     * 
     * @param filePath ABSOLUTE PATH!!
     */
    void setBathymetricFile(const std::string &filePath);
    
    /**
     * 
     * @param moosTime
     * @param theta
     * @param u
     * @param uCov
     */
    void predict(const double &t, const Eigen::Vector2d &u, const Eigen::Matrix2d &uCov);

    /**
     * 
     * @param range
     * @param rangeVar
     */
    void update_range(const Eigen::Vector2d &emitterPos, const double &range, const double &rangeVar);

    /**
     * 
     * @param beam
     * @param beamVar
     */
    void update_sonar(const double &beam, const double &beamVar);

    /**
     * 
     * @param beam
     * @param beamVar
     */
    void update_echosounder(const double &beam, const double &beamVar);
    
    /**
     * Uses Low-variance resampling method
     */
    void resample();
};

#endif // PARTICLE_FILTER_H
