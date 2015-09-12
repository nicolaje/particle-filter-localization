#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#define PARTICLE_NUMBER 500 // Number of particles

#include <random>
#include "../Eigen/Dense"
#include <octomap/octomap.h>
#include <octomap/OcTree.h>

class ParticleFilter {
private:
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    std::uniform_real_distribution<double> uniformDistribution;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    octomap::OcTree *bathyMap;
    
    double maxRange;
    
    bool lastTimeinit=false;
    double lastTime;
    
    std::vector<std::pair<double, double> > wallsPolygon;
    
    Eigen::Matrix<double, 2, PARTICLE_NUMBER> particles;
    double logWeights[PARTICLE_NUMBER];
    
    Eigen::Matrix<double, 2, PARTICLE_NUMBER> particlesSwap;
    double logWeightsSwap[PARTICLE_NUMBER];
    
    double uniforms[PARTICLE_NUMBER];
    
    void normalize();
    
    bool pointIsBetweenWalls(const Eigen::Vector2d &point);
    double getWallsRangeAt(const Eigen::Vector2d &point,const double &theta);
    double getBathyAt(const Eigen::Vector2d &point);

    Eigen::VectorXd drawSample(const Eigen::VectorXd &mean, const Eigen::MatrixXd &cov);
    std::vector<Eigen::VectorXd> drawSamples(const int &nbParticle, const Eigen::VectorXd &mean, const Eigen::MatrixXd &cov);
    
    double mahalanobisDistance(const Eigen::Vector2d &point, const Eigen::Matrix2d &shapeMatrix);
public:

    ParticleFilter();
    
    void init(const Eigen::Vector2d &initPos, const Eigen::Matrix2d &initPosCov);
    
    void setMaxRange(const double &range);
    
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
     * @param beamRange
     * @param beamRangeVar
     * @param beamAngle angle, from the EAST of the GLOBAL FRAME, in trigonometric direction, in DEGREES
     * @param z depth of the robot. Negative means robot is underwater.
     */
    void update_sonar_horizontal(const double &beamRange, const double &beamRangeVar, const double &beamAngle, const double &z);

    /**
     * 
     * @param beamRange
     * @param beamRangeVar
     * @param beamAngle angle of the beam in DEGREES. 0 Means the beam is parallel to the water, on the right of the robot. 180 Means the beam is parallel to the water, on the left of the robot. 90 Means the beams is at the nadir.
     * @param robotYaw yaw of the robot, in DEGREES, counting 0 from east, in a trigonometric fashion.
     * @param z
     */
    void update_sonar_vertical(const double &beamRange, const double &beamRangeVar, const double &beamAngle, const double &robotYaw, const double &z);
    
    /**
     * 
     * @param beam
     * @param beamVar
     * @param z depth of the robot. Negative means robot is underwater.
     */
    void update_echosounder(const double &beam, const double &beamVar, const double &z);
    
    /**
     * Uses Low-variance resampling method
     */
    void resample();
    
    Eigen::Vector2d computeMean();
    Eigen::Matrix2d computeCovariance();

    Eigen::Vector2d computeWeightedMean();
    Eigen::Matrix2d computeWeightedCovariance();
    
};

#endif // PARTICLE_FILTER_H
