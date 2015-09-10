#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#define PARTICLE_NUMBER 500 // Number of particles

#include "../Eigen/Dense"

class ParticleFilter {
private:
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

    double mahalanobisDistance(const Eigen::Vector2d &point, const Eigen::Matrix2d &shapeMatrix);
public:

    ParticleFilter();
    
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
    void predict(const double &t, const double &theta, const double &u, const Eigen::Matrix2d &uCov);

    /**
     * 
     * @param range
     * @param rangeVar
     */
    void update_range(const double &range, const double &rangeVar);

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
};

#endif // PARTICLE_FILTER_H
