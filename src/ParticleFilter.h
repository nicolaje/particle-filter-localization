#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#define PARTICLE_NUMBER 500 // Number of particles

#include <random>
#include "../Eigen/Dense"
#include <octomap/octomap.h>
#include <octomap/OcTree.h>
#include <algorithm>
#include <iostream>
#include <float.h>

    enum sampling_methods{MULTINOMIAL,RESIDUAL,LOW_VARIANCE};
class ParticleFilter {
private:
    double wSlow;
    double wFast;
    double aSlow;
    double aFast;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    std::uniform_real_distribution<double> uniformDistribution;
    std::bernoulli_distribution bernoulliDistribution;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    octomap::OcTree *bathyMap;

    bool isInit;

    double maxRange;

    bool lastTimeinit;
    double lastTime;

#ifndef WALL_MAP
#define WALL_MAP

struct Wall {
        double w[4];

        friend std::istream& operator>>(std::istream& str, Wall& data) {
            std::string line;
            Wall tmp;
            if (std::getline(str, line)) {
                std::stringstream iss(line);
                try {
                    iss >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3];
                    std::swap(data, tmp);
                } catch (std::exception& e) {
                    str.setstate(std::ios::failbit);
                }
            }
            return str;
        }

        double & operator[](std::size_t idx) {
            return w[idx];
        }
    };

    typedef std::vector<Wall> Walls;
#endif
    
    inline double CalcDet(double &ax, double& ay, double& bx, double &by){
    return ax*by - ay*bx;
}
 
    inline void CalcDistanceDirSegment(double& d,double& phi,
                                               double mx, double my, double theta,
                                               double ax, double ay, double bx, double by)
{      // Distance directionnelle du point m au segment [a,b].
    double ma_x=ax-mx;
    double ma_y=ay-my;
    double mb_x=bx-mx;
    double mb_y=by-my;
    double ab_x=bx-ax;
    double ab_y=by-ay;
    double ux=cos(theta);
    double uy=sin(theta);
    double z1=CalcDet(ma_x,ma_y,ux,uy);
    double z2=CalcDet(ux,uy,mb_x,mb_y);
    double z3=CalcDet(ma_x,ma_y,ab_x,ab_y);
    double z4=CalcDet(ux,uy,ab_x,ab_y);
    double z5=std::min(z1,std::min(z2,z3));
    double d1=z3/z4;
    d= (z5 < 0) ? 1000 :d1;
    phi=atan2(-ab_x,ab_y); //phi is the angle of the normal vector of [a,b]
}

    int nbScanline;
    // Number of scanline updates, reseted every "nbScanline" scanline
    int scanLineUpdates;
    
    Walls walls;

    Eigen::Matrix<double, 2, PARTICLE_NUMBER> particles;
    double weights[PARTICLE_NUMBER];

    Eigen::Matrix<double, 2, PARTICLE_NUMBER> particlesSwap;
    double weightsSwap[PARTICLE_NUMBER];

    Eigen::Matrix<double, 2, PARTICLE_NUMBER> particlesInertial;
    
    double uniforms[PARTICLE_NUMBER];

    void normalize();

    bool pointIsBetweenWalls(const Eigen::Vector2d &point);
    double getWallsRangeAt(const Eigen::Vector2d &point, const double &theta);
    double getBathyAt(const Eigen::Vector2d &point);

    Eigen::VectorXd drawSample(const Eigen::VectorXd &mean, const Eigen::MatrixXd &cov);
    std::vector<Eigen::VectorXd> drawSamples(const int &nbParticle, const Eigen::VectorXd &mean, const Eigen::MatrixXd &cov);

    double mahalanobisDistance(const Eigen::Vector2d &point, const Eigen::Matrix2d &shapeMatrix);
public:

    ParticleFilter();

    
    int sampling_method;
    /**
     * Resample every X scanline
     * @param nbScanLine
     */
    void setResampleEvery(const int &nbScanLine);

    void setResampleMethod(const int &method);
    
    void setAlphas(const double &alphaSlow,const double &alphaFast);
    
    Eigen::Matrix<double, 2, PARTICLE_NUMBER> &getParticles();
    
    
    void init(const Eigen::Vector2d &initPos, const Eigen::Matrix2d &initPosCov);

    bool isInitialized();

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

    void update_walls(const double &rho, const double &rhoVar, const double &alpha, const double &alphaVar);

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
     * 
     * @param n NORTHING
     * @param nVar
     * @param e EASTING
     * @param eVar
     */
    void update_GPS(const double &n, const double &nVar, const double &e, const double &eVar);

    /**
     * Uses Multinomial resampling method
     */
    void resample();

    void resampleResidual();
    
    /**
     * 
     */
    void resampleLowVariance();
    
    Eigen::Vector2d computeMean();
    Eigen::Matrix2d computeCovariance();

    Eigen::Vector2d computeWeightedMean();
    Eigen::Matrix2d computeWeightedCovariance();

};



inline double Det2(double &ax, double& ay, double& bx, double &by){
    return ax*by - ay*bx;
}

inline void DistanceDirSegment2(double& d,double& phi,
                                               double mx, double my, double theta,
                                               double ax, double ay, double bx, double by)
{      // Distance directionnelle du point m au segment [a,b].
    double ma_x=ax-mx;
    double ma_y=ay-my;
    double mb_x=bx-mx;
    double mb_y=by-my;
    double ab_x=bx-ax;
    double ab_y=by-ay;
    double ux=cos(theta);
    double uy=sin(theta);
    double z1=Det2(ma_x,ma_y,ux,uy);
    double z2=Det2(ux,uy,mb_x,mb_y);
    double z3=Det2(ma_x,ma_y,ab_x,ab_y);
    double z4=Det2(ux,uy,ab_x,ab_y);
    double z5=std::min(z1,std::min(z2,z3));
    double d1=z3/z4;
    d= (z5 < 0) ? DBL_MAX :d1;
    phi=atan2(-ab_x,ab_y); //phi is the angle of the normal vector of [a,b]
}
#endif // PARTICLE_FILTER_H
