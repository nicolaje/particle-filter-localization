#include "vibes.h"
#include "simulator.h"
#include "ParticleFilter.h"
#include "../Eigen/Dense"
#include <unistd.h>
#include <cmath>

using namespace Eigen;
using namespace std;

int main(int argc, char** argv) {

    // Set Vibes Params
    vibes::beginDrawing();
    vibes::newFigure("test");
    vibes::setFigureProperties("test", vibesParams("x", 0, "y", 0, "width", 1000, "height", 1000));
    vibes::axisLimits(-10, 140, -10, 140);

    // Intantiare the simulator xith the map
    Simulator simu("map.txt", "traj.txt");
    simu.setSonarFreq(5); // 5hz
    simu.setNbOutliers(0); // 5 % outliers

    // Intantiate Localizer object
    ParticleFilter pf;
    Pose &p0 = simu.currentPose();

    Vector2d p(p0.x, p0.y);
    Matrix2d pCov = Matrix2d::Identity();
    
    pf.init(p, pCov);
    cout << "pos init: "<<endl<<p<<endl;
    /*localizer.setInitialPosition(p0.x, p0.y, p0.t);
    localizer.setSpeedNoise(0.02);
    localizer.setHeadingNoise(2*M_PI/180.);
    localizer.setBufferSize(20);
    localizer.setNOutliers(5);*/
    simu.drawMap();
    vibes::newGroup("greenlines", "green");
    vibes::newGroup("redlines", "red");

    vibes::newGroup("particles", vibesParams("color", "red"));

    vibes::newGroup("AUV", vibesParams("color", "black"));
    while (simu.nextStep() != -1) {

        // Current point of the trajectory
        Pose &p = simu.currentPose();

        // Draw the current situation
        //vibes::clearFigure("test");

        vibes::clearGroup("AUV");
        vibes::drawAUV(p.x, p.y, p.theta * 180 / M_PI, 1, vibesParams("group", "AUV"));
        Eigen::Vector2d u(p.speed, p.theta*180/M_PI);
        Eigen::Matrix2d uCov = Matrix2d::Identity();
        uCov(0,0)=1;//pow(0.2/3,2);
        uCov(1,1)=pow(6./3,2);
    
        pf.predict(p.t, u, uCov);
        //localizer.predict(p.speed, p.theta, p.t);

        double d = simu.genSonarDist(p.x, p.y, p.theta, p.t);
        
        Vector2d mean = pf.computeMean();
        
        Matrix2d cov=pf.computeCovariance();
        
        vibes::clearGroup("particles");
        
        Eigen::Matrix<double, 2, PARTICLE_NUMBER> particles=pf.getParticles();
        for(unsigned int i=0;i<PARTICLE_NUMBER;i++)
        {
            vibes::drawCircle(particles.col(i)[0],particles.col(i)[1],0.01,vibesParams("group","particles"));
        }
        
        vibes::drawConfidenceEllipse(mean[0],mean[1],cov(0,0),cov(0,1),cov(1,1),3,vibesParams("group","particles", "color","red"));
        if (d >= 0) { // the measurment is valid
            double alpha = simu.getSonarAngle(p.t) + p.theta; // sonar beam angle
            cout << "beam angle: " << alpha << endl;
            pf.update_walls(d,1,alpha,1,p.theta,1);
            pf.resample();
            //Interval ialpha = Interval(alpha).inflate(2*M_PI/180.0);
            //Interval irho = Interval(d).inflate(0.2);
            //localizer.update(irho, ialpha, p.t);
        }

        // Draw Results
        // Current position with corrections
        // Pour Thomas, quand la boite est vide on la remplace par
        // la valeur de la boite inertielle pure.
        /*if (localizer.X_cur.is_empty()){
            localizer.X_cur[0] = localizer.x_inertial;
            localizer.X_cur[1] = localizer.y_inertial;
            localizer.pos[0] = localizer.X_cur;
        }
        const IntervalVector& tmp = localizer.X_cur;

        vibes::drawBox( tmp[0].lb(), tmp[0].ub(), tmp[1].lb(), tmp[1].ub(), "g" );

        // Current postion without corrections (inertial only)
        Interval &x_cur = localizer.x_inertial;
        Interval &y_cur = localizer.y_inertial;
        vibes::drawBox(x_cur.lb(), x_cur.ub(), y_cur.lb(), y_cur.ub(), "r");*/
        usleep(4000);
    }
    return 0;
}
