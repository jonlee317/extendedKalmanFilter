#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

// for debug
#include <fstream>
using namespace std;

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0187, 0,
			  0, 0.02;
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.1, 0, 0,
			       0, 0.1, 0,
              0, 0, 0.05;
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;
  Hj_ = MatrixXd(3, 4);

  /**
  TODO:
    * Finish initializing the FusionEKF.
  */
  VectorXd x_in = VectorXd(4);

	MatrixXd P = MatrixXd(4, 4);
	P << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

	MatrixXd F = MatrixXd(4, 4);
	F << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

	MatrixXd H = MatrixXd(2, 4);
	H << 1, 0, 0, 0,
			  0, 1, 0, 0;

	MatrixXd R = MatrixXd(2, 2);
  R << 0.015, 0,
			  0, 0.015;

	MatrixXd I = MatrixXd::Identity(2, 2);

	MatrixXd Q = MatrixXd(4, 4);

  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;

  ekf_.Init(x_in, P, F, H, R, Q);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      //VectorXd xloc = measurement_pack.raw_measurements_[0];

      double x_loc = measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]);
      double y_loc = measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]);
      double vx = measurement_pack.raw_measurements_[2]*cos(measurement_pack.raw_measurements_[1]);
      double vy = measurement_pack.raw_measurements_[2]*sin(measurement_pack.raw_measurements_[1]);

      ekf_.x_ << x_loc, y_loc, vx, vy;


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double x_loc = measurement_pack.raw_measurements_[0];
      double y_loc = measurement_pack.raw_measurements_[1];
      ekf_.x_ << x_loc, y_loc, ekf_.x_[2], ekf_.x_[3];

    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;

    ofstream myfile;
    myfile.open ("debug_meas.txt");
    myfile << previous_timestamp_;
    myfile.close();

    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */
   float dt = (measurement_pack.timestamp_ - previous_timestamp_)/ 1000000.0;
   float dt_2 = dt * dt;
   float dt_3 = dt_2 * dt;
   float dt_4 = dt_3 * dt;

   ofstream myfile;
   myfile.open ("debug_meas.txt");
   myfile << dt;
   myfile.close();

   ekf_.F_ = MatrixXd(4, 4);
   // TODO: find out what elapse time dt is
   //measurement_pack.timestamp_;

   ekf_.F_ << 1, 0, dt, 0,
              0, 1, 0, dt,
              0, 0, 1, 0,
              0, 0, 0, 1;

              //set the process covariance matrix Q
  float noise_ax = 8;
  float noise_ay = 8;
	ekf_.Q_ = MatrixXd(4, 4);

	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();
  previous_timestamp_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

     Tools tools;
     Hj_ = tools.CalculateJacobian(measurement_pack.raw_measurements_);
     ekf_.H_ = Hj_;
     ekf_.R_ = R_radar_;
     ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    //double x_loc = measurement_pack.raw_measurements_[0];
    //double y_loc = measurement_pack.raw_measurements_[1];
    //ekf_.x_ << x_loc, y_loc, ekf_.x_[2], ekf_.x_[3];
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
