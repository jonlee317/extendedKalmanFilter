#include "kalman_filter.h"

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_*x_; //TODO:  WHere does u come from?  I guess we don't need it
  P_ = F_*P_*F_.transpose()+Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  MatrixXd PHt = P_*H_.transpose();
  VectorXd y = z - H_*x_;  // H*x is the measurement prediction of z
  MatrixXd S = H_*PHt + R_;
  MatrixXd K = PHt*S.inverse();

  x_ = x_ + K*y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I-K*H_)*P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  double rho = sqrt(x_[0]*x_[0]+x_[1]*x_[1]);
  double theta;
  double rho_dot;
  VectorXd z_pred = VectorXd(3);

  if (rho > 0) {
    theta = atan(x_[1]/x_[0]);
    rho_dot = ((x_[0]*x_[2]+x_[1]*x_[3])/rho);
  } else {
    return;
  }

  z_pred << rho, theta, rho_dot;

  MatrixXd PHt = P_*H_.transpose();
  VectorXd y = z - z_pred;  // z_pred is the measurement prediction of z
  MatrixXd S = H_*PHt + R_;
  MatrixXd K = PHt*S.inverse();

  x_ = x_ + K*y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I-K*H_)*P_;
}
