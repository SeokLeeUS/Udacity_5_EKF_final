#include "kalman_filter.h"
#include <iostream>
#include "tools.h"
Tools tools;
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
* Please note that the Eigen library does not initialize
*   VectorXd or MatrixXd objects with zeros upon creation.
*/
// constructor slee194
KalmanFilter::KalmanFilter() {}
//destructor slee194
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
	* TODO: predict the state
	*/
	cout << "x: " << x_ << endl;
	cout << "F: " << F_ << endl;
	x_ = F_*x_; //predicted state matrix slee194
	cout << "predicted F: " << x_ << endl;
	MatrixXd Ft = F_.transpose();
	cout << "Ft: " << Ft << endl;
	cout << "Q: " << Q_ << endl;
	P_ = F_*P_*Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	/**
	* TODO: update the state by using Kalman Filter equations
	*/
	//slee194
	cout << "size of x, row: " << x_.rows() << " cols: " << x_.cols() << endl;
	cout << "size of H, row: " << H_.rows() << " cols: " << H_.cols() << endl;
	VectorXd z_pred = H_*x_;
	cout << "z_pred: " << z_pred << endl;
	VectorXd y = z - z_pred;
	cout << "y_laser: " << y << endl;
	MatrixXd Ht = H_.transpose();
	cout << "Ht: " << Ht << endl;
	MatrixXd S = H_ * P_ * Ht + R_;
	cout << "S: " << S << endl;
	MatrixXd Si = S.inverse();
	cout << "Si: " << Si << endl;
	MatrixXd PHt = P_*Ht;
	cout << "PHt: " << PHt << endl;
	MatrixXd K = PHt*Si;
	cout << "K: " << K << endl;

	// new estimate

	x_ = x_ + (K*y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K*H_)*P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	* TODO: update the state by using Extended Kalman Filter equations
	*/

	MatrixXd Hj_ = tools.CalculateJacobian(x_);
	const float PI_F = 3.14159265358979f;
	cout << "Hj_:" << Hj_ << endl;
	VectorXd z_pred(3);
	float px = x_[0];
	float py = x_[1];
	float vx = x_[2];
	float vy = x_[3];
	float rho = sqrt(px*px + py*py);
	float phi = atan2(py, px);
	float rho_dot = (px*vx + py*vy) / rho;
	z_pred << rho, phi, rho_dot;
	cout << "z_pred_radar:" << z_pred << endl;

	VectorXd y = z - z_pred;

	if (y[1] >PI_F)
	{
		y[1] = y[1] - (2 * PI_F);
	}
	else if (y[1] <-PI_F)
	{
		y[1] = y[1] + (2 * PI_F);
	}




	cout << "y_radar:" << y << endl;

	MatrixXd Hj_t = Hj_.transpose();

	cout << "Hj_t_radar:" << Hj_t << endl;

	cout << "R_radar_row:" << R_.rows() << " R_radar_column:" << R_.cols() << endl;

	MatrixXd S = Hj_ * P_ * Hj_t + R_;

	cout << "S_radar:" << S << endl;
	MatrixXd Si = S.inverse();
	cout << "Si_radar:" << Si << endl;
	MatrixXd PHt = P_ * Hj_t;
	cout << "PHt_radar:" << PHt << endl;
	MatrixXd K = PHt * Si;
	//new estimate
	cout << "K:" << K << endl;

	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj_) * P_;
}
