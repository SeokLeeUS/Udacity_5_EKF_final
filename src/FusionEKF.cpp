#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
* Constructor.
*/
FusionEKF::FusionEKF() {

	is_initialized_ = false;

	previous_timestamp_ = 0;

	// initializing matrices
	R_laser_ = MatrixXd(2, 2);
	R_radar_ = MatrixXd(3, 3);
	H_laser_ = MatrixXd(2, 4);
	Hj_ = MatrixXd(3, 4);


	Hj_ << 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0;
	//measurement covariance matrix - laser
	R_laser_ << 0.0225, 0,
		0, 0.0225;

	//measurement covariance matrix - radar
	R_radar_ << 0.09, 0, 0,
		0, 0.0009, 0,
		0, 0, 0.09;

	/**
	* TODO: Finish initializing the FusionEKF.
	* TODO: Set the process and measurement noises
	*/

	// process matrix initialization
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
		0, 1, 0, 1,
		0, 0, 1, 0,
		0, 0, 0, 1;


	// initialize variables and matrices (x, F, H_laser,H_jacobian,P, etc.)
	//


	H_laser_ << 1, 0, 0, 0,
		0, 1, 0, 0;

	ekf_.x_ = VectorXd(4);

	ekf_.P_ = MatrixXd(4, 4);

	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;

	ekf_.Q_ = MatrixXd(4, 4);

	ekf_.Q_ << 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
	/**
	* Initialization
	*/
	if (!is_initialized_) {
		/**
		* TODO: Initialize the state ekf_.x_ with the first measurement.
		* TODO: Create the covariance matrix.
		* You'll need to convert radar from polar to cartesian coordinates.
		*/

		// first measurement
		cout << "EKF: " << endl;
		ekf_.x_ = VectorXd(4);
		//ekf_.x_ << 1, 1, 1, 1; why set to 1 for initial measurement? shouldn't it be zero when it's not initialized?
		ekf_.x_ << 0, 0, 0, 0;
		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

			// TODO: Convert radar from polar to cartesian coordinates 
			//         and initialize state.
			float rho_ = measurement_pack.raw_measurements_[0];
			float phi_ = measurement_pack.raw_measurements_[1];
			float rho_dot_ = measurement_pack.raw_measurements_[2];
			/*
			ekf_.x_ <<  rho_*cos(phi_),
			rho_*sin(phi_),
			rho_dot_*cos(phi_),
			rho_dot_*sin(phi_);
			*/

			ekf_.x_ << rho_*cos(phi_),
				rho_*sin(phi_),
				0,
				0;

			// initialization for radar slee194
			cout << "Radar_initialized: " << is_initialized_ << endl;
			//ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
			ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, FusionEKF::Hj_, FusionEKF::R_radar_, ekf_.Q_);
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			// TODO: Initialize state.
			ekf_.x_ << measurement_pack.raw_measurements_[0],
				measurement_pack.raw_measurements_[1],
				0,
				0;

			cout << "Laser_initialized: " << is_initialized_ << endl;
			//ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
			ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, FusionEKF::H_laser_, FusionEKF::R_laser_, ekf_.Q_);
		}

		// done initializing, no need to predict or update
		previous_timestamp_ = measurement_pack.timestamp_;
		is_initialized_ = true;
		return;
	}

	/**
	* Prediction
	*/

	/**
	* TODO: Update the state transition matrix F according to the new elapsed time.
	* Time is measured in seconds.
	* TODO: Update the process noise covariance matrix.
	* Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	*/
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;
	float dt_2 = dt*dt;
	float dt_3 = dt_2*dt;
	float dt_4 = dt_3*dt;

	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;
	// slee194 set the acceleration noise components for process noise     
	float noise_ax = 9;
	float noise_ay = 9;

	// slee194 set the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
		0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
		dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
		0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

	cout << "which sensor before predict?(0:laser,1:radar)" << measurement_pack.sensor_type_ << endl;
	cout << "before prediction, check initialization done:" << is_initialized_ << endl;
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && ekf_.R_.rows() == 2)
	{
		ekf_.R_ = R_radar_;
	}
	else
	{
		ekf_.R_ = R_laser_;
	}

	if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && ekf_.H_.rows() == 3)
	{
		ekf_.H_ = H_laser_;
	}
	else
	{
		ekf_.H_ = Hj_;
	}


	ekf_.Predict();
	cout << "after prediction P_:\n" << ekf_.P_ << "\n" << "size of P_:\n " << "rows: " << ekf_.P_.rows() << " cols: " << ekf_.P_.cols() << endl;
	/**
	* Update
	*/

	/**
	* TODO:
	* - Use the sensor type to perform the update step.
	* - Update the state and covariance matrices.
	*/

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

		// TODO: Radar updates

		//slee194
		cout << "Radar update starts:\t" << "R_radar_row: " << ekf_.R_.rows() << " R_radar_column:" << ekf_.R_.cols() << endl;

		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
		cout << "Radar update success!!" << endl;

	}
	else {
		// TODO: Laser updates
		//slee194
		cout << "Laser update starts:\t" << "R_radar_row: " << ekf_.R_.rows() << " R_radar_column:" << ekf_.R_.cols() << endl;
		ekf_.Update(measurement_pack.raw_measurements_);

	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
