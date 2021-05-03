#include "orient.h"
#include "include/Eigen/Dense"
#include <iostream>

#define PI 3.14159265358979323846264338327950288419716939937510




void RAOrient::solve() {
	inner_orient();
	relative_orient();
	absolute_orient();
}



vector<Vector3d> RAOrient::compute(const vector<vector<double>>& conjugates) {
	double oX = width - origX;
	double oY = height - origY;
	vector<vector<double>> conjugate_in_camera;
	for (size_t i = 0; i < conjugates.size(); i++) {
		vector<double> conjugate = conjugates[i];

		conjugate[0] = (conjugate[0] - oX) * pix_size;
		conjugate[1] = -(conjugate[1] - oY) * pix_size;
		conjugate[2] = (conjugate[2] - oX) * pix_size;
		conjugate[3] = -(conjugate[3] - oY) * pix_size;

		conjugate_in_camera.push_back(conjugate);
	}


	vector<Vector3d> points;
	forward_intersection(conjugate_in_camera, points);

	Matrix3d R = Matrix3d::Identity();
	R(0, 1) = epsz; R(0, 2) = -epsy;
	R(1, 0) = -epsz; R(1, 2) = epsx;
	R(2, 0) = epsy; R(2, 1) = -epsx;

	Vector3d shift = offset + gc;

	vector<Vector3d> res;
	for (size_t i = 0; i < points.size(); i++) {
		res.push_back(lambda * R * (points[i] - mc) + shift);
	}

	return res;
}


void RAOrient::pix2cam(vector<vector<double>>& _conjugate_point) {
	double oX = width - origX;
	double oY = height - origY;

	for (size_t i = 0; i < _conjugate_point.size(); i++) {
		vector<double> conjugate = _conjugate_point[i];
		_conjugate_point[i][0] = (conjugate[0] - oX) * pix_size;
		_conjugate_point[i][1] = -(conjugate[1] - oY) * pix_size;
		_conjugate_point[i][2] = (conjugate[2] - oX) * pix_size;
		_conjugate_point[i][3] = -(conjugate[3] - oY) * pix_size;
	}
}


vector<double> RAOrient::forward_intersection(
	const vector<vector<double>>& conjugates, 
	vector<Vector3d>& points) {
	Matrix<double, 3, 3> R1 = Matrix3d::Identity();
	Matrix<double, 3, 3> R2 = Matrix3d::Identity();
	R1(0, 1) = -k1; R1(0, 2) = -fi1;
	R1(1, 0) = k1;
	R1(2, 0) = fi1;

	R2(0, 1) = -k2; R2(0, 2) = -fi2;
	R2(1, 0) = k2; R2(1, 2) = -w2;
	R2(2, 0) = fi2; R2(2, 1) = w2;
	
	// parallax
	vector<double> q;
	for (int i = 0; i < conjugates.size(); i++) {
		vector<double> tmp = conjugates[i];
		Vector3d pt1(tmp[0], tmp[1], -f);
		Vector3d pt2(tmp[2], tmp[3], -f);

		Vector3d P1 = R1 * pt1;
		Vector3d P2 = R2 * pt2;

		double X1 = P1[0], Y1 = P1[1], Z1 = P1[2];
		double X2 = P2[0], Y2 = P2[1], Z2 = P2[2];

		double N1 = (base_line * Z2) / (X1 * Z2 - X2 * Z1);
		double N2 = (base_line * Z1) / (X1 * Z2 - X2 * Z1);

		double x = N1 * X1;
		double z = N1 * Z1;
		double y = 0.5 * (N1 * Y1 + N2 * Y2);
		q.push_back(N1 * Y1 - N2 * Y2);
		points.push_back(Vector3d(x, y, z));

	}

	return q;
}


void RAOrient::inner_orient() {
	double oX = width - origX;
	double oY = height - origY;

	cout << "inner orientation: " << endl;
	for (size_t i = 0; i < _conjugate_point.size(); i++) {
		vector<double> conjugate = _conjugate_point[i];
		_conjugate_point[i][0] = (conjugate[0] - oX) * pix_size;
		_conjugate_point[i][1] = -(conjugate[1] - oY) * pix_size;
		_conjugate_point[i][2] = (conjugate[2] - oX) * pix_size;
		_conjugate_point[i][3] = -(conjugate[3] - oY) * pix_size;
		cout << i + 1 << ": ";
		for (int j = 0; j < 4; j++) {
			cout << _conjugate_point[i][j] << ", ";
		}
		cout << endl;
	}

	for (size_t i = 0; i < _GCPs.size(); i++) {
		vector<double> conjugate = _GCPs[i];
		_GCPs[i][0] = (conjugate[0] - oX) * pix_size;
		_GCPs[i][1] = -(conjugate[1] - oY) * pix_size;
		_GCPs[i][2] = (conjugate[2] - oX) * pix_size;
		_GCPs[i][3] = -(conjugate[3] - oY) * pix_size;

		double t = _GCPs[i][4];
		_GCPs[i][4] = _GCPs[i][5];
		_GCPs[i][5] = t;
	}
}


double RAOrient::relative_orient(unsigned int max_iteration, double eps) {
	int num = _conjugate_point.size();
	
	Matrix<double, Dynamic, Dynamic> B;
	Matrix<double, Dynamic, Dynamic> L;
	B.resize(num, 5);
	L.resize(num, 1);
	for (int k = 0; k < max_iteration; k++) {
		Matrix<double, 3, 3> R1 = Matrix3d::Identity();
		Matrix<double, 3, 3> R2 = Matrix3d::Identity();
		R1(0, 1) = -k1; R1(0, 2) = -fi1;
		R1(1, 0) = k1;
		R1(2, 0) = fi1;

		R2(0, 1) = -k2; R2(0, 2) = -fi2;
		R2(1, 0) = k2; R2(1, 2) = -w2;
		R2(2, 0) = fi2; R2(2, 1) = w2;
		for (int i = 0; i < num; i++) {
			vector<double> tmp = _conjugate_point[i];
			Vector3d pt1(tmp[0], tmp[1], -f);
			Vector3d pt2(tmp[2], tmp[3], -f);

			Vector3d P1 = R1 * pt1;
			Vector3d P2 = R2 * pt2;

			double X1 = P1[0], Y1 = P1[1], Z1 = P1[2];
			double X2 = P2[0], Y2 = P2[1], Z2 = P2[2];

			B(i, 0) = - base_line * Y2 * X1;
			B(i, 1) = base_line * X1 * Z2;
			B(i, 2) = base_line * X2 * Y1;
			B(i, 3) = base_line * (Y1 * Y2 - f * Z1);
			B(i, 4) = - base_line * X2 * Z1;

			L(i, 0) = - base_line * (Y1 * Z2 - Y2 * Z1);
		}

		Matrix<double, Dynamic, Dynamic> x = (B.transpose() * B).inverse() * (B.transpose() * L);
		fi1 += x(0,0);
		k1 += x(1,0);
		fi2 += x(2,0);
		w2 += x(3,0);
		k2 += x(4,0);

		if (x.norm() < eps)
			break;
	}


	// compute parallax
	Matrix<double, 3, 3> R1 = Matrix3d::Identity();
	Matrix<double, 3, 3> R2 = Matrix3d::Identity();
	R1(0, 1) = -k1; R1(0, 2) = -fi1;
	R1(1, 0) = k1;
	R1(2, 0) = fi1;

	R2(0, 1) = -k2; R2(0, 2) = -fi2;
	R2(1, 0) = k2; R2(1, 2) = -w2;
	R2(2, 0) = fi2; R2(2, 1) = w2;
	double q = 0.0;
	for (int i = 0; i < num; i++) {
		vector<double> tmp = _conjugate_point[i];
		Vector3d pt1(tmp[0], tmp[1], -f);
		Vector3d pt2(tmp[2], tmp[3], -f);

		Vector3d P1 = R1 * pt1;
		Vector3d P2 = R2 * pt2;

		double X1 = P1[0], Y1 = P1[1], Z1 = P1[2];
		double X2 = P2[0], Y2 = P2[1], Z2 = P2[2];

		q += pow(f * (Y1 / Z1 - Y2 / Z2), 2);
	}
	
	return sqrt(q / num);
}



double RAOrient::absolute_orient(unsigned int max_iteration, double eps) {
	int num = _GCPs.size();
	vector<Vector3d> points;
	_q = forward_intersection(_GCPs, points);

	gc = Vector3d(0, 0, 0);
	mc = Vector3d(0, 0, 0);

	// calculate ground and model center
	for (int i = 0; i < points.size(); i++) {
		gc[0] += _GCPs[i][4];
		gc[1] += _GCPs[i][5];
		gc[2] += _GCPs[i][6];

		mc += points[i];
	}

	gc /= points.size();
	mc /= points.size();

	vector<Vector3d> centered_ground;
	vector<Vector3d> centered_model;
	// move origin to center
	for (int i = 0; i < points.size(); i++) {
		Vector3d pg(_GCPs[i][4] - gc[0], _GCPs[i][5] - gc[1], _GCPs[i][6] - gc[2]);
		centered_ground.push_back(pg);
		centered_model.push_back(points[i] - mc);

		//cout << i + 1 << ":\n";
		//cout << pg[0] << "," << pg[1] << "," << pg[2] << endl;
		//cout << centered_model[i][0] << "," << centered_model[i][1] << "," <<centered_model[i][2] << "\n";
		//cout << "-------------------------------------" << endl;
	}



	//// calculate theta and scale
	//double model_length = 0.0;
	//double ground_length = 0.0;
	//double theta = 0.0;
	//for (int i = 0; i < centered_model.size(); i++) {
	//	ground_length += centered_ground[i].norm() / 10;
	//	model_length += centered_model[i].norm() / 10;

	//	double gtheta = atan2f(centered_ground[i][1], centered_ground[i][0]);
	//	double mtheta = atan2f(centered_model[i][1], centered_model[i][0]);
	//	double theta_diff = gtheta - mtheta;
	//	if (theta_diff < 0)
	//		theta_diff += 2 * PI;
	//	cout << theta_diff << endl;
	//	theta += theta_diff;
	//}

	//theta /= centered_ground.size();
	//lambda = ground_length / model_length;
	//

	//// scale the model
	//AngleAxisd Rotation_vector(-0.05, Eigen::Vector3d(0, 0, 1));
	//Matrix3d rot = Rotation_vector.matrix();
	//for (int i = 0; i < centered_model.size(); i++) {
	//	centered_model[i] *= lambda;
	//	centered_model[i] = rot * centered_model[i];

	//	cout  << centered_ground[i][0] << "," << centered_ground[i][1] << "," << centered_ground[i][2] << endl;
	//	cout  << centered_model[i][0] << "," << centered_model[i][1] << "," << centered_model[i][2] << endl;
	//	cout << "-------------------------------------------\n";
	//}


	
	Matrix<double, Dynamic, Dynamic> B;
	Matrix<double, Dynamic, Dynamic> L;
	Matrix<double, Dynamic, Dynamic> x;
	B.resize(num * 3, 7);
	L.resize(num * 3, 1);

	offset = Vector3d(0, 0, 0);
	for (int k = 0; k < max_iteration; k++) {
		Matrix3d R = Matrix3d::Identity();
		R(0, 1) = epsz; R(0, 2) = -epsy;
		R(1, 0) = -epsz; R(1, 2) = epsx;
		R(2, 0) = epsy; R(2, 1) = -epsx;


		for (int i = 0; i < num; i++) {
			Vector3d ground_tilde = lambda * R * centered_model[i] + offset;

			Vector3d Pt(centered_ground[i][0], centered_ground[i][1], centered_ground[i][2]);

			B(i * 3, 0) = 1; B(i * 3, 1) = 0; B(i * 3, 2) = 0;
			B(i * 3 + 1, 0) = 0; B(i * 3 + 1, 1) = 1; B(i * 3 + 1, 2) = 0;
			B(i * 3 + 2, 0) = 0; B(i * 3 + 2, 1) = 0; B(i * 3 + 2, 2) = 1;

			B(i * 3, 3) = centered_model[i][0];
			B(i * 3 + 1, 3) = centered_model[i][1];
			B(i * 3 + 2, 3) = centered_model[i][2];

			B(i * 3, 4) = 0;
			B(i * 3 + 1, 5) = 0;
			B(i * 3 + 2, 6) = 0;

			B(i * 3, 5) = -centered_model[i][2];
			B(i * 3, 6) = centered_model[i][1];
			B(i * 3 + 1, 4) = centered_model[i][2];
			B(i * 3 + 1, 6) = -centered_model[i][0];
			B(i * 3 + 2, 4) = -centered_model[i][1];
			B(i * 3 + 2, 5) = centered_model[i][0];


			L(i * 3, 0) = (Pt - ground_tilde)[0];
			L(i * 3 + 1, 0) = (Pt - ground_tilde)[1];
			L(i * 3 + 2, 0) = (Pt - ground_tilde)[2];

			
		}
		
		x = (B.transpose() * B).inverse() * (B.transpose() * L);

		offset[0] += x(0, 0);
		offset[1] += x(1, 0);
		offset[2] += x(2, 0);
		lambda += x(3, 0);
		epsx += x(4, 0);
		epsy += x(5, 0);
		epsz += x(6, 0);
		if (x.norm() < eps)
			break;
	}

	auto v = B * x - L;
	double mse = pow(v.norm(), 2) / centered_ground.size();
	return mse;
}