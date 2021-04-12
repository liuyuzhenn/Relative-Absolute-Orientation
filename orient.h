#pragma once
#include <vector>
#include "include/Eigen/Dense"
using namespace std;
using namespace Eigen;

class RAOrient
{
public:
	/* 
	conjugate_points: conjugate points in 2 photos (xleft, yleft, xright, yright)
	GCPs: ground control points
	*/
	RAOrient(const vector<vector<double>>& conjugate_points, const vector<vector<double>>& GCPs) :
		_conjugate_point(conjugate_points), _GCPs(GCPs){}

	void solve();

	// for computing unknown points
	vector<Vector3d> compute(const vector<vector<double>>& conjugate_points);

private:
	void inner_orient();
	double relative_orient(unsigned int max_iteration = 20, double eps = 1e-8);
	double absolute_orient(unsigned int max_iteration = 20, double eps = 1e-8);
	vector<double> forward_intersection(const vector<vector<double>>& conjugate_points, vector<Vector3d>& points);


	// camera parameters
	double pix_size = 3.900000000001569e-6;
	double width = 14592;
	double height = 25728;
	double origX = 7295.5;
	double origY = 12863.5;

	double f = 0.092;
	double base_line = 500;
	vector<vector<double>> _conjugate_point, _GCPs;// xl,yl,xr,yr,X,Y,Z

	// relative element
	double fi1 = 0, k1 = 0, fi2 = 0, w2 = 0, k2 = 0;

	// absolute element
	Vector3d offset;
	double lambda = 1;
	double epsx = 0, epsy = 0, epsz = 0;

	// absolute orientation parallax
	vector<double> _q;
	
	Vector3d gc, mc;
};

