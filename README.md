# Relative-Absolute-Orientation

### Description	

​	Use relative-absolute-orientation method to recover poses and triangulate object points.



### Class RAOrient

- RAOrient(const vector<vector<double>>& conjugate_points, const vector<vector<double>>& GCPs
  - input correspondence and their object points' coordinate
- void solve()
- vector<Vector3d> compute(const vector<vector<double>>& conjugate_points)
  - input: correspondences
  - output: object points



1. Instantiate a RAOrient object using known correspondence points and their object points;
2. Call solve() to solve the problem;
3. Call compute() to compute object points.





### Dependency

- Eigen