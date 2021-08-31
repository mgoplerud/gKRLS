#ifndef calc_kernel
#define calc_kernel

double kern_gauss(Eigen::MatrixXd X_one, 
                  Eigen::MatrixXd X_two, 
                  double bandwidth);

#endif