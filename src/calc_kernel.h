#ifndef calc_kernel
#define calc_kernel

double kern_gauss(Eigen::VectorXd X_one, 
                  Eigen::VectorXd X_two, 
                  double bandwidth,
                  bool raw);

#endif