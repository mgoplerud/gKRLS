#ifndef calc_kernel
#define calc_kernel

double kern_gauss(const Eigen::VectorXd& X_one,
                  const Eigen::VectorXd& X_two,
                  double bandwidth,
                  bool raw);

#endif