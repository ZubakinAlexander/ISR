#ifndef __INTUTILS_H
#define __INTUTILS_H_
#include <functional>

double integral_qag(std::function<double(double)>&, double, double, double&);
double integral_qags(std::function<double(double)>&, double, double, double&);

#endif
