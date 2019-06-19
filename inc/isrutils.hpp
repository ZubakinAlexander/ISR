#ifndef __ISRUTILS_H
#define __ISRUTILS_H
#include <functional>

double integral_radpol(double, double, double, int);

double integral_KuraevFadin(double, const std::function<double(double)>&,
                            double);

double kernel_KuraevFadin(double, double);

#endif
