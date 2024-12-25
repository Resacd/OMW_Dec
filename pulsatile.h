#pragma once

#ifndef _PALSATILE_H_
#define _PALSATILE_H_

#include <cmath>

class Pulsatile {

public:

    double period = 0.9;

    inline double InletVelequation1(double x) {
        return (-8.579e+08 * pow(x, 8) + 4.059e+08 * pow(x, 7) - 7.817e+07 * pow(x, 6) +
            7.994e+06 * pow(x, 5) - 4.835e+05 * pow(x, 4) + 1.782e+04 * pow(x, 3) -
            321.6 * pow(x, 2) + 3.215 * x + 0.1423);
    }

    inline double InletVelequation2(double x) {
        return (-4.32e+06 * pow(x, 8) + 7.94e+06 * pow(x, 7) - 6.298e+06 * pow(x, 6) +
            2.82e+06 * pow(x, 5) - 7.809e+05 * pow(x, 4) + 1.37e+05 * pow(x, 3) -
            1.484e+04 * pow(x, 2) + 901.9 * x - 22.87);
    }

    inline double InletVelequation3(double x) {
        return (-5.123 * pow(x, 4) + 14.15 * pow(x, 3) - 14.21 * pow(x, 2) + 5.865 * x - 0.5696);
    }

};

#endif