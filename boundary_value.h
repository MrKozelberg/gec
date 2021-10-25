//
// Created by mrk on 25.10.2021.
//

#ifndef GEC_BOUNDARY_VALUE_H
#define GEC_BOUNDARY_VALUE_H

#include <cmath>

/**
 * @brief   Parent class for classes that provide parameterization of ionospheric sources
 *
 *          This is not a boundary value, it is something like a boundary value
 */
class ParentBoundValue {
protected:
    static constexpr double phi_0 = 15.0; ///<    [kV]
    static constexpr double theta_0 = 20.0; ///<  [deg]
public:
    ParentBoundValue() = default;
};
/**
 * @brief   Class for zero parameterization of phi_s
 *          We mainly use this parameterization
 */
class [[maybe_unused]] ZeroPhiS : public ParentBoundValue {
public:
    static double phi_s(...) {
        return 0.0;
    }
};

/**
 * @brief Volland's parameterization of ionospheric sources
 *
 * @warning Mb it doesn't work, it was writen off some article
 */
class [[maybe_unused]] VollandPhiS : public ParentBoundValue {
public:

    /**
     * @param theta it is a longitude in deg
     * @return Some k from the article
     */
    static double k(double theta) {
        try {
            if (theta < theta_0) {
                throw theta;
            }
        } catch (double i) {
            std::cout << "k(theta) has a wrong argument" << std::endl;
            exit(-2);
        }
        if (theta < 30.0) {
            return 1.0;
        } else if (theta >= 30.0 && theta < 90.0) {
            return (1 + 2 * sin(3.0 * theta)) / 2;
        }
        return 0.0;
    }

    /**
     * @param theta It is a geomagnetic longitude
     * @param lambda it is a latitude in deg
     * @return The value of ionospheric sources potential
     */
    static double phi_s(double theta, double lambda) {
        return phi_0 * sin(lambda) * ((theta < theta_0) ?
                                      sin(theta) / sin(theta_0) :
                                      k(theta) * pow(sin(theta_0) / sin(theta), 4));
    }
};

#endif //GEC_BOUNDARY_VALUE_H
