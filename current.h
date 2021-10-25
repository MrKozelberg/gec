//
// Created by mrk on 25.10.2021.
//

#ifndef GEC_CURRENT_H
#define GEC_CURRENT_H

#include <cmath>

/**
 * @brief   Parent class for classes of parameterization of source currents
 */
class ParentCurrent {
protected:
    static constexpr double j_0 = 6.4e-9;
public:
    ParentCurrent() = default;
    double j[steps]{};
};

/**
 * @brief   Simplest parameterization of current
 *
 *          This class provides function that creates array of current value оn the considered grid
 * @tparam  Alt Height grid
 */
template<class Alt>
class StepCurrent : public ParentCurrent {
public:
    static double current_func(double z, ...) {
        return (z >= 5.0 and z <= 10.0) ? j_0 : 0.0;
    }

    StepCurrent(...) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func(a.altitude[i]);
        }
    }
};

/**
 * @brief Zero parameterization of current, it is needed for testing
 *
 * This class provides function that creates array of current value оn the considered grid
 * @tparam Alt Height grid
 */
template<class Alt>
class [[maybe_unused]] ZeroCurrent : public ParentCurrent {
public:
    static double current_func(...) {
        return 0.0;
    }

    ZeroCurrent(...) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func();
        }
    }
};

/**
 * @brief Simple parameterization of current
 *
 * This class provides function that creates array of current value оn the considered grid
 * @tparam Alt Height grid
 */
template<class Alt>
class [[maybe_unused]] SimpleGeoCurrent : public ParentCurrent {
public:
    static double current_func(double z, double lat, ...) {
        return (std::abs(lat) <= 5 and std::abs(z-7.5) <= 2.5) ? j_0 : 0.0;
    }

    explicit SimpleGeoCurrent(unsigned lat, ...) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func(a.altitude[i], lat);
        }
    }
};

/**
 * @brief Parameterization of current which we mainly work with
 * This class provides function that creates array of current value оn the considered grid
 * @tparam Alt Height grid
 */
template<class Alt>
class GeoCurrent : public ParentCurrent {
private:
    static constexpr double cape_0 = 1'000; // J/kg
public:
    /**
     * @brief current_func
     * @param z     in km
     * @param lat   in deg
     * @param lon   in deg
     * @param cbot  in m
     * @param ctop  in m
     * @param cape  in J/kg
     * @return
     */
    double current_func(double z, double lat, double lon, double cbot, double ctop, double cape) {
        if (cape >= cape_0 and z*1000>=cbot and z*1000<=ctop){
            return j_0;
        } else{
            return 0.0;
        }
    }

    GeoCurrent(unsigned lat, unsigned lon, double cbot, double ctop, double cape) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func(a.altitude[i], lat, lon, cbot, ctop, cape);
        }
    }
};

#endif //GEC_CURRENT_H
