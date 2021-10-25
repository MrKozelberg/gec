//
// Created by mrk on 25.10.2021.
//

#ifndef GEC_AREA_H
#define GEC_AREA_H

#include <cmath>

/**
 * @brief   Parent class for classes of area parameterization
 *
 *          It is possible to work with different parameterization
 */
class ParentArea {
protected:
    static constexpr double earth_radius2 = 6370.0 * 6370.0; ///< km^2
public:
    ParentArea() = default;
    ~ParentArea() = default;
};

/**
 * @brief This class provides a function that calculates area of geographical cell
 */
class GeoArea : public ParentArea {
public:
    GeoArea() : ParentArea() {};

    /**
     * @brief This compute areas of geographical cells (like on a globe)
     *
     * It takes arguments in deg
     */
    static double area(double lat_1, double lat_2, double lon_1, double lon_2) {
        assert(lat_2 > lat_1 and lon_2 > lon_1);
        assert(lat_2 >=-90 and lat_2 <=90);
        assert(lat_1 >=-90 and lat_1 <=90);
        assert(lon_2 >=0 and lon_2 <=360);
        assert(lon_1 >=0 and lon_1 <=360);
        return earth_radius2 * (lon_2 - lon_1) / 180.0 * M_PI * (std::sin(lat_2 / 180.0 * M_PI) - std::sin(lat_1 / 180.0 * M_PI));
    }
    //    static double area(size_t n, size_t N, size_t m, size_t M, double d_lat, double d_lon) {
    //        double lat_n = -90.0 + d_lat * double(n);
    //        if (n != N - 1 and m != M - 1) {
    //            return fabs(earth_radius2 * M_PI / 180.0 * d_lon *
    //                        (sin(M_PI / 180.0 * (lat_n + d_lat)) - sin(M_PI / 180.0 * lat_n)));
    //        } else {
    //            if (m == M - 1) {
    //                return fabs(earth_radius2 * M_PI / 180.0 * (360.0 - double(m) * d_lon) *
    //                            (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
    //            } else {
    //                return fabs(earth_radius2 * M_PI / 180.0 * d_lon *
    //                            (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
    //            }
    //        }
    //    }
};

#endif //GEC_AREA_H
