//
// Created by mrk on 25.10.2021.
//

#ifndef GEC_CONSTANTS_H
#define GEC_CONSTANTS_H

constexpr size_t lat_dim = 180;
constexpr size_t lon_dim = 360;
constexpr double delta_lat = 180.0 / double(lat_dim); ///< dimensions of latitude-longitude cell by latitude
constexpr double delta_lon = 360.0 / double(lon_dim); ///< dimensions of latitude-longitude cell by longitude
constexpr double z_1 = 21.0; ///< the break of height grid (integer)
constexpr double z_max = 70.0;
constexpr size_t n_1 = 11; ///< points per kilometer lower than z_1
constexpr size_t n_2 = 19; ///< points per upper atmosphere
constexpr size_t steps = n_1 * size_t(z_1) + n_2; ///< total number of points

#endif //GEC_CONSTANTS_H
