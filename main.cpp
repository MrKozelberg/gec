#include <cmath>
#include <ctime>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

#include "conductivity.h"
#include "geomagnetic.h"
#include "cnpy/cnpy.h"
#include "constants.h"
#include "hg.h"
#include "area.h"
#include "current.h"
#include "boundary_value.h"

template<typename T>
T sqr(T x) {
    return x * x;
}

template<typename T>
T cub(T x) {
    return x * x * x;
}

/**
 * @brief This is a container of column data
 *
 * Atmosphere is divorced into columns (add the scheme!)
 */
struct Column {
    double area{};
    double altitude[steps]{}; ///< array of height points
    double sigma[steps]{}; ///< conductivity
    double j_s[steps]{}; ///< source current
    double phi_s{}; ///< additive IP from different ionospheric sources

    // test
    double ctop{};
    double cbot{};
    double cape{};
    bool isThereSource{};
};

/**
 * @brief This class provides possibility to calculate IP and the electric potential vs altitude
 *
 * @warning The electric potential vs altitude is temporarily unavailable
 */
class GECModel {
protected:
    std::vector<Column> model;
private:
    double IP = 0.0;
    bool isIPCalculated = false;

    /**
     * @brief IP calculator
     *
     * Calculation of integrals with the help of Simpson's rule
     * On each step integrands is approximated by parabola with coefficients A, B, C
     */
    double calc_IP_Simp() {
        double up = 0.0;
        double down = 0.0;
        for (size_t i = 0; i < model.capacity(); ++i) {
            double int_curr_by_cond = 0.0;
            double int_inv_cond = 0.0;
            double x1 = 0, x2 = 0, x3 = 0, y1 = 0, y2 = 0, y3 = 0, denom = 0;
            double A = 0, B = 0, C = 0;
            for (size_t k = 0; k + 2 < steps; k+=2){
                x1 = model[i].altitude[k], x2 = model[i].altitude[k+1], x3 = model[i].altitude[k+2];
                y1 = 1 / model[i].sigma[k], y2 = 1 / model[i].sigma[k+1], y3 = 1 / model[i].sigma[k+2];
                denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
                B = (sqr(x3) * (y1 - y2) + sqr(x2) * (y3 - y1) + sqr(x1) * (y2 - y3)) / denom;
                C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
                int_inv_cond += A * (cub(x3) - cub(x1)) / 3 + B * (sqr(x3) - sqr(x1)) / 2 + C * (x3 - x1);
                y1 = model[i].j_s[k] / model[i].sigma[k], y2 = model[i].j_s[k+1] / model[i].sigma[k+1], y3 = model[i].j_s[k] / model[i].sigma[k];
                denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
                B = (sqr(x3) * (y1 - y2) + sqr(x2) * (y3 - y1) + sqr(x1) * (y2 - y3)) / denom;
                C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
                int_curr_by_cond += A * (cub(x3) - cub(x1)) / 3 + B * (sqr(x3) - sqr(x1)) / 2 + C * (x3 - x1);
            }
            up += model[i].area * (int_curr_by_cond - model[i].phi_s) / int_inv_cond;
            down += model[i].area / int_inv_cond;
        }
        return up / down;
    }

    double calc_IP_trap() {
        double up = 0.0;
        double down = 0.0;
        double h = 0.0;

        /**
         * @warning .capacity() shows how much memory has been reserved for this vector and .size() shows how much memory is being used
         */
        for (size_t i = 0; i < model.capacity(); ++i) {
            volatile double int_curr_by_cond = 0.0;
            volatile double int_inv_cond = 0.0;
            for (size_t q = 1; q < steps; ++q) {
                h = model[i].altitude[q] - model[i].altitude[q - 1];
                int_curr_by_cond += (model[i].j_s[q] / model[i].sigma[q] + model[i].j_s[q - 1] / model[i].sigma[q - 1]) * h / 2.0;
                int_inv_cond += (1 / model[i].sigma[q] + 1 / model[i].sigma[q - 1]) * h / 2.0;
            }
            up = up + model[i].area * (int_curr_by_cond - model[i].phi_s) / int_inv_cond;
            down = down + model[i].area / int_inv_cond;
        }
        return up / down;
    }

    /// test
    double calc_IP_test(){
        double sigma_0 = 6.0e-14;
        double H_0 = 6.0;
        double j_0 = 6.4e-9;
        double cape_0 = 1e3;
        double sum_up = 0.0, sum_down = 0.0;
        for(size_t i = 0; i < model.capacity(); ++i) {
            // here the area already takes into account the factor alpha
            if (model[i].cape >= cape_0 and model[i].isThereSource) {
                sum_up += model[i].area * j_0 * H_0 / sigma_0 * (std::exp(-model[i].cbot / 1000 / H_0) - std::exp(-model[i].ctop / 1000 / H_0));
            }
            sum_down += model[i].area;
        }
        return sum_up / sum_down;
    }

public:
    double phi[steps-1]{}; ///< potential vs height (from 1/n_1 km to ...)
    GECModel() {
        assert(steps - n_2 % 2 != 0);
    }
    double get_IP() {
        if (not isIPCalculated) {
            IP = calc_IP_Simp();
            isIPCalculated = true;
        }
        return IP;
    }
    /**
     * @brief get_phi_trap  computes an array of potential values
     * @param i             the number of the considered cell
     * @param filename      output
     */
    void get_phi_trap(size_t i, const std::string& filename) {
        double h = 0.0;
        double c2 = 0.0;
        double c1 = 0.0;
        double i1 = 0.0;
        double i2 = 0.0;

        /// calculation of c1, c2
        for (size_t k = 1; k < steps; ++k) {
            h = model[i].altitude[k] - model[i].altitude[k - 1];
            c1 += (1 / model[i].sigma[k] + 1 / model[i].sigma[k - 1]) * h / 2;
            c2 += (model[i].j_s[k] / model[i].sigma[k] + model[i].j_s[k - 1] / model[i].sigma[k - 1]) * h / 2.0;
        }
        //        std::cout << c1 << "\t\t" << c2 << std::endl;
        std::ofstream fout(filename);
        if (!fout.is_open()) {
            std::cout << "Impossible to find a file" << "\t" << filename << std::endl;
            exit(-1);
        }

        /// calculation of i1, i2
        for (size_t p = 1; p < steps; ++p) {
            h = model[i].altitude[p] - model[i].altitude[p - 1];
            i1 += (1 / model[i].sigma[p] + 1 / model[i].sigma[p - 1]) * h / 2;
            i2 += (model[i].j_s[p] / model[i].sigma[p] + model[i].j_s[p - 1] / model[i].sigma[p - 1]) * h / 2.0;
            //            std::cout << p << "\t" << i1 << "\t\t\t" << i2 << "\t\t" << i2 - i1/c1*(c2-get_IP()) << std::endl;
            phi[p-1] = i2 - i1/c1*(c2-get_IP());
        }
        for (size_t j = 0; j < steps-1; ++j){
            fout << model[i].altitude[j+1] << " " << phi[j] << std::endl;
        }
        fout.close();
    }
};

/**
 * @brief The ParentLatLonGrid class    Parent class for latitude-longitude partitioning classes, using 1x1 deg grid
 *                                      Works with NPZ-files taken from colleagues
 *
 * @warning                             py_array[n,m,k] is cpp_array[n*180*360 + m*360 + k]
 *                                      CHECK datatypes of npy-files!
 *                                      '<f8' is an equivalent of double
 *                                      '<f4' is an equivalent of float
 */
class ParentLatLonGrid : public GECModel {
protected:
    static constexpr double earth_radius2 = 6371.0 * 6371.0; ///< Earth radius in km^2
public:
    ParentLatLonGrid() : GECModel() {};

    /**
     * @brief lat_arg   calculates arguments by latitude for some functions
     *                  This function converts the cell number by latitude into the real latitude value of the considered cell
     * @param n         the number of the considered cell by latitude
     * @return          latitude in gedrees
     */
    static double lat_arg(size_t n) {
        return -90.0 + 0.5 + delta_lat * double(n);
    }

    /**
     * @brief lon_arg   calculates arguments by longitude for some functions
     * @param m         the number of the considered cell by longitude
     * @return          longitude in gedrees
     */
    static double lon_arg(size_t m) {
        return 0.5 + delta_lon * double(m);
    }

    /**
     * @brief lon_arg_m calculates arguments by geomagnetic longitude for some functions
     * @param lon_arg_  geodesic longitude in degrees
     * @param lat_arg_  geodesic latitude in degrees
     * @param year_     the considered year in the following format:
     *                  {year}.{the day number from the year beginning by the number of days per this year}
     * @return          geomagnetic longitude in degrees
     */
    static double lon_arg_m(double lon_arg_, double lat_arg_/*, double year_*/) {
        double lat_m, lon_m, alt_m;
        gdz_to_mag(/*year_,*/ lat_arg_, lon_arg_, 10.0, lat_m, lon_m, alt_m);
        return lon_m;
    }

    /**
     * @brief lat_arg_m calculates arguments by geomagnetic latitude for some functions
     * @param lon_arg_  geodesic longitude in gedrees
     * @param lat_arg_  geodesic latitude in gedrees
     * @param year_     the considered year in the following format:
     *                  {year}.{the day number from the year beginning by the number of days per this year}
     * @return          geomargetic latitude in degrees
     */
    static double lat_arg_m(double lon_arg_, double lat_arg_/*, double year_*/) {
        double lat_m, lon_m, alt_m;
        gdz_to_mag(/*year_,*/ lat_arg_, lon_arg_, 10.0, lat_m, lon_m, alt_m);
        return lat_m;
    }
};

template<class Alt, class PhiS, class Cond>
/**
 * @brief   The DataProc class      Child class for latitude-longitude partitioning considering date and time
 *                                  Its main purpose is to process data
 * @warning Cond should be Sigmd<Alt>!
 */
class DataProc : public ParentLatLonGrid {
private:
    static constexpr double coef = 1.97; ///< some proportionality coefficient that defines what area occupied by sources
    double set_year = 0.0;
    /**
     * @brief setYear   Initializes year in Fortran code if this has not been done before
     * @param year      Considered year
     *                  This is necessary for the transformation of coordinates
     */
    void setYear(double year){
        if (year ==! set_year) {
            init_dtd(year);
            set_year = year;
        }
    }
public:
    static_assert(std::is_same<Alt, typename Cond::AltClass>::value,  "Error, Alt and AltClass must be same types");
    /**
     * @brief DataProc
     * @param year  {year}.{the day number from the year beginning by the number of days per this year}
     * @param hour  from 0 to 48
     * @param input filename
     */
    explicit DataProc(double year, size_t hour, double xi, std::string input) : ParentLatLonGrid() {
        setYear(year);
        // input DATA
        cnpy::npz_t data = cnpy::npz_load(std::move(input));
        cnpy::NpyArray cape_arr = data["cape"];
        cnpy::NpyArray cbot_arr = data["cbot"];
        cnpy::NpyArray ctop_arr = data["ctop"];
        /**
         * @brief The portion of the area of the grid column occupied by GEC generators
         *
         * It equals the ratio of non-cumulative "rainc" to "pw"
         */
        cnpy::NpyArray alpha_arr = data["alpha"];

        auto *cape = cape_arr.data<float>(); ///< mb it is double, check it
        auto *cbot = cbot_arr.data<double>();
        auto *ctop = ctop_arr.data<double>();
        auto *alpha = alpha_arr.data<double>();

        model.reserve(2 * lat_dim * lon_dim);

        /**
         * This variable serves as the lower boundary of the cells when calculating their areas;
         * when finding the values of functions from coordinates, use a special function for the argument
         */
        double lat_n = -90.0;
        double S;
        for (size_t n = 0; n < lat_dim; ++n) {
            double lon_n = 0.0;
            for (size_t m = 0; m < lon_dim; ++m) {
                Alt alt;
                Cond cond(lat_arg_m(lon_arg(m), lat_arg(n)/*, year*/),xi);
                ZeroCurrent<Alt> zero_j;
                GeoCurrent<Alt> geo_j(n, m, cbot[hour * 180 * 360 + n * 360 + m],
                        ctop[hour * 180 * 360 + n * 360 + m],
                        cape[hour * 180 * 360 + n * 360 + m]);

                if (n == lat_dim - 1 and m == lon_dim - 1) {
                    S = GeoArea::area(lat_n, 90.0, lon_n, 360.0);
                } else if (n == lat_dim - 1) {
                    S = GeoArea::area(lat_n, 90.0, lon_n, lon_n + delta_lon);
                } else if (m == lon_dim - 1) {
                    S = GeoArea::area(lat_n, lat_n + delta_lat, lon_n,360.0);
                } else {
                    S = GeoArea::area(lat_n, lat_n + delta_lat, lon_n, lon_n + delta_lon);
                }

                /// Part without sources
                size_t n1 = n * 2 * lon_dim + 2 * m;
                model[n1].area = S * (1 - coef * alpha[hour * 360 * 180 + n * 360 + m]);
                std::copy(std::begin(alt.altitude), std::end(alt.altitude), std::begin(model[n1].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma), std::begin(model[n1].sigma));
                std::copy(std::begin(zero_j.j), std::end(zero_j.j), std::begin(model[n1].j_s));
                model[n1].phi_s = PhiS::phi_s(lat_arg(n), lon_arg(m));

                // test
                model[n1].cbot = cbot[hour * 180 * 360 + n * 360 + m];
                model[n1].ctop = ctop[hour * 180 * 360 + n * 360 + m];
                model[n1].cape = cape[hour * 180 * 360 + n * 360 + m];\
                model[n1].isThereSource = false;

                /// Part with sources
                size_t n2 = n * 2 * lon_dim + 2 * m + 1;
                model[n2].area = S * coef * alpha[hour * 360 * 180 + n * 360 + m];
                std::copy(std::begin(alt.altitude), std::end(alt.altitude), std::begin(model[n2].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma), std::begin(model[n2].sigma));
                std::copy(std::begin(geo_j.j), std::end(geo_j.j), std::begin(model[n2].j_s));
                model[n2].phi_s = PhiS::phi_s(lat_arg(n), lon_arg(m));

                // test
                model[n2].cbot = cbot[hour * 180 * 360 + n * 360 + m];
                model[n2].ctop = ctop[hour * 180 * 360 + n * 360 + m];
                model[n2].cape = cape[hour * 180 * 360 + n * 360 + m];
                model[n2].isThereSource = true;

                lon_n += delta_lon;
            }
            lat_n += delta_lat;
        }
    }
};

template<class Alt, class PhiS>
/**
 * @brief The DataProcTest1 class   It includes exponential conductivity, step-current (the area of sources defines from data)
 *                                  and it has spherical geometry
 */
class DataProcTest1 : public ParentLatLonGrid {
private:
    static constexpr double coef = 1; ///< some proportionality coefficient that defines what area occupied by sources
public:
    /**
     * @brief DataProcTest1 the constructor
     * @param hour          from 0 to 48
     * @param path          output
     */
    explicit DataProcTest1(size_t hour, std::string path) : ParentLatLonGrid() {
        // input DATA
        cnpy::npz_t data = cnpy::npz_load(std::move(path));
        cnpy::NpyArray cape_arr = data["cape"];
        cnpy::NpyArray cbot_arr = data["cbot"];
        cnpy::NpyArray ctop_arr = data["ctop"];
        /**
         * @brief alpha_arr The portion of the area of the grid column occupied by GEC generators
         *                  It equals the ratio of non-cumulative "rainc" to "pw"
         */
        cnpy::NpyArray alpha_arr = data["alpha"];

        auto *cape = cape_arr.data<double>(); ///< as it is processed data (not float)
        auto *cbot = cbot_arr.data<double>();
        auto *ctop = ctop_arr.data<double>();
        auto *alpha = alpha_arr.data<double>();

        model.reserve(2 * lat_dim * lon_dim);
        /**
         * @brief lat_n serves as the lower boundary of the cells when calculating their areas
         */
        double lat_n = -90.0;
        double S;
        for (size_t n = 0; n < lat_dim; ++n) {
            double lon_n = 0.0;
            for (size_t m = 0; m < lon_dim; ++m) {
                Alt alt;
                ExpCond<Alt> cond{};
                ZeroCurrent<Alt> zero_j;
                GeoCurrent<Alt> geo_j(n, m, cbot[hour * 180 * 360 + n * 360 + m],
                        ctop[hour * 180 * 360 + n * 360 + m],
                        cape[hour * 180 * 360 + n * 360 + m]);

                //                if(n * 360 + m == 90 * 360 + 90){
                //                    std::cout << cbot[hour * 180 * 360 + n * 360 + m] << "\t" << ctop[hour * 180 * 360 + n * 360 + m] << "\t" << cape[hour * 180 * 360 + n * 360 + m] << std::endl;
                //                }

                if (n == lat_dim - 1 and m == lon_dim - 1) {
                    S = GeoArea::area(lat_n, 90.0, lon_n, 360.0);
                } else if (n == lat_dim - 1) {
                    S = GeoArea::area(lat_n, 90.0, lon_n, lon_n + delta_lon);
                } else if (m == lon_dim - 1) {
                    S = GeoArea::area(lat_n, lat_n + delta_lat, lon_n,360.0);
                } else {
                    S = GeoArea::area(lat_n, lat_n + delta_lat, lon_n, lon_n + delta_lon);
                }

                /// Part without sources
                size_t n1 = n * 2 * lon_dim + 2 * m;
                model[n1].area = S * (1 - coef * alpha[hour * 360 * 180 + n * 360 + m]);
                std::copy(std::begin(alt.altitude), std::end(alt.altitude), std::begin(model[n1].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma), std::begin(model[n1].sigma));
                std::copy(std::begin(zero_j.j), std::end(zero_j.j), std::begin(model[n1].j_s));
                model[n1].phi_s = PhiS::phi_s(lat_arg(n), lon_arg(m));

                /// Part with sources
                size_t n2 = n * 2 * lon_dim + 2 * m + 1;
                model[n2].area = S * coef * alpha[hour * 360 * 180 + n * 360 + m];
                std::copy(std::begin(alt.altitude), std::end(alt.altitude), std::begin(model[n2].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma), std::begin(model[n2].sigma));
                std::copy(std::begin(geo_j.j), std::end(geo_j.j), std::begin(model[n2].j_s));
                model[n2].phi_s = PhiS::phi_s(lat_arg(n), lon_arg(m));

                //                if(alpha[hour * 360 * 180 + n * 360 + m] ==! 0){
                //                    std::cout << n << "\t" << m << std::endl;
                //                }

                lon_n += delta_lon;
            }
            lat_n += delta_lat;
        }
        //        std::cout << alpha[hour * 360 * 180 + 90 * 360 + 90] << std::endl;
    }
};

/**
 * @brief   This is very simplistic GEC model, using just to test the programme
 *
 *          1x1, exp-sigma, simple-geo-j, zero-phi_s
 */
template<class Alt>
/**
 * @brief The Test1 class
 */
class Test1 : public ParentLatLonGrid {
public:
    Test1() : ParentLatLonGrid() {
        model.reserve(lat_dim * lon_dim);
        double lat_n = -90.0;
        for (size_t n = 0; n < lat_dim; ++n) {
            //std::cout << n << "\t" << lat_n << "\t" << lat_arg(n, delta_lat) << "\n";
            for (size_t m = 0; m < lon_dim; ++m) {
                Alt alt{};
                ExpCond<Alt> cond{};
                StepCurrent<Alt> j_test{};
                ZeroCurrent<Alt> j_zero{};

                size_t n1 = n * lon_dim + m;
                model[n1].area = 1.0;
                std::copy(std::begin(alt.altitude), std::end(alt.altitude), std::begin(model[n1].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma), std::begin(model[n1].sigma));
                model[n1].phi_s = ZeroPhiS::phi_s();

                if (std::fabs(lat_arg(n)) < 5.0) {
                    std::copy(std::begin(j_test.j), std::end(j_test.j), std::begin(model[n1].j_s));
                } else {
                    std::copy(std::begin(j_zero.j), std::end(j_zero.j), std::begin(model[n1].j_s));
                }
            }
            lat_n += delta_lat;
        }
        //        std::ofstream basicOfstream("test/test1.txt");
        //        if (!basicOfstream.is_open()) {
        //            std::cout << "Impossible to find a file" << std::endl;
        //            exit(-1);
        //        }
        //        for (size_t k = 0; k < steps; ++k){
        //            basicOfstream << model[90*360 + 90].altitude[k] << "\t" << model[90*360 + 90].j_s[k] / model[90*360 + 90].sigma[k] << "\n";
        //        }
        //        basicOfstream.close();
    }
};

template<class Alt>
class Test2 : public ParentLatLonGrid {
public:
    Test2() : ParentLatLonGrid() {
        model.reserve(lat_dim * lon_dim);
        /*volatile*/ double lat_n = -90.0;
        for (size_t n = 0; n < lat_dim; ++n) {
            /*volatile*/ double lon_n = 0.0;
            //std::cout << n << "\t" << lat_n << "\t" << lat_arg(n, delta_lat) << "\n";
            for (size_t m = 0; m < lon_dim; ++m) {
                Alt alt{};
                ExpCond<Alt> cond{};
                StepCurrent<Alt> j_test{};
                ZeroCurrent<Alt> j_zero{};

                size_t n1 = n * lon_dim + m;
                if (n == lat_dim - 1 and m == lon_dim - 1) {
                    model[n1].area = GeoArea::area(lat_n, 90.0, lon_n, 360.0);
                } else if (n == lat_dim - 1) {
                    model[n1].area = GeoArea::area(lat_n, 90.0, lon_n, lon_n + delta_lon);
                } else if (m == lon_dim - 1) {
                    model[n1].area = GeoArea::area(lat_n, lat_n + delta_lat, lon_n,360.0);
                } else {
                    model[n1].area = GeoArea::area(lat_n, lat_n + delta_lat, lon_n, lon_n + delta_lon);
                }
                std::copy(std::begin(alt.altitude), std::end(alt.altitude), std::begin(model[n1].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma), std::begin(model[n1].sigma));
                model[n1].phi_s = ZeroPhiS::phi_s();

                if (std::fabs(lat_arg(n)) < 5.0) {
                    std::copy(std::begin(j_test.j), std::end(j_test.j), std::begin(model[n1].j_s));
                } else {
                    std::copy(std::begin(j_zero.j), std::end(j_zero.j), std::begin(model[n1].j_s));
                }
                lon_n += delta_lon;
            }
            lat_n += delta_lat;
        }
    }
};

template<class Alt>
/**
 * @brief The Test3 class with angles
 */
class Test3 : public ParentLatLonGrid {
public:
    Test3() : ParentLatLonGrid() {
        model.reserve(lat_dim * lon_dim);
        double lat_n = -90.0;
        for (size_t n = 0; n < lat_dim; ++n) {
            double lon_n = 0.0;
            //            std::cout << lat_arg(n) / 180.0 * M_PI << "\n";
            for (size_t m = 0; m < lon_dim; ++m) {
                Alt alt{};
                ExpCosCond<Alt> cond(lat_arg(n) / 180.0 * M_PI);
                StepCurrent<Alt> j_test{};
                ZeroCurrent<Alt> j_zero{};

                size_t n1 = n * lon_dim + m;
                if (n == lat_dim - 1 and m == lon_dim - 1) {
                    model[n1].area = GeoArea::area(lat_n, 90.0, lon_n, 360.0);
                } else if (n == lat_dim - 1) {
                    model[n1].area = GeoArea::area(lat_n, 90.0, lon_n, lon_n + delta_lon);
                } else if (m == lon_dim - 1) {
                    model[n1].area = GeoArea::area(lat_n, lat_n + delta_lat, lon_n,360.0);
                } else {
                    model[n1].area = GeoArea::area(lat_n, lat_n + delta_lat, lon_n, lon_n + delta_lon);
                }
                std::copy(std::begin(alt.altitude), std::end(alt.altitude), std::begin(model[n1].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma), std::begin(model[n1].sigma));
                model[n1].phi_s = ZeroPhiS::phi_s();

                if (std::fabs(lat_arg(n)) < 5.0) {
                    std::copy(std::begin(j_test.j), std::end(j_test.j), std::begin(model[n1].j_s));
                } else {
                    std::copy(std::begin(j_zero.j), std::end(j_zero.j), std::begin(model[n1].j_s));
                }
                lon_n += delta_lon;
            }
            lat_n += delta_lat;
        }
    }
};

int main(/*int argc, char* argv[]*/) {
/*    assert(argc == 2);
    std::string input(argv[1]);

    std::string file_out(input,9,14);
    file_out = "ipne/IP" + file_out + ".txt";

    std::ofstream fout(file_out);
    if (!fout.is_open()) {
        std::cout << "Imposible to find a file" << "\t" << file_out << std::endl;
    }

    for (int hour = 0; hour < 49; hour++) {
        DataProc<LinHG, ZeroPhiS, Conductivity<LinHG>> m(2015.9, hour, 0.5, input);
        fout << hour << "\t" << m.get_IP() << std::endl;
    }

    fout.close();*/

    std::ofstream fout("/home/mrk/gec/ip/IP-2016-FULL-NEW.txt");
    if (!fout.is_open()) {
        std::cout << "Imposible to find a file" << std::endl;
    }
    std::string input = "/home/mrk/gec/data/DATA-2016-FULL-NEW.npz";
    for (size_t hour = 0; hour < 25; hour++) {
        DataProc<LinHG, ZeroPhiS, Conductivity<LinHG>> m(2015.9, hour, 0.5, input);
        fout << hour << "\t" << m.get_IP() << std::endl;
    }
    fout.close();

//        DataProc<LinHG, ZeroPhiS, ExpCond<LinHG>> m(2015.9, 18, 0.5, "data/DATA-2015-12-31-00-NEW.npz");
//        std::cout << m.get_IP() << std::endl;
    //    m.get_phi_trap(90 * 2 * lon_dim + 2 * 90, "plots/90_90_0.txt");
    //    m.get_phi_trap(90 * 2 * lon_dim + 2 * 90 + 1, "plots/90_90_1.txt");

    /**
     * Test 1, the analytic answer is 8736.8, while the programme gives 8904.91 because of numerical integration
     */
    //    Test1<LinHG> m{};
    //    std::cout << m.get_IP() << std::endl;

    /**
     * Test 2, spherical geometry
     */
    //    Test2<LinHG> m{};
    //    std::cout << m.get_IP() << std::endl;

    /**
     * Test 3, angles
     */
    //    Test3<LinHG> m{};
    //    std::cout << m.get_IP() << std::endl;

    /**
     * DataProc Testing
     */
    //    DataProcTest1<LinHG,ZeroPhiS> m1(0, "data/test-01.npz");
    //    std::cout << m1.get_IP() << std::endl;
    //    DataProcTest1<LinHG,ZeroPhiS> m2(0, "data/test-02.npz");
    //    std::cout << m2.get_IP() << std::endl;
    //    DataProcTest1<LinHG,ZeroPhiS> m3(0, "data/test-03.npz");
    //    std::cout << m3.get_IP() << std::endl;

    //    m1.get_phi_trap(90 * 2 * lon_dim + 2 * 90, "plots/90_90_0.txt");
    //    m1.get_phi_trap(90 * 2 * lon_dim + 2 * 90 + 1, "plots/90_90_1.txt");
    return EXIT_SUCCESS;
}
