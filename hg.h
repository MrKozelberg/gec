//
// Created by mrk on 25.10.2021.
//

#ifndef GEC_HG_H
#define GEC_HG_H

#include "constants.h"
#include <cmath>

/**
 * @brief (Parent Height Grid) Parent class for classes of height grids
 *
 * It is possible to create several parameterization of height grid
 */
class ParentHG {
public:
    double altitude[steps]{};
    ParentHG() = default;
    ~ParentHG() = default;
};

/**
 * @brief   Linear height grid with a break of a derivative
 *
 *          (Add a picture!)
 */
class LinHG: public ParentHG {
public:
    LinHG(): ParentHG() {
        for (size_t i = 0; i < steps; i++) {
            if (i <= n_1 * size_t(z_1)) {
                altitude[i] = double(i) / double(n_1);
            } else {
                altitude[i] = (double(i) - double(n_1) * z_1) * (z_max - z_1) / (n_2 - 1) + z_1;
            };
        }
    }
};

#endif //GEC_HG_H
