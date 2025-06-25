#ifndef BLACK_SCHOLES_HPP
#define BLACK_SCHOLES_HPP

#include <cmath>

namespace options {

struct OptionParams {
    double S;
    double K;
    double T;
    double r;
    double sigma;
};

struct BlackScholesResult {
    double call_price;
    double put_price;
};

inline double normalCDF(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

inline double normalPDF(double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2 * M_PI);
}

inline BlackScholesResult price(const OptionParams& params) {
    double d1 = (std::log(params.S / params.K) + (params.r + 0.5 * params.sigma * params.sigma) * params.T) / (params.sigma * std::sqrt(params.T));
    double d2 = d1 - params.sigma * std::sqrt(params.T);

    double call_price = params.S * normalCDF(d1) - params.K * std::exp(-params.r * params.T) * normalCDF(d2);
    double put_price = params.K * std::exp(-params.r * params.T) * normalCDF(-d2) - params.S * normalCDF(-d1);

    return {call_price, put_price};
}

}

#endif
