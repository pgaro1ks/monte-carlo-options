#ifndef GREEKS_HPP
#define GREEKS_HPP

#include "black_scholes.cpp"

namespace options {

struct Greeks {
    double delta;
    double gamma;
    double theta;
    double vega;
    double rho;
};

class GreeksCalculator {
public:
    static Greeks computeAnalytical(const OptionParams& params) {
        double d1 = (std::log(params.S / params.K) + (params.r + 0.5 * params.sigma * params.sigma) * params.T) / (params.sigma * std::sqrt(params.T));
        double d2 = d1 - params.sigma * std::sqrt(params.T);
        
        double Nd1 = normalCDF(d1);
        double Nd2 = normalCDF(d2);
        double nd1 = normalPDF(d1);
        
        double delta = Nd1;
        double gamma = nd1 / (params.S * params.sigma * std::sqrt(params.T));
        
        double theta_call = -((params.S * nd1 * params.sigma) / (2 * std::sqrt(params.T))) - params.r * params.K * std::exp(-params.r * params.T) * Nd2;
        double theta_put = -((params.S * nd1 * params.sigma) / (2 * std::sqrt(params.T))) + params.r * params.K * std::exp(-params.r * params.T) * normalCDF(-d2);
        
        double vega = params.S * std::sqrt(params.T) * nd1;
        double rho = params.K * params.T * std::exp(-params.r * params.T) * Nd2;
        
        return {delta, gamma, (theta_call + theta_put) / 2 / 365.0, vega / 100.0, rho / 100.0};
    }
    
    static Greeks computeFiniteDifference(const OptionParams& params, double h = 0.01) {
        OptionParams p_up = params;
        OptionParams p_down = params;
        OptionParams p_vol_up = params;
        
        p_up.S = params.S * (1.0 + h);
        p_down.S = params.S * (1.0 - h);
        p_vol_up.sigma = params.sigma * (1.0 + h);
        
        auto [call_up, _] = price(p_up);
        auto [call_down, __] = price(p_down);
        auto [call_vol_up, ___] = price(p_vol_up);
        
        auto [call_curr, ____] = price(params);
        
        double delta = (call_up - call_down) / (2 * params.S * h);
        double gamma = (call_up - 2 * call_curr + call_down) / (params.S * params.S * h * h);
        double vega = (call_vol_up - call_curr) / (params.sigma * h);
        
        return {delta, gamma, 0.0, vega / 100.0, 0.0};
    }
};

}

#endif
