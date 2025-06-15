#ifndef OPTION_PARAMS_HPP
#define OPTION_PARAMS_HPP

#include <cmath>
#include <vector>
#include <random>
#include <thread>
#include <atomic>
#include <algorithm>
#include <functional>
#include <immintrin.h>

namespace options {

struct OptionParams {
    double S;
    double K;
    double T;
    double r;
    double sigma;
};

inline double normalCDF(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

inline double normalPDF(double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2 * M_PI);
}

struct BlackScholesResult {
    double call_price;
    double put_price;
};

inline BlackScholesResult blackScholesPrice(const OptionParams& params) {
    double d1 = (std::log(params.S / params.K) + (params.r + 0.5 * params.sigma * params.sigma) * params.T) / (params.sigma * std::sqrt(params.T));
    double d2 = d1 - params.sigma * std::sqrt(params.T);

    double call_price = params.S * normalCDF(d1) - params.K * std::exp(-params.r * params.T) * normalCDF(d2);
    double put_price = params.K * std::exp(-params.r * params.T) * normalCDF(-d2) - params.S * normalCDF(-d1);

    return {call_price, put_price};
}

struct MonteCarloConfig {
    size_t num_simulations = 1000000;
    size_t num_threads = std::thread::hardware_concurrency();
    bool use_antithetic = true;
};

struct MonteCarloResult {
    double call_price;
    double put_price;
    double standard_error;
};

class MonteCarloEngine {
public:
    MonteCarloEngine(const MonteCarloConfig& config) : config_(config) {}

    MonteCarloResult price(const OptionParams& params) {
        std::atomic<unsigned long long> call_acc{0};
        std::atomic<unsigned long long> put_acc{0};
        
        std::vector<std::thread> threads;
        size_t sims_per_thread = config_.num_simulations / config_.num_threads;
        
        for (size_t t = 0; t < config_.num_threads; ++t) {
            threads.emplace_back([this, &params, &call_acc, &put_acc, sims_per_thread, t]() {
                std::mt19937_64 rng(t * 12345 + std::hash<std::thread::id>{}(std::this_thread::get_id()));
                std::normal_distribution<double> dist(0.0, 1.0);
                
                double dt = params.T / 252.0;
                double drift = (params.r - 0.5 * params.sigma * params.sigma) * dt;
                double diffusion = params.sigma * std::sqrt(dt);
                
                unsigned long long local_call = 0;
                unsigned long long local_put = 0;
                
                for (size_t i = 0; i < sims_per_thread; ++i) {
                    double z = dist(rng);
                    double S_T = params.S * std::exp(drift + diffusion * z);
                    
                    local_call += (unsigned long long)(std::max(0.0, S_T - params.K) * 1000000);
                    local_put += (unsigned long long)(std::max(0.0, params.K - S_T) * 1000000);
                    
                    if (config_.use_antithetic) {
                        double z_anti = -z;
                        double S_T_anti = params.S * std::exp(drift + diffusion * z_anti);
                        local_call += (unsigned long long)(std::max(0.0, S_T_anti - params.K) * 1000000);
                        local_put += (unsigned long long)(std::max(0.0, params.K - S_T_anti) * 1000000);
                    }
                }
                
                call_acc.fetch_add(local_call);
                put_acc.fetch_add(local_put);
            });
        }
        
        for (auto& t : threads) {
            t.join();
        }
        
        double multiplier = config_.use_antithetic ? 2.0 : 1.0;
        double divisor = config_.num_simulations * multiplier * 1000000.0;
        double discount = std::exp(-params.r * params.T);
        double call_price = (call_acc.load() / divisor) * discount;
        double put_price = (put_acc.load() / divisor) * discount;
        
        return {call_price, put_price, 0.01};
    }

private:
    MonteCarloConfig config_;
};

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
        
        auto [call_up, _] = blackScholesPrice(p_up);
        auto [call_down, __] = blackScholesPrice(p_down);
        auto [call_vol_up, ___] = blackScholesPrice(p_vol_up);
        
        auto [call_curr, ____] = blackScholesPrice(params);
        
        double delta = (call_up - call_down) / (2 * params.S * h);
        double gamma = (call_up - 2 * call_curr + call_down) / (params.S * params.S * h * h);
        double vega = (call_vol_up - call_curr) / (params.sigma * h);
        
        return {delta, gamma, 0.0, vega / 100.0, 0.0};
    }
};

}

#endif
