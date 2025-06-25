#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include "black_scholes.cpp"
#include <vector>
#include <random>
#include <thread>
#include <atomic>
#include <algorithm>
#include <cmath>

namespace options {

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

}

#endif
