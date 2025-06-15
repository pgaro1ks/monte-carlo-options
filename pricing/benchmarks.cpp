#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <thread>
#include "option_pricing.hpp"

using namespace options;

void benchmarkBlackScholes() {
    OptionParams params{100.0, 100.0, 1.0, 0.05, 0.2};
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 100000; ++i) {
        blackScholesPrice(params);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Black-Scholes: " << duration.count() << "ms for 100k iterations\n";
}

void benchmarkMonteCarlo() {
    OptionParams params{100.0, 100.0, 1.0, 0.05, 0.2};
    MonteCarloConfig config{100000, std::thread::hardware_concurrency(), true};
    MonteCarloEngine engine(config);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    auto result = engine.price(params);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Monte Carlo (" << config.num_simulations << " sims, " << config.num_threads << " threads): " << duration.count() << "ms\n";
    std::cout << "  Call: " << std::fixed << std::setprecision(2) << result.call_price << "\n";
    std::cout << "  Put:  " << result.put_price << "\n";
}

void benchmarkGreeks() {
    OptionParams params{100.0, 100.0, 1.0, 0.05, 0.2};
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 10000; ++i) {
        GreeksCalculator::computeAnalytical(params);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Greeks (analytical): " << duration.count() << "ms for 10k iterations\n";
}

void benchmarkThreadScaling() {
    OptionParams params{100.0, 100.0, 1.0, 0.05, 0.2};
    
    std::vector<size_t> thread_counts{1, 2, 4, 8, std::thread::hardware_concurrency()};
    
    std::cout << "\nThread Scaling (1M simulations):\n";
    
    for (size_t threads : thread_counts) {
        MonteCarloConfig config{100000, threads, true};
        MonteCarloEngine engine(config);
        
        auto start = std::chrono::high_resolution_clock::now();
        auto result = engine.price(params);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "  " << threads << " threads: " << duration.count() << "ms (call: " << result.call_price << ")\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== Option Pricing Benchmarks ===\n\n";
    
    benchmarkBlackScholes();
    benchmarkGreeks();
    benchmarkMonteCarlo();
    benchmarkThreadScaling();
    
    return 0;
}
