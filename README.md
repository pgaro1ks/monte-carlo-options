# Monte Carlo Option Pricing Engine

High-performance C++ library for options pricing using Monte Carlo simulation with Black-Scholes analytical pricing and Greeks calculation.

## Features

- Black-Scholes analytical pricing
- Monte Carlo simulation with multi-threading
- Antithetic variates for variance reduction
- Greeks calculation (Delta, Gamma, Theta, Vega, Rho)
- SIMD optimization support
- Thread scaling benchmarks

## Build

```bash
mkdir build && cd build
cmake ..
make
```

## Run Benchmarks

```bash
./benchmarks
```

## Usage

```cpp
#include "pricing/option_pricing.hpp"

using namespace options;

OptionParams params{100.0, 100.0, 1.0, 0.05, 0.2};

auto [call, put] = blackScholesPrice(params);

MonteCarloConfig config{100000, 4, true};
MonteCarloEngine engine(config);
auto result = engine.price(params);

auto greeks = GreeksCalculator::computeAnalytical(params);
```

## Performance

- Black-Scholes: ~1ms for 100k iterations
- Monte Carlo: ~50ms for 1M simulations (4 threads)
- Greeks: ~2ms for 10k iterations
