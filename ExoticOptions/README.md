# Option Pricing and Hedging in C++

This repository contains a set of C++ implementations for pricing and hedging derivative securities, based on the assignment specifications shown in the project brief.
It was prepared as part of a computational finance / derivatives pricing assignment.  
The exact implementation details may vary depending on the numerical schemes used in each file.

## Overview

The project covers four main topics in Computational Finance:

1. **Binomial pricing of vanilla put options**
   - European put
   - American put

2. **Monte Carlo pricing of path-dependent options**
   - Lookback call option with floating strike

3. **Monte Carlo pricing of barrier options**
   - European barrier knock-out put

4. **Hedging cost simulation**
   - Stop-loss hedging for a European put
   - Delta hedging for a European put

The implementations are designed to compare pricing methods, exercise features, and hedging performance under different market assumptions.

## Project Structure

The repository will contain the expected following source files:

```text
LatticeEurPut.cpp
LatticeAmerPut.cpp
LookbackCallMC.cpp
BarrierKnockoutPut.cpp
StopLossPut.cpp
DeltaHedgingPut.cpp
```

### File descriptions

#### `LatticeEurPut.cpp`
Computes the price of a **European put option** using the **Binomial model**.

**Inputs**
- `S`: current underlying asset price
- `K`: strike price
- `r`: continuously compounded risk-free interest rate
- `T`: time to maturity
- `v`: volatility
- `N`: number of time steps

---

#### `LatticeAmerPut.cpp`
Computes the price of an **American put option** using the **Binomial model**.

**Inputs**
- `S`: current underlying asset price
- `K`: strike price
- `r`: continuously compounded risk-free interest rate
- `T`: time to maturity
- `v`: volatility
- `N`: number of time steps

---

#### `LookbackCallMC.cpp`
Computes the price of a **Lookback call option with floating strike** using the **Monte Carlo method**.

**Inputs**
- `S`: current underlying asset price
- `r`: continuously compounded risk-free interest rate
- `T`: time to maturity
- `v`: volatility
- `N`: number of time steps
- `num_sims`: number of Monte Carlo simulations

---

#### `BarrierKnockoutPut.cpp`
Computes the price of a **European barrier knock-out put option** using the **Monte Carlo method**.

**Inputs**
- `S`: current underlying asset price
- `Sb`: barrier level
- `K`: strike price
- `r`: continuously compounded risk-free interest rate
- `T`: time to maturity
- `v`: volatility
- `N`: number of time steps
- `num_sims`: number of Monte Carlo simulations

---

#### `StopLossPut.cpp`
Estimates the **cost of stop-loss hedging** for a **European put option** by simulation.

**Inputs**
- `S`: current underlying asset price
- `K`: strike price
- `μ`: annual rate of return of the underlying asset
- `r`: continuously compounded risk-free interest rate
- `T`: time to maturity
- `v`: volatility
- `N`: number of time steps
- `num_sims`: number of simulations

---

#### `DeltaHedgingPut.cpp`
Estimates the **cost of delta hedging** for a **European put option** by simulation.

**Inputs**
- `S`: current underlying asset price
- `K`: strike price
- `μ`: annual rate of return of the underlying asset
- `r`: continuously compounded risk-free interest rate
- `T`: time to maturity
- `v`: volatility
- `N`: number of time steps
- `num_sims`: number of simulations

---

## Methods Used

### 1. Binomial Model
Used for pricing European and American put options on a discrete time lattice.

Key idea was to:
- Build a recombining tree for the stock price
- Compute option payoffs at maturity
- Work backward through the tree using risk-neutral valuation
- For the American put, compare continuation value with early exercise value at each node

---

### 2. Monte Carlo Simulation
Used for pricing path-dependent and barrier options.

Key idea again:
- Simulate many stock price paths under a geometric Brownian motion framework
- Compute the payoff for each path
- Discount the average payoff back to present value

Applications in this project:
- **Lookback call with floating strike**: payoff depends on the minimum asset price along the path
- **Barrier knock-out put**: payoff becomes zero if the barrier is hit during the life of the option

---

### 3. Hedging Simulations
Used to compare practical hedging strategies for a European put.

- **Stop-loss hedging**: a simpler rule-based hedging approach
- **Delta hedging**: a more refined dynamic hedging strategy based on option sensitivity

The goal is to compare:
- hedging cost
- replication accuracy
- efficiency of each method

---

## Example Assignment Scenarios

### Binomial pricing comparison
Use:
- `S = 35`
- `K = 40`
- `r = 3%`
- `T = 9/12 = 0.75`
- `v = 25%`

Compare the European and American put prices.

Expected qualitative result:
- The **American put** should be worth **at least as much as** the European put
- In many cases, it is worth more because of the early exercise feature

---

### Lookback call vs at-the-money European call
Use:
- `S = 40`
- `r = 5%`
- `T = 6/12 = 0.5`
- `v = 30%`
- `N = 1000`
- `num_sims = 10000`

Compare the floating-strike lookback call with an at-the-money European call.

Expected qualitative result:
- The **lookback call** is typically more expensive because its payoff benefits from the best price observed along the path

---

### Barrier knock-out put vs standard European put
Use:
- `S = 40`
- `Sb = 55`
- `K = 30`
- `r = 5%`
- `T = 6/12 = 0.5`
- `v = 30%`
- `N = 1000`
- `num_sims = 15000`

Compare the knock-out put with a standard European put.

Expected qualitative result:
- The **barrier knock-out put** should be cheaper because the option can expire worthless if the barrier is breached before maturity

---

### Stop-loss hedging vs delta hedging
Use:
- `S = 27`
- `K = 30`
- `μ = 20%`
- `r = 5%`
- `T = 6/12 = 0.5`
- `v = 30%`
- `N = 80`
- `num_sims = 13000`

Compare both hedging methods in terms of cost and accuracy.

Expected qualitative result:
- **Delta hedging** is usually more accurate than stop-loss hedging
- **Stop-loss hedging** may be simpler, but often produces larger replication errors

---

## Build Instructions

You can compile each file separately with `g++`.

### Example
```bash
g++ -O2 -std=c++17 LatticeEurPut.cpp -o LatticeEurPut
g++ -O2 -std=c++17 LatticeAmerPut.cpp -o LatticeAmerPut
g++ -O2 -std=c++17 LookbackCallMC.cpp -o LookbackCallMC
g++ -O2 -std=c++17 BarrierKnockoutPut.cpp -o BarrierKnockoutPut
g++ -O2 -std=c++17 StopLossPut.cpp -o StopLossPut
g++ -O2 -std=c++17 DeltaHedgingPut.cpp -o DeltaHedgingPut
```

Run a program with:

```bash
./LatticeEurPut
```

If your code is written to read from standard input, you can enter parameters interactively or redirect them from a file.

---

## Suggested Improvements

Possible extensions for this project include:
- adding Black–Scholes closed-form benchmarks where applicable
- computing confidence intervals for Monte Carlo prices
- improving random number generation
- supporting command-line arguments instead of interactive input
- exporting results to CSV for analysis and plotting
- adding convergence studies with respect to `N` and `num_sims`

---

## Learning Goals

This project demonstrates my tryout on:
- implementation of lattice methods for derivatives
- Monte Carlo pricing of path-dependent contracts
- pricing impact of early exercise and barrier conditions
- comparison of hedging strategies in incomplete/discrete settings
- practical C++ programming for quantitative finance

---


