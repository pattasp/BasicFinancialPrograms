#Basic EU Call/Put Options formulas from BSM PDE

import numpy as np
from scipy.stats import norm

def black_scholes(S, K, T, r, sigma):
    """
    European call/put prices (no dividends).

    S: spot price
    K: strike
    T: time to maturity (in years)
    r: risk-free rate (cont. comp.)
    sigma: volatility
    """
    if sigma <= 0 or T <= 0:
        raise ValueError("sigma and T must be > 0")

    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    C = S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)           # Call
    P = K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)          # Put
    print(f"\n--- Black–Scholes (No Dividends) ---")
    print(f"Spot Price (S): {S:,.2f}")
    print(f"Strike (K):     {K:,.2f}")
    print(f"Time (T):       {T:.3f} years")
    print(f"r:              {r:.4f}")
    print(f"sigma:          {sigma:.4f}")
    print(f"Call Price (C): {C:,.4f}")
    print(f"Put Price (P):  {P:,.4f}")
    print("------------------------------------\n")
    
    return C, P

def black_scholes_div(S,K,T,r,sigma,q): #Here 'q' represents the dividend yield
      """
      European call/put prices (with dividends).

      S: spot price
      K: strike
      T: time to maturity (in years)
      r: risk-free rate (cont. comp.)
      sigma: volatility
      q :dividend yield
      """
      if sigma <= 0 or T <= 0:
          raise ValueError("sigma and T must be > 0")

      d1 = (np.log(S / K) + (r -q + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
      d2 = d1 - sigma * np.sqrt(T)

      C = S *  np.exp(-q * T)*norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)           # Call
      P = K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)          # Put
      
      print(f"\n--- Black–Scholes with Dividends ---")
      print(f"Spot Price (S): {S:,.2f}")
      print(f"Strike (K):     {K:,.2f}")
      print(f"Time (T):       {T:.3f} years")
      print(f"r:              {r:.4f}")
      print(f"sigma:          {sigma:.4f}")
      print(f"Dividend (q):   {q:.4f}")
      print(f"Call Price (C): {C:,.4f}")
      print(f"Put Price (P):  {P:,.4f}")
      print("------------------------------------\n")
      
      return C, P
