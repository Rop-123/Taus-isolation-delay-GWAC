"""SEIQR model with onset-to-isolation delay (Team TAUS).

States: S, E, I (infectious, not yet isolated), Q (isolated infectious), R.
  lambda = beta*I/N + beta*(1-eps)*Q/N
  dS=-lambda S; dE=lambda S - sigma E; dI=sigma E-(tau+gamma) I
  dQ=tau I - gamma_iso Q; dR=gamma I + gamma_iso Q
tau = 1/delay.  Cumulative symptomatic cases C: dC = sigma E.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

N_POP = 1e7


def jacobian(beta, sigma, gamma, tau, eps, gamma_iso):
    return np.array([[-sigma, beta, beta * (1 - eps)],
                     [sigma, -(tau + gamma), 0.0],
                     [0.0, tau, -gamma_iso]])


def growth_rate(beta, sigma, gamma, tau, eps, gamma_iso):
    return np.max(np.linalg.eigvals(jacobian(beta, sigma, gamma, tau, eps, gamma_iso)).real)


def re_value(beta, sigma, gamma, tau, eps, gamma_iso):
    """Effective reproduction number at the disease-free state (S ~ N)."""
    return beta * (1.0 / (tau + gamma) + tau / (tau + gamma) * (1 - eps) / gamma_iso)


def calibrate_beta(r, sigma, gamma, tau0, eps, gamma_iso):
    """beta such that the baseline (tau0) system grows at rate r."""
    f = lambda b: growth_rate(b, sigma, gamma, tau0, eps, gamma_iso) - r
    return brentq(f, 1e-6, 20.0)


def initial_state(beta, sigma, gamma, tau, eps, gamma_iso, inc0, n=N_POP):
    """Dominant eigenvector of the linearised (E,I,Q) system, scaled so sigma*E = inc0."""
    w, v = np.linalg.eig(jacobian(beta, sigma, gamma, tau, eps, gamma_iso))
    vec = np.abs(v[:, np.argmax(w.real)].real)
    vec = vec * (inc0 / (sigma * vec[0]))
    E, I, Q = vec
    return np.array([n - E - I - Q, E, I, Q, 0.0, 0.0])


def simulate(beta, sigma, gamma, tau, eps, gamma_iso, y0, t_end, n=N_POP):
    def rhs(t, y):
        S, E, I, Q, R, C = y
        lam = beta * I / n + beta * (1 - eps) * Q / n
        return [-lam * S, lam * S - sigma * E, sigma * E - (tau + gamma) * I,
                tau * I - gamma_iso * Q, gamma * I + gamma_iso * Q, sigma * E]
    t = np.arange(0, t_end + 1)
    sol = solve_ivp(rhs, (0, t_end), y0, t_eval=t, method="LSODA", rtol=1e-8, atol=1e-6)
    return sol.t, sol.y


def outcomes(t, y, horizon=365):
    """cumulative symptomatic cases at horizon; duration = first day E+I+Q < 1 (nan if never)."""
    C = y[5]
    cum = C[horizon] - C[0]
    active = y[1] + y[2] + y[3]
    below = np.where(active < 1.0)[0]
    dur = float(t[below[0]]) if len(below) else np.nan
    inc = np.diff(C)
    return cum, dur, inc


def discrete_matrix(beta, sigma, gamma, tau, eps, gamma_iso, dt):
    """Expected one-step map on (E,I,Q) of the binomial-chain model below."""
    es, et, eg, egi = np.exp(-sigma * dt), np.exp(-tau * dt), np.exp(-gamma * dt), np.exp(-gamma_iso * dt)
    return np.array([[es, beta * dt, beta * dt * (1 - eps)],
                     [1 - es, et * eg, 0.0],
                     [0.0, 1 - et, egi]])


def calibrate_beta_discrete(r, sigma, gamma, tau0, eps, gamma_iso, dt):
    f = lambda b: np.max(np.linalg.eigvals(discrete_matrix(b, sigma, gamma, tau0, eps, gamma_iso, dt)).real) - np.exp(r * dt)
    return brentq(f, 1e-6, 20.0)


def initial_state_discrete(beta, sigma, gamma, tau, eps, gamma_iso, inc0, dt, n=N_POP):
    w, v = np.linalg.eig(discrete_matrix(beta, sigma, gamma, tau, eps, gamma_iso, dt))
    vec = np.abs(v[:, np.argmax(w.real)].real)
    vec = vec * (inc0 * dt / ((1 - np.exp(-sigma * dt)) * vec[0]))
    E, I, Q = vec
    return np.array([n - E - I - Q, E, I, Q, 0.0, 0.0])


def stochastic(beta, sigma, gamma, tau, eps, gamma_iso, y0, t_end, rng, dt=0.25, n=N_POP):
    """Binomial-chain stochastic version; vectorised over replicates (params may be arrays)."""
    reps = np.size(beta)
    S, E, I, Q = [np.full(reps, np.round(v), dtype=np.int64) for v in y0[:4]]
    C = np.zeros(reps)
    ext = np.full(reps, np.nan)
    steps = int(t_end / dt)
    for k in range(steps):
        lam = beta * I / n + beta * (1 - eps) * Q / n
        newE = rng.binomial(S, 1 - np.exp(-lam * dt))
        newI = rng.binomial(E, 1 - np.exp(-sigma * dt))
        newQ = rng.binomial(I, 1 - np.exp(-tau * dt))
        remI = rng.binomial(I - newQ, 1 - np.exp(-gamma * dt))
        remQ = rng.binomial(Q, 1 - np.exp(-gamma_iso * dt))
        S = S - newE
        E = E + newE - newI
        I = I + newI - newQ - remI
        Q = Q + newQ - remQ
        if k * dt < 365:
            C += newI
        done = np.isnan(ext) & ((E + I + Q) == 0)
        ext[done] = (k + 1) * dt
    return C, ext
