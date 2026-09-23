"""
A toy malaria outbreak as a Hawkes process.

Companion code for the blog post "Hawkes processes in epidemiology".

The model
---------
Conditional intensity (cases per day):

    lambda*(t) = mu(t) + sum_{t_i < t} eta * g(t - t_i)

Background (imported cases), a seasonal cycle with a one-year period:

    mu(t) = m * (1 + a * cos(2*pi*(t - peak) / 365))      with 0 <= a < 1

Kernel (local transmission), a Rayleigh curve shifted by a fixed delay,
similar in shape to Unwin et al. (2021):

    g(u) = ((u - delay) / sigma^2) * exp(-(u - delay)^2 / (2 sigma^2))  for u > delay
    g(u) = 0                                                            otherwise

g is a probability density, so eta is the branching ratio.

The numbers below are illustrative, not estimates for any real setting.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy import stats

rng = np.random.default_rng(2026)

YEAR = 365.0
T = 5 * YEAR          # length of the observation window (days)
DELAY = 15.0          # fixed delay before a case can cause another (days)

TRUE = dict(
    m=0.10,           # average imported cases per day (about 3 per month)
    a=0.8,            # strength of the seasonal cycle in importations
    peak=200.0,       # day of the year when importations peak
    eta=0.6,          # branching ratio: local cases caused per case
    sigma=20.0,       # spread of the transmission delay (days)
)

IMPORTED_COLOUR = "#E07B39"
LOCAL_COLOUR = "#3A6EA5"


# ---------------------------------------------------------------------------
# Model pieces
# ---------------------------------------------------------------------------

def mu(t, m, a, peak):
    """Background rate of imported cases at time t."""
    return m * (1 + a * np.cos(2 * np.pi * (t - peak) / YEAR))


def mu_integral(t0, t1, m, a, peak):
    """Integral of mu from t0 to t1 (closed form)."""
    k = 2 * np.pi / YEAR
    return m * (t1 - t0) + m * a / k * (np.sin(k * (t1 - peak)) - np.sin(k * (t0 - peak)))


def g(u, sigma):
    """Delay density: zero before DELAY, then a Rayleigh curve."""
    v = np.asarray(u, dtype=float) - DELAY
    out = np.zeros_like(v)
    pos = v > 0
    out[pos] = v[pos] / sigma**2 * np.exp(-v[pos] ** 2 / (2 * sigma**2))
    return out


def g_cdf(u, sigma):
    """Integral of g from 0 to u."""
    v = np.asarray(u, dtype=float) - DELAY
    out = np.zeros_like(v)
    pos = v > 0
    out[pos] = 1 - np.exp(-v[pos] ** 2 / (2 * sigma**2))
    return out


def intensity(t_grid, times, m, a, peak, eta, sigma):
    """lambda*(t) on a grid of times, given the case times."""
    lam = mu(t_grid, m, a, peak)
    for ti in times:
        lam = lam + eta * g(t_grid - ti, sigma)
    return lam


# ---------------------------------------------------------------------------
# Simulation (the family-tree method)
# ---------------------------------------------------------------------------

def simulate(m, a, peak, eta, sigma, T):
    """
    Simulate by building the family tree directly.

    1. Draw imported cases from the seasonal background, by thinning a
       constant-rate Poisson process (the background is at most m * (1 + a)).
    2. Give every case a Poisson(eta) number of children, each after a delay
       drawn from g. Repeat generation by generation.

    Returns case times, whether each case was imported, and each case's parent
    (-1 for imported cases), all sorted by time.
    """
    # 1. Imported cases
    mu_max = m * (1 + a)
    n_candidates = rng.poisson(mu_max * T)
    candidates = rng.uniform(0, T, n_candidates)
    keep = rng.uniform(0, 1, n_candidates) < mu(candidates, m, a, peak) / mu_max
    times = list(np.sort(candidates[keep]))
    parent = [-1] * len(times)

    # 2. Local cases, one generation at a time
    current = list(range(len(times)))
    while current:
        next_generation = []
        for idx in current:
            n_children = rng.poisson(eta)
            # Rayleigh delay: DELAY + sigma * sqrt(-2 log U)
            delays = DELAY + sigma * np.sqrt(-2 * np.log(rng.uniform(size=n_children)))
            for d in delays:
                t_child = times[idx] + d
                if t_child < T:            # children after T are never observed
                    times.append(t_child)
                    parent.append(idx)
                    next_generation.append(len(times) - 1)
        current = next_generation

    times = np.array(times)
    parent = np.array(parent)
    order = np.argsort(times)
    new_index = np.empty_like(order)
    new_index[order] = np.arange(len(order))
    parent_sorted = np.where(parent[order] >= 0, new_index[parent[order]], -1)
    return times[order], parent_sorted == -1, parent_sorted


# ---------------------------------------------------------------------------
# Likelihood and fitting
# ---------------------------------------------------------------------------

def unpack(theta):
    """Map unconstrained parameters to the model's parameters."""
    log_m, logit_a, peak, log_eta, log_sigma = theta
    return (np.exp(log_m), 1 / (1 + np.exp(-logit_a)), peak % YEAR,
            np.exp(log_eta), np.exp(log_sigma))


def triggering_terms(times, sigma):
    """Matrix K[i, j] = g(t_i - t_j) for j < i (zero otherwise)."""
    diffs = times[:, None] - times[None, :]
    K = g(np.where(diffs > 0, diffs, 0.0), sigma)
    K[diffs <= 0] = 0.0
    return K


def neg_log_likelihood(theta, times, T):
    m, a, peak, eta, sigma = unpack(theta)
    K = triggering_terms(times, sigma)
    lam = mu(times, m, a, peak) + eta * K.sum(axis=1)
    if np.any(lam <= 0):
        return np.inf
    compensator = mu_integral(0, T, m, a, peak) + eta * g_cdf(T - times, sigma).sum()
    return -(np.log(lam).sum() - compensator)


def fit(times, T, n_starts=10):
    """Maximum likelihood, repeated from several random starting points."""
    best = None
    for _ in range(n_starts):
        start = np.array([
            np.log(rng.uniform(0.02, 0.3)),
            rng.normal(0, 1),
            rng.uniform(0, YEAR),
            np.log(rng.uniform(0.1, 0.9)),
            np.log(rng.uniform(5, 50)),
        ])
        result = minimize(neg_log_likelihood, start, args=(times, T), method="Nelder-Mead",
                          options=dict(maxiter=5000, xatol=1e-6, fatol=1e-6))
        if best is None or result.fun < best.fun:
            best = result
    return unpack(best.x), best


# ---------------------------------------------------------------------------
# Classification and goodness of fit
# ---------------------------------------------------------------------------

def prob_imported(times, m, a, peak, eta, sigma):
    """P(case i is imported) = mu(t_i) / lambda*(t_i)."""
    K = triggering_terms(times, sigma)
    background = mu(times, m, a, peak)
    return background / (background + eta * K.sum(axis=1))


def rescaled_gaps(times, m, a, peak, eta, sigma):
    """Time-rescaling: gaps in Lambda(t_i) should be Exp(1) if the model is right."""
    Lam = np.array([
        mu_integral(0, t, m, a, peak) + eta * g_cdf(t - times[times < t], sigma).sum()
        for t in times
    ])
    return np.diff(np.concatenate([[0.0], Lam]))


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def style(ax):
    ax.spines[["top", "right"]].set_visible(False)


def figure_poisson_vs_hawkes(path):
    """Figure 1: Poisson vs Hawkes with the same background rate (exponential kernel)."""
    T_demo, base, alpha, beta = 200.0, 0.1, 0.8 * 0.25, 0.25   # branching ratio 0.8

    # Poisson
    n = rng.poisson(base * T_demo)
    poisson_times = np.sort(rng.uniform(0, T_demo, n))

    # Hawkes with an exponential kernel, by Ogata's thinning: the kernel only
    # decays, so the current intensity is an upper bound until the next event.
    hawkes_times, t = [], 0.0
    while True:
        past = np.array(hawkes_times)
        bound = base + (alpha * np.exp(-beta * (t - past)).sum() if len(past) else 0)
        t += rng.exponential(1 / bound)
        if t > T_demo:
            break
        lam_t = base + (alpha * np.exp(-beta * (t - past)).sum() if len(past) else 0)
        if rng.uniform() < lam_t / bound:
            hawkes_times.append(t)
    hawkes_times = np.array(hawkes_times)

    grid = np.linspace(0, T_demo, 2000)
    lam_hawkes = base + np.array([
        alpha * np.exp(-beta * (s - hawkes_times[hawkes_times < s])).sum() for s in grid
    ])

    fig, axes = plt.subplots(2, 1, figsize=(9, 4.5), sharex=True,
                             gridspec_kw=dict(height_ratios=[3, 1]))
    ax = axes[0]
    ax.plot(grid, lam_hawkes, color=LOCAL_COLOUR, lw=1.2, label="Hawkes intensity")
    ax.axhline(base, color="grey", ls="--", lw=1, label="Poisson intensity (constant)")
    ax.set_ylabel("events per day")
    ax.legend(frameon=False, loc="upper left")
    style(ax)

    ax = axes[1]
    ax.eventplot([hawkes_times, poisson_times], lineoffsets=[1, 0], linelengths=0.7,
                 colors=[LOCAL_COLOUR, "grey"])
    ax.set_yticks([1, 0], ["Hawkes", "Poisson"])
    ax.set_xlabel("day")
    style(ax)
    fig.suptitle("Same background rate: Poisson events are scattered, Hawkes events come in bursts",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def figure_outbreak(times, imported, path):
    """Figure 2: the simulated outbreak, coloured by true origin."""
    grid = np.linspace(0, T, 4000)
    lam = intensity(grid, times, **TRUE)
    background = mu(grid, TRUE["m"], TRUE["a"], TRUE["peak"])

    fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex=True,
                             gridspec_kw=dict(height_ratios=[3, 1]))
    ax = axes[0]
    ax.fill_between(grid, 0, background, color=IMPORTED_COLOUR, alpha=0.35,
                    label=r"background $\mu(t)$ (imported)")
    ax.fill_between(grid, background, lam, color=LOCAL_COLOUR, alpha=0.35,
                    label="triggered part (local transmission)")
    ax.plot(grid, lam, color="black", lw=0.8, label=r"total intensity $\lambda^{\ast}(t)$")
    ax.set_ylabel("cases per day")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    style(ax)

    ax = axes[1]
    ax.eventplot([times[imported], times[~imported]], lineoffsets=[1, 0], linelengths=0.7,
                 colors=[IMPORTED_COLOUR, LOCAL_COLOUR])
    ax.set_yticks([1, 0], ["imported", "local"])
    ax.set_xlabel("day")
    style(ax)
    fig.suptitle(f"Simulated outbreak: {imported.sum()} imported and {(~imported).sum()} local cases",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def figure_classification(p_imp, imported, path):
    """Figure 3: estimated probability of being imported, split by true label."""
    fig, ax = plt.subplots(figsize=(7, 4))
    jitter = rng.uniform(-0.15, 0.15, len(p_imp))
    ax.scatter(np.where(imported, 1, 0) + jitter, p_imp, s=10, alpha=0.6,
               c=np.where(imported, IMPORTED_COLOUR, LOCAL_COLOUR))
    ax.axhline(0.5, color="grey", ls="--", lw=1)
    ax.set_xticks([0, 1], ["truly local", "truly imported"])
    ax.set_ylabel("estimated P(imported)")
    ax.set_ylim(-0.02, 1.02)
    style(ax)
    fig.suptitle("Classifying cases from timing alone", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def figure_qq(gaps, path):
    """Figure 4: Q-Q plot of rescaled gaps against Exp(1)."""
    gaps = np.sort(gaps)
    n = len(gaps)
    theoretical = stats.expon.ppf((np.arange(1, n + 1) - 0.5) / n)
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.scatter(theoretical, gaps, s=8, color=LOCAL_COLOUR, alpha=0.7)
    lim = max(theoretical.max(), gaps.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="grey", ls="--", lw=1)
    ax.set_xlabel("Exp(1) quantiles")
    ax.set_ylabel("rescaled gaps")
    ax.set_title("Time-rescaling check", fontsize=11)
    style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Run everything
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import os
    out = "figures"
    os.makedirs(out, exist_ok=True)

    figure_poisson_vs_hawkes(f"{out}/poisson_vs_hawkes.png")

    times, imported, parent = simulate(T=T, **TRUE)
    print(f"Simulated {len(times)} cases: {imported.sum()} imported, "
          f"{(~imported).sum()} local ({(~imported).mean():.0%} local)")
    figure_outbreak(times, imported, f"{out}/simulated_outbreak.png")

    (m, a, peak, eta, sigma), result = fit(times, T)
    print("\nParameter   true     fitted")
    for name, fitted in zip(["m", "a", "peak", "eta", "sigma"], [m, a, peak, eta, sigma]):
        print(f"{name:<10}{TRUE[name]:>7.3f}  {fitted:>9.3f}")

    p_imp = prob_imported(times, m, a, peak, eta, sigma)
    predicted_imported = p_imp > 0.5
    accuracy = (predicted_imported == imported).mean()
    auc = stats.mannwhitneyu(p_imp[imported], p_imp[~imported]).statistic / (
        imported.sum() * (~imported).sum())
    print(f"\nExpected number imported (sum of probabilities): {p_imp.sum():.0f} "
          f"(true: {imported.sum()})")
    print(f"Accuracy at a 0.5 threshold: {accuracy:.0%}")
    print(f"AUC: {auc:.2f}")
    figure_classification(p_imp, imported, f"{out}/classification.png")

    gaps = rescaled_gaps(times, m, a, peak, eta, sigma)
    ks = stats.kstest(gaps, "expon")
    print(f"\nTime-rescaling KS test: statistic = {ks.statistic:.3f}, p = {ks.pvalue:.2f}")
    figure_qq(gaps, f"{out}/time_rescaling_qq.png")
