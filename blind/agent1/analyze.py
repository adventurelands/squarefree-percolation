"""Analysis helpers: effective exponent alpha_eff(t) and run summaries."""
import json
import sys
import numpy as np


def load(fname):
    with open(fname) as f:
        if fname.endswith(".jsonl"):
            return [json.loads(line) for line in f if line.strip()]
        return json.load(f)


def avg_runs(results, key):
    """Average rho(t) over runs matching (construction, p). key=(cons,p)."""
    sel = [r for r in results if (r["construction"], round(r["p"], 5)) == key]
    if not sel:
        return None, None, 0
    t = np.array(sel[0]["t"], dtype=float)
    rhos = np.array([r["rho"] for r in sel])
    return t, rhos, len(sel)


def alpha_eff(t, rho, half_window=0.25):
    """Local slope -d ln rho / d ln t via linear fit in sliding log-window.

    half_window: half-width of the fit window in decades.
    Returns (t_mid, alpha) arrays. Points with rho<=0 dropped.
    """
    good = rho > 0
    lt, lr = np.log(t[good]), np.log(rho[good])
    out_t, out_a = [], []
    hw = half_window * np.log(10)
    for i in range(len(lt)):
        m = (lt >= lt[i] - hw) & (lt <= lt[i] + hw)
        if m.sum() < 6:
            continue
        A = np.polyfit(lt[m], lr[m], 1)
        out_t.append(np.exp(lt[i]))
        out_a.append(-A[0])
    return np.array(out_t), np.array(out_a)


def summary_table(fname):
    results = load(fname)
    keys = sorted({(r["construction"], round(r["p"], 5)) for r in results})
    print(f"{'cons':4} {'p':8} {'n':2} {'rho(1e3)':>10} {'rho(1e4)':>10} "
          f"{'rho(T)':>10} {'a_eff@1e3':>10} {'a_eff@T/3':>10} {'a_eff@T':>9}")
    for key in keys:
        t, rhos, n = avg_runs(results, key)
        rho = rhos.mean(axis=0)
        T = t[-1]

        def rho_at(tt):
            i = np.argmin(np.abs(t - tt))
            return rho[i]
        ta, aa = alpha_eff(t, rho)
        def a_at(tt):
            if len(ta) == 0:
                return np.nan
            i = np.argmin(np.abs(np.log(ta) - np.log(tt)))
            return aa[i]
        print(f"{key[0]:4} {key[1]:<8} {n:2} {rho_at(1e3):10.3e} "
              f"{rho_at(1e4):10.3e} {rho[-1]:10.3e} {a_at(1e3):10.3f} "
              f"{a_at(T/3):10.3f} {a_at(T):9.3f}")


if __name__ == "__main__":
    for f in sys.argv[1:]:
        print(f"== {f}")
        summary_table(f)
