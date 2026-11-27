/**
 * Coupling between the capillary models and the 2D foam simulation.
 *
 * In the foam model the capillary pressure is a Leverett-type law
 *     Pc(Sw) = sigma * sqrt(phi/k) * 0.022 * (1-Sw-Sgr)^c_cap / (Sw-Swc)
 * and it enters the saturation equation through the capillary diffusion
 *     D_cap = -lambda_g * fw * dPc/dSw ,
 * which controls the width of the foam/shock front.
 *
 * Modes
 *  - 'none'    : unchanged foam model (sigma from the "foam" section).
 *  - 'lwb'     : Lucas-Washburn-Bosanquet.
 *        (a) static: sigma -> gamma*cos(theta) of the capillary section (wetting-aware Leverett).
 *        (b) dynamic: D_cap is multiplied by f(t) = min(1, v_LWB(t)/v_LW(t)).
 *            The inertial (Bosanquet) regime delays the capillary action at early times:
 *            v_LW -> infinity as t -> 0 while v_LWB stays finite (f -> 0), and f -> 1
 *            once t >> tau = rho r^2 / (8 eta).
 *  - 'inverse' : parameters identified by the inverse fit (gamma, theta, ...) replace the
 *        file values; static coupling sigma -> gamma_fit*cos(theta_fit). No dynamic factor
 *        (the fit calibrates the wetting parameters from measured data).
 */
import { solveLWB, lucasWashburn, logspace, inertialVelocity } from './lwb.js';

export const CAP_MODES = ['none', 'lwb', 'inverse'];

export function buildCoupling(mode, capNorm, fitted, foamParams, Tmax) {
  if (mode === 'none') {
    return { mode, sigmaScale: 1, sigmaEff: foamParams.sigma, dyn: null, info: 'Sin acoplamiento capilar (modelo base).' };
  }
  const p = { ...capNorm };
  if (mode === 'inverse') {
    for (const [k, v] of Object.entries(fitted || {})) {
      p[k] = v;
      if (k === 'theta_deg') p.theta = (v * Math.PI) / 180;
    }
  }
  const sigmaEff = p.gamma * Math.cos(p.theta);
  const out = {
    mode,
    sigmaEff,
    sigmaScale: sigmaEff / foamParams.sigma,
    dyn: null,
    info: `σ_ef = γ·cosθ = ${sigmaEff.toExponential(3)} N/m (σ del archivo foam: ${foamParams.sigma})`,
  };
  if (mode === 'lwb') {
    const tau = (p.rho * p.r * p.r) / (8 * p.eta);
    const tmin = Math.min(1e-9, tau * 1e-3);
    const t = logspace(tmin, Math.max(Tmax, tau * 100), 120);
    const { x } = solveLWB(p, t);
    // velocity ratio from the numerical x(t) and analytical LW derivative v_LW = x_LW/(2t)
    const v = t.map((tt, i) => {
      const a = Math.max(0, i - 1);
      const b = Math.min(t.length - 1, i + 1);
      return (x[b] - x[a]) / (t[b] - t[a]);
    });
    const f = t.map((tt, i) => Math.min(1, Math.max(0, v[i] / (lucasWashburn(p, tt) / (2 * tt)))));
    out.dyn = { t, f };
    out.info += `; factor dinámico LWB con τ = ${tau.toExponential(2)} s, v0 = ${inertialVelocity(p).toExponential(2)} m/s`;
  }
  return out;
}

/** Linear interpolation of the dynamic factor in log(t). */
export function dynFactor(dyn, t) {
  if (!dyn) return 1;
  const { t: ts, f } = dyn;
  if (t <= ts[0]) return f[0];
  if (t >= ts[ts.length - 1]) return f[f.length - 1];
  let lo = 0, hi = ts.length - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (ts[m] <= t) lo = m; else hi = m;
  }
  const w = (Math.log(t) - Math.log(ts[lo])) / (Math.log(ts[hi]) - Math.log(ts[lo]));
  return f[lo] + w * (f[hi] - f[lo]);
}
