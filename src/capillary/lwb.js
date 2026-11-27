/**
 * Lucas-Washburn-Bosanquet (LWB) capillary penetration model.
 *
 * ODE (cylindrical tube of radius r, SI units):
 *   d/dt( pi r^2 rho x x' ) + 8 pi eta x x' = 2 pi r gamma cos(theta) - pi r^2 rho g x sin(beta)
 *
 * Expanding and dividing by pi r^2 rho x gives the explicit acceleration:
 *   x'' = [ 2 gamma cos(theta)/r - rho g x sin(beta) - 8 eta x x'/r^2 - rho x'^2 ] / (rho x)
 *
 * Singularity at x = 0:  we start at a tiny x0 with the analytical inertial
 * velocity  v0 = sqrt(2 gamma cos(theta) / (rho r))  (x = v0 t as t -> 0),
 * i.e. t_start = x0 / v0.
 *
 * Solver: adaptive Dormand-Prince Runge-Kutta 5(4).
 */

export const G = 9.80665;

/** Validate and normalise capillary parameters (SI; angles in degrees in the file). */
export function normalizeCapillary(c = {}) {
  const errors = [];
  const num = (k, req = true) => {
    const v = c[k];
    if (v === undefined || v === null) {
      if (req) errors.push(`capillary.${k} es obligatorio`);
      return undefined;
    }
    if (typeof v !== 'number' || !Number.isFinite(v)) errors.push(`capillary.${k} debe ser un número`);
    return v;
  };
  const rho = num('rho');
  const eta = num('eta');
  const gamma = num('gamma');
  const theta_deg = num('theta_deg');
  const r = num('r');
  const beta_deg = c.beta_deg ?? 0;
  const tmax = c.tmax ?? 1.0;
  const n_points = Math.round(c.n_points ?? 200);
  for (const [k, v] of Object.entries({ rho, eta, gamma, r, tmax })) {
    if (v !== undefined && !(v > 0)) errors.push(`capillary.${k} debe ser > 0`);
  }
  if (theta_deg !== undefined && !(theta_deg >= 0 && theta_deg < 90)) {
    errors.push('capillary.theta_deg debe estar en [0, 90): para ángulos ≥ 90° no hay penetración espontánea');
  }
  if (!(n_points >= 10 && n_points <= 5000)) errors.push('capillary.n_points debe estar entre 10 y 5000');
  if (errors.length) return { errors };
  return {
    errors,
    value: {
      rho, eta, gamma, r, tmax, n_points,
      theta: (theta_deg * Math.PI) / 180,
      theta_deg,
      beta: (beta_deg * Math.PI) / 180,
      beta_deg,
      gravity: c.gravity === undefined ? Math.abs(beta_deg) > 0 : !!c.gravity,
    },
  };
}

/** Inertial initial velocity v0 = sqrt(2 gamma cos(theta) / (rho r)). */
export function inertialVelocity(p) {
  return Math.sqrt((2 * p.gamma * Math.cos(p.theta)) / (p.rho * p.r));
}

/** Classical Lucas-Washburn: x = sqrt(gamma r cos(theta) t / (2 eta)). */
export function lucasWashburn(p, t) {
  return Math.sqrt((p.gamma * p.r * Math.cos(p.theta) * t) / (2 * p.eta));
}

/** Viscous relaxation time tau = rho r^2 / (8 eta): inertial -> viscous transition. */
export function transitionTime(p) {
  return (p.rho * p.r * p.r) / (8 * p.eta);
}

function accel(x, v, p) {
  const g = p.gravity ? G : 0;
  const drive = (2 * p.gamma * Math.cos(p.theta)) / p.r - p.rho * g * x * Math.sin(p.beta);
  return (drive - (8 * p.eta * x * v) / (p.r * p.r) - p.rho * v * v) / (p.rho * x);
}

// Dormand-Prince coefficients
const C = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
const A = [
  [],
  [1 / 5],
  [3 / 40, 9 / 40],
  [44 / 45, -56 / 15, 32 / 9],
  [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
  [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
  [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
];
const B5 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0];
const B4 = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40];

function f(y, p) {
  return [y[1], accel(y[0], y[1], p)];
}

/** Integrate from (t0, y0) to t1 with adaptive DP45. Returns y(t1). */
function integrate(t0, y0, t1, p, h0, rtol, atol) {
  let t = t0;
  let y = y0.slice();
  let h = Math.min(h0, t1 - t0);
  let steps = 0;
  while (t < t1) {
    if (++steps > 2e6) throw new Error('LWB: demasiados pasos del integrador');
    if (t + h > t1) h = t1 - t;
    const k = [];
    for (let s = 0; s < 7; s++) {
      const ys = y.map((yi, i) => {
        let acc = 0;
        for (let j = 0; j < s; j++) acc += A[s][j] * k[j][i];
        return yi + h * acc;
      });
      k.push(f(ys, p));
    }
    const y5 = y.map((yi, i) => yi + h * B5.reduce((a, b, j) => a + b * k[j][i], 0));
    const y4 = y.map((yi, i) => yi + h * B4.reduce((a, b, j) => a + b * k[j][i], 0));
    let err = 0;
    for (let i = 0; i < 2; i++) {
      const sc = atol[i] + rtol * Math.max(Math.abs(y[i]), Math.abs(y5[i]));
      err = Math.max(err, Math.abs(y5[i] - y4[i]) / sc);
    }
    if (err <= 1 && y5.every(Number.isFinite) && y5[0] > 0) {
      t += h;
      y = y5;
    }
    const fac = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * Math.pow(err, -0.2)));
    h *= fac;
    if (h < 1e-18) throw new Error('LWB: paso del integrador demasiado pequeño');
  }
  return y;
}

/**
 * Solve LWB at the given increasing times.
 * Returns { t, x, v } arrays. Times before t_start (= x0/v0) use x = v0 t.
 */
export function solveLWB(p, times, { x0 = 1e-9, rtol = 1e-8 } = {}) {
  const v0 = inertialVelocity(p);
  const tStart = x0 / v0;
  let t = tStart;
  let y = [x0, v0];
  const xs = [];
  const vs = [];
  const atol = [1e-14, 1e-12];
  let h = Math.max(tStart * 0.1, 1e-12);
  for (const tt of times) {
    if (tt <= tStart) {
      xs.push(v0 * tt);
      vs.push(v0);
      continue;
    }
    if (tt > t) y = integrate(t, y, tt, p, h, rtol, atol);
    t = Math.max(t, tt);
    h = Math.max(1e-12, (tt - tStart) * 1e-3);
    xs.push(y[0]);
    vs.push(y[1]);
  }
  return { t: times.slice(), x: xs, v: vs };
}

export function logspace(a, b, n) {
  const la = Math.log10(a);
  const lb = Math.log10(b);
  return Array.from({ length: n }, (_, i) => Math.pow(10, la + ((lb - la) * i) / (n - 1)));
}

/** Full forward simulation with LW comparison and regime diagnostics. */
export function simulateCapillary(p) {
  const tmin = Math.min(1e-6, p.tmax * 1e-4);
  const t = logspace(tmin, p.tmax, p.n_points);
  const { x, v } = solveLWB(p, t);
  const x_lw = t.map((tt) => lucasWashburn(p, tt));
  // local log-log slope d ln x / d ln t : 1 = inertial, 0.5 = viscous
  const slope = x.map((_, i) => {
    const a = Math.max(0, i - 1);
    const b = Math.min(x.length - 1, i + 1);
    return (Math.log(x[b]) - Math.log(x[a])) / (Math.log(t[b]) - Math.log(t[a]));
  });
  const v0 = inertialVelocity(p);
  const tau = transitionTime(p);
  const h_jurin = p.gravity && Math.sin(p.beta) > 0
    ? (2 * p.gamma * Math.cos(p.theta)) / (p.rho * G * p.r * Math.sin(p.beta)) : null;
  return {
    t, x_lwb: x, v_lwb: v, x_lw, slope,
    summary: {
      v0_inertial: v0,
      tau_transition: tau,
      jurin_height: h_jurin,
      x_final_lwb: x[x.length - 1],
      x_final_lw: x_lw[x_lw.length - 1],
    },
  };
}

// ---------------------------------------------------------------------------
// Inverse capillary: identify unknown parameters from measured x(t)
// ---------------------------------------------------------------------------

const TRANSFORMS = {
  r: { to: (v) => Math.log(v), from: (z) => Math.exp(z), get: (p) => p.r, set: (p, v) => { p.r = v; } },
  eta: { to: (v) => Math.log(v), from: (z) => Math.exp(z), get: (p) => p.eta, set: (p, v) => { p.eta = v; } },
  gamma: { to: (v) => Math.log(v), from: (z) => Math.exp(z), get: (p) => p.gamma, set: (p, v) => { p.gamma = v; } },
  theta_deg: {
    to: (v) => Math.log(v / (89.9 - v)),
    from: (z) => 89.9 / (1 + Math.exp(-z)),
    get: (p) => p.theta_deg,
    set: (p, v) => { p.theta_deg = v; p.theta = (v * Math.PI) / 180; },
  },
};

function nelderMead(fn, x0, { iters = 300, step = 0.3, tol = 1e-12 } = {}) {
  const n = x0.length;
  let simplex = [x0.slice()];
  for (let i = 0; i < n; i++) {
    const p = x0.slice();
    p[i] += step;
    simplex.push(p);
  }
  let vals = simplex.map(fn);
  for (let it = 0; it < iters; it++) {
    const order = vals.map((v, i) => i).sort((a, b) => vals[a] - vals[b]);
    simplex = order.map((i) => simplex[i]);
    vals = order.map((i) => vals[i]);
    if (Math.abs(vals[n] - vals[0]) < tol) break;
    const cen = Array(n).fill(0);
    for (let i = 0; i < n; i++) simplex[i].forEach((v, j) => { cen[j] += v / n; });
    const pt = (c) => cen.map((v, j) => v + c * (simplex[n][j] - v));
    const xr = pt(-1);
    const fr = fn(xr);
    if (fr < vals[0]) {
      const xe = pt(-2);
      const fe = fn(xe);
      if (fe < fr) { simplex[n] = xe; vals[n] = fe; } else { simplex[n] = xr; vals[n] = fr; }
    } else if (fr < vals[n - 1]) {
      simplex[n] = xr; vals[n] = fr;
    } else {
      const xc = pt(fr < vals[n] ? -0.5 : 0.5);
      const fc = fn(xc);
      if (fc < Math.min(fr, vals[n])) { simplex[n] = xc; vals[n] = fc; } else {
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[i].map((v, j) => simplex[0][j] + 0.5 * (v - simplex[0][j]));
          vals[i] = fn(simplex[i]);
        }
      }
    }
  }
  const best = vals.indexOf(Math.min(...vals));
  return { x: simplex[best], value: vals[best] };
}

/**
 * Inverse problem: fit parameters listed in `fit` (subset of r, eta, gamma, theta_deg)
 * so that the LWB model reproduces measured penetration x_i(t_i).
 * Cost = sum (ln x_model - ln x_meas)^2  (relative error, valid across both regimes).
 */
export function inverseCapillary(base, data, fit = ['r']) {
  const errors = [];
  if (!data || !Array.isArray(data.t) || !Array.isArray(data.x) || data.t.length !== data.x.length) {
    errors.push('inverse.t e inverse.x deben ser arreglos de igual longitud');
  } else if (data.t.length < Math.max(3, fit.length + 1)) {
    errors.push('Se necesitan más puntos medidos que parámetros a identificar');
  } else if (data.t.some((v, i) => !(v > 0) || !(data.x[i] > 0) || (i && v <= data.t[i - 1]))) {
    errors.push('t debe ser estrictamente creciente y t, x > 0');
  }
  for (const k of fit) if (!TRANSFORMS[k]) errors.push(`Parámetro no ajustable: ${k} (use r, eta, gamma, theta_deg)`);
  if (errors.length) return { errors };

  const idx = data.t.map((_, i) => i);
  const tm = idx.map((i) => data.t[i]);
  const xm = idx.map((i) => data.x[i]);
  const build = (z) => {
    const p = { ...base };
    fit.forEach((k, i) => TRANSFORMS[k].set(p, TRANSFORMS[k].from(z[i])));
    return p;
  };
  const cost = (z) => {
    try {
      const p = build(z);
      if (p.r < 1e-9 || p.r > 1) return 1e12;
      const { x } = solveLWB(p, tm, { rtol: 1e-6 });
      let s = 0;
      for (let i = 0; i < x.length; i++) s += (Math.log(x[i]) - Math.log(xm[i])) ** 2;
      return Number.isFinite(s) ? s : 1e12;
    } catch {
      return 1e12;
    }
  };
  const z0 = fit.map((k) => TRANSFORMS[k].to(TRANSFORMS[k].get(base)));
  const res = nelderMead(cost, z0);
  const pBest = build(res.x);
  const fitted = {};
  fit.forEach((k) => { fitted[k] = TRANSFORMS[k].get(pBest); });
  const model = solveLWB(pBest, tm).x;
  return {
    errors,
    fitted,
    rmse_log: Math.sqrt(res.value / tm.length),
    t: tm, x_measured: xm, x_fit: model,
    note: fit.length > 1
      ? 'Con varios parámetros, en el régimen viscoso solo es identificable el grupo gamma·r·cos(theta)/eta; use datos en el régimen inercial o fije parámetros conocidos.'
      : undefined,
  };
}
