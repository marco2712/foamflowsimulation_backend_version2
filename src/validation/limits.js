/**
 * Operating limits from "parametros_de_Foamflow2d.pdf" (Tablas 1 y 2).
 * Each check returns { key, level: 'ok'|'caution'|'fail', message }.
 *   ok      -> "Rango funciona bien"
 *   caution -> "Rango precautorio" (or outside recommended range)
 *   fail    -> "Límite de fallo": blocks starting the simulation
 */

const inR = (v, a, b) => v >= a && v <= b;

export const REQUIRED_KEYS = [
  'phi1', 'phi2', 'Swc', 'Sgr', 'mu_w', 'mu_g', 'k1', 'k2', 'Sw_star', 'A', 'Kc', 'theta_s',
  'u1', 'u2', 'd', 'sigma', 'c_cap', 'Sw_minus', 'Sw_plus', 'L', 'Nx', 'Nz', 'dt', 'Tmax',
];

export function validateFoamParams(p) {
  const out = [];
  const add = (key, level, message) => out.push({ key, level, message });

  // Basic sanity
  for (const k of REQUIRED_KEYS) {
    if (typeof p[k] !== 'number' || !Number.isFinite(p[k])) add(k, 'fail', `${k} falta o no es numérico`);
  }
  if (out.length) return out;
  for (const k of ['phi1', 'phi2', 'mu_w', 'mu_g', 'k1', 'k2', 'u1', 'u2', 'd', 'L', 'Nx', 'Nz', 'dt', 'Tmax']) {
    if (!(p[k] > 0)) add(k, 'fail', `${k} debe ser > 0`);
  }
  if (out.length) return out;

  // CFL advectivo (velocidad intersticial u/phi)
  const dx = p.L / p.Nx;
  const cfl = Math.max(p.u1 / p.phi1, p.u2 / p.phi2) * p.dt / dx;

  // Saturación de inyección Sw-
  if (p.Sw_minus < p.Swc) {
    add('Sw_minus', 'fail', 'Sw- < Swc: físicamente imposible (por debajo de la saturación connata)');
  } else if (p.Sw_minus <= p.Sw_star) {
    add('Sw_minus', 'caution', 'Sw- ≤ Sw*: el motor solo genera espuma en la inyección si Sw- > Sw* (nD_inj = tanh(A(Sw- − Sw*)))');
  } else if (inR(p.Sw_minus, 0.25, 0.35)) {
    add('Sw_minus', 'ok', 'Sw- en rango (0.25–0.35)');
  } else {
    add('Sw_minus', 'caution', `Sw- = ${p.Sw_minus} fuera de 0.25–0.35 (activo sobre Sw* para inyección de espuma)`);
  }

  // Saturación inicial Sw+
  if (p.Sw_plus <= p.Sw_star) add('Sw_plus', 'fail', 'Sw+ ≤ Sw*: el reservorio inicial carece de banco de agua');
  else if (inR(p.Sw_plus, 0.65, 0.80)) add('Sw_plus', 'ok', 'Sw+ en rango (0.65–0.80)');
  else add('Sw_plus', 'caution', 'Sw+ fuera de 0.65–0.80 (precautorio 0.50–0.65)');

  // Saturación crítica Sw*
  if (!(p.Sw_star < p.Sw_plus)) add('Sw_star', 'fail', 'Sw* ≥ Sw+: la espuma no se forma en el medio inicial');
  else if (inR(p.Sw_star, 0.35, 0.42)) add('Sw_star', 'ok', 'Sw* en rango (0.35–0.42)');
  else add('Sw_star', 'caution', 'Sw* fuera de 0.35–0.42');

  // Pendiente de transición A
  if (p.A > 1000 || p.A < 20) add('A', 'fail', 'A > 1000 (difusión artificial) o A < 20 (frente excesivamente difuso)');
  else if (inR(p.A, 100, 500)) add('A', 'ok', 'A en rango (100–500)');
  else add('A', 'caution', 'A fuera de 100–500 (precautorio hasta 800)');

  // Contraste de permeabilidad
  const kr = Math.max(p.k1 / p.k2, p.k2 / p.k1);
  if (kr > 50) add('k1/k2', 'fail', `Contraste K1/K2 = ${kr.toFixed(1)} > 50: gradiente vertical excesivo`);
  else if (kr <= 10) add('k1/k2', 'ok', `Contraste K1/K2 = ${kr.toFixed(2)} (1–10)`);
  else add('k1/k2', 'caution', `Contraste K1/K2 = ${kr.toFixed(1)} (precautorio 10–25)`);

  // Porosidad
  for (const k of ['phi1', 'phi2']) {
    if (p[k] <= 0.05) add(k, 'fail', `${k} ≤ 0.05: acelera artificialmente la advección nodal`);
    else if (inR(p[k], 0.15, 0.35)) add(k, 'ok', `${k} en rango (0.15–0.35)`);
    else add(k, 'caution', `${k} fuera de 0.15–0.35`);
  }

  // Velocidad de inyección
  for (const k of ['u1', 'u2']) {
    if (p[k] > 1e-3 && cfl > 1.0) add(k, 'fail', `${k} > 1e-3 m/s sin reducir dt: viola la condición CFL`);
    else if (inR(p[k], 1e-5, 1e-4)) add(k, 'ok', `${k} en rango (1e-5–1e-4 m/s)`);
    else add(k, 'caution', `${k} fuera de 1e-5–1e-4 m/s`);
  }

  // Viscosidades
  const ratio = p.mu_w / p.mu_g;
  if (inR(ratio, 20, 100)) add('mu_w/mu_g', 'ok', `μw/μg = ${ratio.toFixed(1)} (20–100)`);
  else add('mu_w/mu_g', 'caution', `μw/μg = ${ratio.toFixed(1)} fuera de 20–100`);

  // Malla
  const cells = p.Nx * p.Nz;
  if (cells > 50000) add('mesh', 'fail', `Malla ${p.Nx}×${p.Nz} = ${cells} celdas > 50 000: la UI se congela`);
  else if (cells <= 8000) add('mesh', 'ok', `Malla ${p.Nx}×${p.Nz} = ${cells} celdas (≤ 8 000)`);
  else add('mesh', 'caution', `Malla ${p.Nx}×${p.Nz} = ${cells} celdas (puede haber lag)`);

  // CFL
  if (cfl > 1.0) add('CFL', 'fail', `CFL = ${cfl.toFixed(3)} > 1: divergencia numérica (NaN/Inf). Reduzca u o dt`);
  else if (cfl > 0.85) add('CFL', 'caution', `CFL = ${cfl.toFixed(3)} (precautorio 0.85–0.95)`);
  else add('CFL', 'ok', `CFL = ${cfl.toFixed(3)}`);

  return out;
}

export function summarize(report) {
  return {
    errors: report.filter((r) => r.level === 'fail'),
    warnings: report.filter((r) => r.level === 'caution'),
    ok: report.filter((r) => r.level === 'ok'),
  };
}
