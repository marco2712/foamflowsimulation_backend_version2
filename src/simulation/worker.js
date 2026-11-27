import { parentPort, workerData } from 'worker_threads';
import { fw, D_cap, nD_LE } from './constitutive.js';
import { dynFactor } from '../capillary/coupling.js';
import { advectionUpwindX, divDiffusion, applyBCs } from './operators.js';

let paused = false;
let running = true;

// Capillary coupling (see capillary/coupling.js): sigma is replaced by the effective gamma*cos(theta)
const cap = workerData.cap || { sigmaScale: 1, dyn: null };
const params = { ...workerData.params, sigma: workerData.params.sigma * cap.sigmaScale };
const {
  Nx, Nz, Sw_plus, Sw_minus, Sw_star, A, u1, u2, d,
  Swc, Sgr, phi1, phi2, k1, k2, dt, Tmax, L, frameEvery,
} = params;

// Control de cada cuántos pasos se envía un frame (múltiplo de 10, mínimo 10)
const FRAME_EVERY = Math.max(10, Math.round((frameEvery ?? 100) / 10) * 10);

const dx = L / Nx;
const dz = (2 * d) / Nz;
const invDx2 = 1.0 / (dx * dx);
const invDz2 = 1.0 / (dz * dz);
const thetaRaw = Number(params.theta_s);
const thetaS = Number.isFinite(thetaRaw) ? thetaRaw : 1;
const numCells = (Nz + 1) * (Nx + 1);

const Sw      = new Float64Array(numCells);
const nD      = new Float64Array(numCells);
const phi_2d  = new Float64Array(numCells);
const u_2d    = new Float64Array(numCells);
const isL1Arr = new Uint8Array(numCells);
const M_phi   = new Float64Array(numCells);
const w_area  = new Float64Array(numCells);

// Initialize geometry (z runs from -d to +d)
for (let j = 0; j <= Nz; j++) {
  const wz    = (j === 0 || j === Nz) ? 0.5 : 1.0;
  const z_val = -d + j * dz;
  const isL1  = z_val > 0.0;
  for (let i = 0; i <= Nx; i++) {
    const wx  = (i === 0 || i === Nx) ? 0.5 : 1.0;
    const idx = j * (Nx + 1) + i;
    isL1Arr[idx] = isL1 ? 1 : 0;
    phi_2d[idx]  = isL1 ? phi1 : phi2;
    u_2d[idx]    = isL1 ? u1   : u2;
    w_area[idx]  = wz * wx * dx * dz;
    M_phi[idx]   = phi_2d[idx] * w_area[idx];
    Sw[idx]      = Sw_plus;
  }
}

// High-resolution 1D lookup tables (4096 points with linear interpolation) for local-equilibrium
// fw(Sw), D_cap(Sw), nD_LE(Sw), and |dfw/dSw| in Layer 1 and Layer 2.
const LUT_N = 4096;
const swMin = Swc + 1e-6;
const swMax = 1.0 - Sgr - 1e-6;
const swSpan = Math.max(1e-9, swMax - swMin);
const invSwStep = (LUT_N - 1) / swSpan;

const fwLut1  = new Float64Array(LUT_N);
const fwLut2  = new Float64Array(LUT_N);
const dcLut1  = new Float64Array(LUT_N);
const dcLut2  = new Float64Array(LUT_N);
const ndLut   = new Float64Array(LUT_N);
const dfwLut1 = new Float64Array(LUT_N);
const dfwLut2 = new Float64Array(LUT_N);

for (let m = 0; m < LUT_N; m++) {
  const s = swMin + (swSpan * m) / (LUT_N - 1);
  const ndVal = Math.max(0, Math.min(nD_LE(s, params), 1));
  ndLut[m] = ndVal;
  fwLut1[m] = fw(s, ndVal, k1, params);
  fwLut2[m] = fw(s, ndVal, k2, params);
  dcLut1[m] = Math.max(0, D_cap(s, ndVal, k1, phi1, params));
  dcLut2[m] = Math.max(0, D_cap(s, ndVal, k2, phi2, params));

  const h = 1e-6;
  const sp = Math.min(swMax, s + h);
  const sm = Math.max(swMin, s - h);
  const denom = Math.max(1e-12, sp - sm);
  const np = nD_LE(sp, params);
  const nm = nD_LE(sm, params);
  dfwLut1[m] = Math.abs(fw(sp, np, k1, params) - fw(sm, nm, k1, params)) / denom;
  dfwLut2[m] = Math.abs(fw(sp, np, k2, params) - fw(sm, nm, k2, params)) / denom;
}

function evalLut(lut, s) {
  const pos = (s - swMin) * invSwStep;
  if (pos <= 0) return lut[0];
  if (pos >= LUT_N - 1) return lut[LUT_N - 1];
  const i = pos | 0;
  const frac = pos - i;
  return lut[i] + frac * (lut[i + 1] - lut[i]);
}

// Initial nD = nD_LE(Sw_plus) everywhere
for (let idx = 0; idx < numCells; idx++) {
  nD[idx] = evalLut(ndLut, Sw[idx]);
}

const nD_inj = Sw_minus > Sw_star ? Math.tanh(A * (Sw_minus - Sw_star)) : 0.0;
applyBCs(Sw, nD, Nx, Nz, Sw_minus, nD_inj);

const UMBRAL = (Sw_minus + Sw_plus) / 2.0;
const Nt = Math.floor(Tmax / dt);
let step = 1;
let t = 0;

const layer1Rows = [], layer2Rows = [];
for (let j = 0; j <= Nz; j++) {
  (-d + j * dz) > 0 ? layer1Rows.push(j) : layer2Rows.push(j);
}

// Extract row of each layer for 1D profiles (0 = interface z=0, 1 = outer boundary; default 0.5 = mid-layer)
const rL1 = Math.max(0, Math.min(1, Number.isFinite(params.zExtractL1) ? params.zExtractL1 : 0.5));
const rL2 = Math.max(0, Math.min(1, Number.isFinite(params.zExtractL2) ? params.zExtractL2 : 0.5));
const midL1 = layer1Rows[Math.round(rL1 * Math.max(0, layer1Rows.length - 1))] ?? Math.floor((Nz * 3) / 4);
const midL2 = layer2Rows[Math.round((1 - rL2) * Math.max(0, layer2Rows.length - 1))] ?? Math.floor((Nz * 1) / 4);

function getFrontPos(arr, thresh, rows) {
  const pos = [];
  for (const j of rows) {
    let found = false;
    for (let i = 0; i <= Nx; i++) {
      const idx = j * (Nx + 1) + i;
      if (arr[idx] >= thresh) {
        if (i === 0) {
          pos.push(0);
        } else {
          const y0 = arr[j * (Nx + 1) + i - 1], y1 = arr[idx];
          const dy = y1 - y0;
          const xf = Math.abs(dy) > 1e-12 ? (i - 1) * dx + ((thresh - y0) * dx) / dy : i * dx;
          pos.push(Math.max(0, Math.min(xf, L)));
        }
        found = true;
        break;
      }
    }
    if (!found) pos.push(L);
  }
  return pos.length ? pos.reduce((a, b) => a + b, 0) / pos.length : 0;
}

const hist_t = [], hist_fp1 = [], hist_fp2 = [];
const hist_fp1_nD = [], hist_fp2_nD = [];
const hist_trans_z = [], hist_rec_pct = [];

const Fw = new Float64Array(numCells);
const Dc = new Float64Array(numCells);
const Sw_new = new Float64Array(numCells);
const diffGeomFactor = invDx2 + (thetaS > 0 ? invDz2 : 0.0);
const CFL_SAFE = 0.65;
const MAX_SUBSTEPS = 80;

function computeRequiredSubsteps(fCap) {
  let maxRate = 0.0;
  for (let idx = 0; idx < numCells; idx++) {
    const s = Sw[idx];
    const isL1 = isL1Arr[idx] === 1;
    const phi = phi_2d[idx];
    const dfw = evalLut(isL1 ? dfwLut1 : dfwLut2, s);
    const dc = fCap * evalLut(isL1 ? dcLut1 : dcLut2, s);
    const rate = (u_2d[idx] * dfw) / (phi * dx) + (2.0 * dc * diffGeomFactor) / phi;
    if (rate > maxRate) maxRate = rate;
  }
  const nSub = Math.ceil((dt * maxRate) / CFL_SAFE);
  return Math.max(1, Math.min(MAX_SUBSTEPS, nSub));
}

function advanceSubstep(dtSub, fCap) {
  for (let idx = 0; idx < numCells; idx++) {
    const s = Sw[idx];
    const isL1 = isL1Arr[idx] === 1;
    Fw[idx] = evalLut(isL1 ? fwLut1 : fwLut2, s);
    Dc[idx] = fCap * evalLut(isL1 ? dcLut1 : dcLut2, s);
  }

  const adv  = advectionUpwindX(Fw, u_2d, Nx, Nz, dx);
  const diff = divDiffusion(Sw, Dc, Nx, Nz, dx, dz, thetaS);

  for (let idx = 0; idx < numCells; idx++) {
    const rhs = (-adv[idx] + diff[idx]) * w_area[idx];
    const sNext = Math.max(swMin, Math.min(Sw[idx] + (dtSub * rhs) / (M_phi[idx] + 1e-30), swMax));
    Sw_new[idx] = sNext;
  }

  for (let idx = 0; idx < numCells; idx++) {
    const sNext = Sw_new[idx];
    Sw[idx] = sNext;
    nD[idx] = evalLut(ndLut, sNext);
  }
  applyBCs(Sw, nD, Nx, Nz, Sw_minus, nD_inj);
}

function runBatch() {
  if (!running) return;
  if (paused) { setTimeout(runBatch, 30); return; }

  const batchSize = Math.min(50, Nt - step + 1);

  for (let b = 0; b < batchSize; b++) {
    const tStart = (step - 1) * dt;
    t = step * dt;

    const fCapInit = dynFactor(cap.dyn, Math.max(tStart + 0.5 * dt, 1e-12));
    const nSub = computeRequiredSubsteps(fCapInit);
    const dtSub = dt / nSub;

    for (let s = 0; s < nSub; s++) {
      const tMid = tStart + (s + 0.5) * dtSub;
      const fCap = dynFactor(cap.dyn, Math.max(tMid, 1e-12));
      advanceSubstep(dtSub, fCap);
    }

    // Crossflow at z=0 interface
    const midJ = Math.floor(Nz / 2);
    let q_cross = 0;
    if (thetaS > 0) {
      for (let i = 0; i <= Nx; i++) {
        const idxLo = midJ * (Nx + 1) + i;
        const idxHi = (midJ + 1) * (Nx + 1) + i;
        const Di = (2 * Dc[idxHi] * Dc[idxLo]) / (Dc[idxHi] + Dc[idxLo] + 1e-30);
        q_cross += ((Di * (Sw[idxLo] - Sw[idxHi])) / dz) * dx;
      }
    }

    // Swept area
    let swept = 0;
    for (let idx = 0; idx < numCells; idx++) {
      if (Sw[idx] > Sw_minus + 0.01) swept += dx * dz;
    }
    const rec_pct = (100 * swept) / (L * 2 * d);

    if (step % 10 === 0) {
      hist_t.push(t);
      hist_fp1.push(getFrontPos(Sw, UMBRAL, layer1Rows));
      hist_fp2.push(getFrontPos(Sw, UMBRAL, layer2Rows));
      hist_fp1_nD.push(getFrontPos(nD, 0.5, layer1Rows));
      hist_fp2_nD.push(getFrontPos(nD, 0.5, layer2Rows));
      hist_trans_z.push(q_cross);
      hist_rec_pct.push(rec_pct);
    }

    step++;
  }

  if (step % FRAME_EVERY < batchSize || step > Nt) {
    const SwL1 = [], SwL2 = [], nDL1 = [], nDL2 = [];
    for (let i = 0; i <= Nx; i++) {
      SwL1.push(Sw[midL1 * (Nx + 1) + i]);
      SwL2.push(Sw[midL2 * (Nx + 1) + i]);
      nDL1.push(nD[midL1 * (Nx + 1) + i]);
      nDL2.push(nD[midL2 * (Nx + 1) + i]);
    }

    parentPort.postMessage({
      type: 'frame', step, t,
      Sw: Array.from(Sw),
      nD: Array.from(nD),
      profiles1D: { SwL1, SwL2, nDL1, nDL2 },
      metricsBatch: {
        hist_t:       [...hist_t],       hist_fp1:    [...hist_fp1],
        hist_fp2:     [...hist_fp2],     hist_fp1_nD: [...hist_fp1_nD],
        hist_fp2_nD:  [...hist_fp2_nD],  hist_trans_z:[...hist_trans_z],
        hist_rec_pct: [...hist_rec_pct],
      },
    });
    hist_t.length = hist_fp1.length = hist_fp2.length = 0;
    hist_fp1_nD.length = hist_fp2_nD.length = 0;
    hist_trans_z.length = hist_rec_pct.length = 0;
  }

  if (step > Nt) {
    running = false;
    parentPort.postMessage({ type: 'done', summary: {} });
  } else {
    setTimeout(runBatch, 0);
  }
}

parentPort.on('message', ({ command }) => {
  if (command === 'pause')  paused = true;
  if (command === 'resume') paused = false;
  if (command === 'stop')   running = false;
});

runBatch();
