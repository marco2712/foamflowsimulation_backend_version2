import express from 'express';
import path from 'path';
import { fileURLToPath } from 'url';
import WebSocket from 'ws';
import { state, resetState } from '../store/simState.js';
import { validateFoamParams, summarize, REQUIRED_KEYS } from '../validation/limits.js';
import { normalizeCapillary } from '../capillary/lwb.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const router = express.Router();

const CAPILLARY_KEYS = ['rho', 'eta', 'gamma', 'theta_deg', 'r', 'beta_deg', 'tmax', 'n_points', 'gravity'];

function coerceNumbers(obj) {
  if (!obj || typeof obj !== 'object') return obj;
  const out = {};
  for (const [k, v] of Object.entries(obj)) {
    if (typeof v === 'string' && v.trim() !== '' && Number.isFinite(Number(v))) {
      out[k] = Number(v);
    } else {
      out[k] = v;
    }
  }
  return out;
}

function broadcast(serverWss, data) {
  if (!serverWss) return;
  const json = JSON.stringify(data);
  serverWss.clients.forEach((client) => {
    if (client.readyState === WebSocket.OPEN) client.send(json);
  });
}

/**
 * Parameter file format (JSON):
 * {
 *   "foam":         { phi1, phi2, ... }
 *   "capillary":    { rho, eta, gamma, theta_deg, r, [beta_deg], [tmax], [n_points] }
 *   "measurements": { "t": [...], "x": [...] }   // optional (alias: "inverse")
 * }
 * Also supports flat JSON files where foam or capillary keys are placed at the root.
 */
router.post('/upload', (req, res) => {
  if (state.status === 'running') {
    return res.status(400).json({ error: 'Detenga la simulación antes de cargar otro archivo' });
  }
  let doc = req.body?.content;
  if (typeof doc === 'string') {
    try {
      doc = JSON.parse(doc.replace(/^\uFEFF/, ''));
    } catch (e) {
      return res.status(400).json({ error: `JSON inválido: ${e.message}` });
    }
  }
  if (!doc || typeof doc !== 'object' || Array.isArray(doc)) {
    return res.status(400).json({ error: 'El archivo debe contener un objeto JSON' });
  }

  let rawFoam = doc.foam;
  let rawCapillary = doc.capillary;
  let rawMeas = doc.measurements || doc.inverse;

  // Support flat JSON files if the user did not nest under "foam" / "capillary"
  if (!rawFoam && !rawCapillary) {
    const foamKeysFound = Object.keys(doc).filter((k) => k in state.defaults);
    const capKeysFound = Object.keys(doc).filter((k) => CAPILLARY_KEYS.includes(k));
    if (foamKeysFound.length > 0) {
      rawFoam = {};
      foamKeysFound.forEach((k) => { rawFoam[k] = doc[k]; });
    }
    if (capKeysFound.length > 0) {
      rawCapillary = {};
      capKeysFound.forEach((k) => { rawCapillary[k] = doc[k]; });
    }
    if (Array.isArray(doc.t) && Array.isArray(doc.x)) {
      rawMeas = { t: doc.t, x: doc.x };
    }
  }

  if (!rawFoam && !rawCapillary) {
    return res.status(400).json({ error: 'El archivo debe incluir la sección "foam" y/o "capillary" (o sus parámetros)' });
  }

  const uploadedAt = new Date().toISOString();
  const response = {
    valid: true,
    uploadedAt,
    foam: null,
    capillary: null,
    inverse: null,
    report: [],
    errors: [],
    warnings: [],
  };

  let foamValid = false;
  if (rawFoam) {
    const cleanFoam = coerceNumbers(rawFoam);
    const unknown = Object.keys(cleanFoam).filter((k) => !(k in state.defaults));
    unknown.forEach((k) => response.warnings.push({ key: k, level: 'caution', message: `Parámetro desconocido "${k}" ignorado` }));
    const params = {};
    for (const k of Object.keys(state.defaults)) {
      if (k in cleanFoam) params[k] = cleanFoam[k];
    }
    for (const k of REQUIRED_KEYS) {
      if (!(k in params)) {
        params[k] = state.defaults[k];
        response.warnings.push({ key: k, level: 'caution', message: `${k} no está en el archivo: se usó el valor por defecto ${state.defaults[k]}` });
      }
    }
    for (const k of ['zExtractL1', 'zExtractL2']) {
      if (!(k in params)) params[k] = state.defaults[k];
    }
    const report = validateFoamParams(params);
    const s = summarize(report);
    response.report = report;
    response.errors.push(...s.errors);
    response.warnings.push(...s.warnings);
    response.foam = params;
    foamValid = s.errors.length === 0;
  }

  let capValid = false;
  if (rawCapillary) {
    const cleanCap = coerceNumbers(rawCapillary);
    const n = normalizeCapillary(cleanCap);
    n.errors.forEach((m) => response.errors.push({ key: 'capillary', level: 'fail', message: m }));
    response.capillary = cleanCap;
    capValid = n.errors.length === 0;
  }

  if (rawMeas) {
    const tArr = Array.isArray(rawMeas.t) ? rawMeas.t.map(Number) : rawMeas.t;
    const xArr = Array.isArray(rawMeas.x) ? rawMeas.x.map(Number) : rawMeas.x;
    if (!rawCapillary) {
      response.errors.push({ key: 'measurements', level: 'fail', message: 'La sección "measurements" requiere "capillary"' });
    } else {
      response.inverse = { t: tArr, x: xArr };
    }
  }

  response.valid = response.errors.length === 0;

  // Always reset previous simulation frames/history and update state with the newly uploaded file
  resetState();
  if (response.foam) state.params = response.foam;
  state.capillary = capValid ? response.capillary : null;
  state.inverse = capValid ? response.inverse : null;
  state.capFit = null;
  state.capModel = 'none';
  state.paramsLoaded = foamValid;
  state.capillaryLoaded = capValid;
  state.lastReport = { errors: response.errors, warnings: response.warnings, report: response.report };
  state.uploadedAt = uploadedAt;
  state.fileName = req.body?.fileName || null;

  broadcast(req.app?.locals?.wss, {
    type: 'params_updated',
    valid: response.valid,
    uploadedAt,
    fileName: state.fileName,
    foam: response.foam,
    capillary: response.capillary,
    inverse: response.inverse,
    report: state.lastReport,
  });

  // Return 200 even when validation errors exist so the frontend always displays the newly uploaded values
  // and highlights offending parameters, while blocking simulation start via state.paramsLoaded = false.
  res.set('Cache-Control', 'no-store');
  return res.json(response);
});

router.get('/template', (req, res) => {
  res.set('Cache-Control', 'no-store');
  const file = path.resolve(__dirname, '../../examples/parametros_ejemplo.json');
  res.download(file, 'parametros_ejemplo.json');
});

export default router;
