import express from 'express';
import { Worker } from 'worker_threads';
import path from 'path';
import { fileURLToPath } from 'url';
import { state, resetState } from '../store/simState.js';
import { wss } from '../index.js';
import WebSocket from 'ws';
import { buildCoupling, CAP_MODES } from '../capillary/coupling.js';
import { normalizeCapillary } from '../capillary/lwb.js';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const router = express.Router();

let worker = null;

function broadcast(data) {
  const json = JSON.stringify(data);
  wss.clients.forEach((client) => {
    if (client.readyState === WebSocket.OPEN) {
      client.send(json);
    }
  });
}

router.post('/start', (req, res) => {
  if (!state.paramsLoaded) {
    return res.status(400).json({ error: 'Cargue primero un archivo de parámetros válido (sin errores en la sección "foam")' });
  }
  if (state.status === 'running') {
    return res.status(400).json({ error: 'Simulation already running' });
  }

  const capModel = req.body?.capModel || 'none';
  if (!CAP_MODES.includes(capModel)) return res.status(400).json({ error: `capModel inválido: ${capModel}` });
  let coupling;
  try {
    let capNorm = null;
    if (capModel !== 'none') {
      if (!state.capillary) return res.status(400).json({ error: 'El archivo no tiene sección "capillary" válida' });
      capNorm = normalizeCapillary(state.capillary).value;
    }
    if (capModel === 'inverse' && !state.capFit) {
      return res.status(400).json({ error: 'Ejecute primero el ajuste de capilaridad inversa' });
    }
    coupling = buildCoupling(capModel, capNorm, state.capFit, state.params, state.params.Tmax);
  } catch (e) {
    return res.status(500).json({ error: e.message });
  }
  state.capModel = capModel;
  resetState();
  state.status = 'running';
  broadcast({ type: 'status', status: 'running' });

  worker = new Worker(path.resolve(__dirname, '../simulation/worker.js'), {
    workerData: { params: state.params, cap: coupling },
  });

  worker.on('message', (msg) => {
    if (msg.type === 'frame') {
      const { t, metricsBatch } = msg;

      if (metricsBatch && metricsBatch.hist_t.length > 0) {
        state.hist_t.push(...metricsBatch.hist_t);
        state.hist_fp1.push(...metricsBatch.hist_fp1);
        state.hist_fp2.push(...metricsBatch.hist_fp2);
        state.hist_fp1_nD.push(...metricsBatch.hist_fp1_nD);
        state.hist_fp2_nD.push(...metricsBatch.hist_fp2_nD);
        state.hist_trans_z.push(...metricsBatch.hist_trans_z);
        state.hist_rec_pct.push(...metricsBatch.hist_rec_pct);
      }

      state.lastFrame = { Sw: msg.Sw, nD: msg.nD };
      state.progress = t / state.params.Tmax;
      broadcast(msg);
    } else if (msg.type === 'done') {
      state.status = 'done';
      worker = null;
      broadcast({ type: 'done', summary: msg.summary });
    }
  });

  worker.on('error', (err) => {
    console.error('Worker error:', err);
    state.status = 'idle';
    broadcast({ type: 'status', status: 'idle' });
    worker = null;
  });

  worker.on('exit', (code) => {
    if (code !== 0) {
      console.error('Worker exited with code', code);
      state.status = 'idle';
      broadcast({ type: 'status', status: 'idle' });
      worker = null;
    }
  });

  res.json({
    message: 'Simulation started',
    params: state.params,
    coupling: { mode: coupling.mode, sigmaEff: coupling.sigmaEff, info: coupling.info },
  });
});

router.post('/pause', (req, res) => {
  if (worker && state.status === 'running') {
    worker.postMessage({ command: 'pause' });
    state.status = 'paused';
    broadcast({ type: 'status', status: 'paused' });
    res.json({ message: 'Simulation paused' });
  } else {
    res.status(400).json({ error: 'Simulation not running' });
  }
});

router.post('/resume', (req, res) => {
  if (worker && state.status === 'paused') {
    worker.postMessage({ command: 'resume' });
    state.status = 'running';
    broadcast({ type: 'status', status: 'running' });
    res.json({ message: 'Simulation resumed' });
  } else {
    res.status(400).json({ error: 'Simulation not paused' });
  }
});

router.post('/stop', (req, res) => {
  if (worker) {
    worker.postMessage({ command: 'stop' });
    setTimeout(() => {
      if (worker) {
        worker.terminate();
        worker = null;
      }
    }, 100);
    state.status = 'idle';
    broadcast({ type: 'status', status: 'idle' });
    worker = null;
    res.json({ message: 'Simulation stopped' });
  } else {
    state.status = 'idle';
    broadcast({ type: 'status', status: 'idle' });
    res.json({ message: 'No active simulation, reset to idle' });
  }
});

export default router;
