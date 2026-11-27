import express from 'express';
import { state } from '../store/simState.js';
import { normalizeCapillary, simulateCapillary, inverseCapillary } from '../capillary/lwb.js';

const router = express.Router();

/** Forward LWB simulation. Uses the uploaded file's "capillary" section unless a body is given. */
router.post('/simulate', (req, res) => {
  const raw = req.body?.capillary ?? state.capillary;
  if (!raw) return res.status(400).json({ error: 'No hay sección "capillary" cargada' });
  const n = normalizeCapillary(raw);
  if (n.errors.length) return res.status(422).json({ errors: n.errors });
  try {
    res.json(simulateCapillary(n.value));
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

/** Inverse capillary: fit unknown parameters from measured x(t). */
router.post('/inverse', (req, res) => {
  const raw = req.body?.capillary ?? state.capillary;
  const inv = req.body?.inverse ?? state.inverse;
  if (!raw || !inv) return res.status(400).json({ error: 'Faltan las secciones "capillary" y "measurements"' });
  const n = normalizeCapillary(raw);
  if (n.errors.length) return res.status(422).json({ errors: n.errors });
  try {
    const out = inverseCapillary(n.value, inv, req.body?.fit || inv.fit || ['r']);
    if (out.errors?.length) return res.status(422).json({ errors: out.errors });
    state.capFit = out.fitted; // used by the foam simulation in 'inverse' coupling mode
    res.json(out);
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

export default router;
