import express from 'express';
import cors from 'cors';
import { WebSocketServer } from 'ws';
import http from 'http';
import { state } from './store/simState.js';
import simulationRoutes from './routes/simulation.js';
import exportRoutes from './routes/export.js';
import paramsRoutes from './routes/params.js';
import capillaryRoutes from './routes/capillary.js';

const app = express();
app.use(cors());
app.use(express.json({ limit: '5mb' }));

app.use('/api/simulate', simulationRoutes);
app.use('/api/export', exportRoutes);
app.use('/api/params', paramsRoutes);
app.use('/api/capillary', capillaryRoutes);

app.get('/api/params/current', (req, res) => {
  res.set('Cache-Control', 'no-store');
  res.json({
    loaded: state.paramsLoaded,
    valid: state.lastReport ? state.lastReport.errors.length === 0 : state.paramsLoaded,
    fileName: state.fileName || null,
    uploadedAt: state.uploadedAt || null,
    params: (state.paramsLoaded || state.uploadedAt) ? state.params : null,
    capillary: state.capillary,
    inverse: state.inverse,
    report: state.lastReport || null,
  });
});

const server = http.createServer(app);
export const wss = new WebSocketServer({ server, path: '/ws' });
app.locals.wss = wss;

wss.on('connection', (ws) => {
  console.log('Client connected to WebSocket');

  ws.send(JSON.stringify({
    type: 'init',
    status: state.status,
    params: (state.paramsLoaded || state.uploadedAt) ? state.params : null,
    capillary: state.capillary,
    inverse: state.inverse,
    progress: state.progress,
  }));

  ws.on('close', () => {
    console.log('Client disconnected');
  });
});

const PORT = process.env.PORT || 3001;
server.listen(PORT, () => {
  console.log(`Server listening on port ${PORT}`);
});
