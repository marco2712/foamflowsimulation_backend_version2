# Foam Flow Simulation — Backend (Version 2)

Node.js / Express / WebSocket simulation backend for **2D Stratified Porous-Media Foam Flow (Local Equilibrium)** coupled with **Lucas–Washburn–Bosanquet (LWB)** and **Inverse Capillary** models.

---

## Summary of Changes in Version 2

1. **File-Based Parameter Ingestion & Operating-Limit Validation (`src/routes/params.js`, `src/validation/limits.js`)**
   - Replaced manual parameter entry with structured JSON parameter file uploads (`POST /api/params/upload`) and a downloadable reference template (`GET /api/params/template` -> `examples/parametros_ejemplo.json`).
   - Supports both structured JSON (`foam`, `capillary`, `measurements`) and flat root-level parameter files, including automatic numeric coercion and UTF-8 BOM handling.
   - Validates physical and computational parameters against the operating limits specification (`ok`, `caution`, and blocking `fail` levels covering saturations $S_w^-$, $S_w^+$, $S_w^*$, foam transition slope $A$, permeability contrast $k_1/k_2$, porosities $\phi_1,\phi_2$, injection velocities $u_1,u_2$, viscosity ratio $\mu_w/\mu_g$, mesh size $N_x \times N_z$, and advective CFL).
   - Always returns parsed values and detailed diagnostics so the frontend can display newly uploaded values and highlight invalid keys while preventing invalid simulations from starting (`state.paramsLoaded = false` when blocking errors exist).
   - Broadcasts `params_updated` events over WebSocket and resets server-side simulation buffers upon upload.

2. **Lucas–Washburn–Bosanquet (LWB) Forward Solver (`src/capillary/lwb.js`, `POST /api/capillary/simulate`)**
   - Solves the full second-order ODE for capillary meniscus penetration $x(t)$ in a cylindrical tube of radius $r$ (SI units), including optional hydrostatic gravity for inclined/vertical tubes ($\beta$):
     $$\frac{d}{dt}\left(\pi r^2 \rho x \frac{dx}{dt}\right) + 8\pi \eta x \frac{dx}{dt} = 2\pi r \gamma \cos\theta - \pi r^2 \rho g x \sin\beta$$
   - **Singularity handling at $t = 0$:** Starts integration at $x_0 = 10^{-9}\text{ m}$ with the exact analytical early-time inertial velocity $v_0 = \sqrt{\frac{2\gamma\cos\theta}{\rho r}}$ at $t_0 = x_0 / v_0$ (using $x(t) = v_0 t$ for $t \le t_0$).
   - **Integrator:** Adaptive Dormand–Prince Runge–Kutta 5(4) (`RK45`) with relative/absolute error control.
   - **Validation outputs:** Computes classical diffusive Lucas–Washburn $x_{LW}(t) = \sqrt{\frac{\gamma r \cos\theta}{2\eta} t}$, local log–log regime slope $\frac{d\ln x}{d\ln t}$ (verifying the transition from $1.0$ inertial to $0.5$ viscous around $\tau = \frac{\rho r^2}{8\eta}$), meniscus velocity $v(t)$, and Jurin equilibrium height when gravity is active.

3. **Inverse Capillary Parameter Identification (`src/capillary/lwb.js`, `POST /api/capillary/inverse`)**
   - Identifies unknown capillary parameters (any subset of `r`, `eta`, `gamma`, `theta_deg` chosen at runtime from the UI) from experimental/synthetic penetration data `measurements: { t: [...], x: [...] }`.
   - Uses constrained parameter transformations (log-space for positive quantities and logit-space for $\theta \in (0^\circ, 89.9^\circ)$) with Nelder–Mead simplex optimization minimizing $\sum_i (\ln x_{\text{model}}(t_i) - \ln x_{\text{meas},i})^2$ across inertial and viscous regimes.
   - Stores fitted parameters in server state (`state.capFit`) so they can be coupled directly into the 2D foam simulation.

4. **Capillary-to-Foam 2D Coupling (`src/capillary/coupling.js`, `src/routes/simulation.js`)**
   - Connects the selected capillary model (`none`, `lwb`, or `inverse`) to the 2D Leverett capillary pressure and diffusion operator $D_{\text{cap}}(S_w, n_D) = -\lambda_g f_w \frac{dP_c}{dS_w}$:
     - **`none`**: Uses the baseline `foam.sigma` parameter directly.
     - **`lwb`**: Scales Leverett capillary pressure using effective wetting tension $\sigma_{\text{eff}} = \gamma \cos\theta$ and modulates $D_{\text{cap}}$ dynamically by $f(t) = \min\!\left(1, \frac{v_{\text{LWB}}(t)}{v_{\text{LW}}(t)}\right)$ to capture early-time inertial capillary delay.
     - **`inverse`**: Applies $\sigma_{\text{eff}} = \gamma_{\text{fit}} \cos\theta_{\text{fit}}$ calibrated from the inverse capillary fit.

5. **Automatic Nonlinear CFL + 2D Parabolic Sub-Stepping & Fast LUT Acceleration (`src/simulation/worker.js`)**
   - Prevents odd–even sawtooth numerical oscillations in $S_w$ and foam texture $n_D$ on fine meshes (e.g., short domains $L = 0.25\text{ m}$) or high capillary pressures.
   - Precomputes 4096-point linear-interpolation lookup tables (`Float64Array`) per layer for $f_w(S_w)$, $D_{\text{cap}}(S_w)$, $n_{D,\text{LE}}(S_w)$, and $\left|\frac{df_w}{dS_w}\right|$.
   - Dynamically subdivides each macro time step $\Delta t$ into $N_{\text{sub}}$ sub-steps so the combined nonlinear Buckley–Leverett advective speed and 2D explicit capillary diffusion satisfy $\text{CFL}_{\text{total}} \le 0.65$ at all times.
   - Honors `zExtractL1` and `zExtractL2` ratios when extracting 1D layer profiles.

6. **Developer Workflow (`package.json`)**
   - Configured `npm start` to run `node --watch src/index.js` for automatic backend reloading during development.

---

## Project Structure

```text
src/
├── index.js                  # Express HTTP server + WebSocket (/ws) setup (default port 3001)
├── capillary/
│   ├── lwb.js                # LWB RK45 ODE solver, analytical LW, and Nelder–Mead inverse solver
│   └── coupling.js           # Static & dynamic capillary coupling for the 2D foam solver
├── routes/
│   ├── params.js             # Parameter JSON upload, validation, and template download
│   ├── capillary.js          # Forward LWB and inverse capillary REST endpoints
│   ├── simulation.js         # Start / pause / resume / stop worker thread & WS broadcast
│   └── export.js             # CSV export of time-series history
├── simulation/
│   ├── constitutive.js       # Corey relative permeabilities, MRF, nD_LE, Pc, and D_cap
│   ├── operators.js          # 1st-order upwind advection, harmonic-mean 2D diffusion, BCs
│   └── worker.js             # Worker-thread solver with LUT acceleration & adaptive sub-stepping
├── store/
│   └── simState.js           # In-memory simulation state and immutable fallback defaults
└── validation/
    └── limits.js             # Physical and numerical operating-limit rules
examples/
└── parametros_ejemplo.json   # Complete sample input file (foam + capillary + measurements)
```

---

## Getting Started

```bash
npm install
npm start       # Starts server on http://localhost:3001 with --watch
```

To run without file watching in production:

```bash
npm run serve
```

---

## REST & WebSocket API

| Method & Path | Description |
|---|---|
| `POST /api/params/upload` | Body: `{ content, fileName? }`. Parses & validates JSON parameters and broadcasts `params_updated`. |
| `GET /api/params/current` | Returns currently loaded `params`, `capillary`, `inverse` measurements, and validation report. |
| `GET /api/params/template` | Downloads `examples/parametros_ejemplo.json`. |
| `POST /api/capillary/simulate` | Runs the forward LWB vs. Lucas–Washburn ODE solver on the loaded `capillary` parameters. |
| `POST /api/capillary/inverse` | Body: `{ fit?: ["r" \| "eta" \| "gamma" \| "theta_deg"] }`. Runs inverse parameter identification against `measurements`. |
| `POST /api/simulate/start` | Body: `{ capModel?: "none" \| "lwb" \| "inverse" }`. Starts the 2D foam worker thread. |
| `POST /api/simulate/pause` | Pauses the running 2D simulation worker. |
| `POST /api/simulate/resume` | Resumes a paused 2D simulation worker. |
| `POST /api/simulate/stop` | Stops and terminates the active 2D simulation worker. |
| `GET /api/export/csv` | Exports accumulated 2D front positions, crossflow, and recovery history as CSV. |
| `WS /ws` | Streams `init`, `status`, `frame`, `done`, and `params_updated` messages in real time. |

---

## Parameter File Format

```json
{
  "foam": {
    "phi1": 0.25, "phi2": 0.25,
    "Swc": 0.20, "Sgr": 0.18,
    "mu_w": 0.001, "mu_g": 0.00002,
    "k1": 2e-12, "k2": 1e-12,
    "Sw_star": 0.37, "A": 400, "Kc": 200, "theta_s": 0.00032,
    "u1": 0.0000293, "u2": 0.00001465,
    "d": 0.005, "sigma": 0.03, "c_cap": 0.01,
    "Sw_minus": 0.372, "Sw_plus": 0.72,
    "L": 1.0, "Nx": 200, "Nz": 40,
    "dt": 0.2, "Tmax": 5000,
    "zExtractL1": 0.5, "zExtractL2": 0.5
  },
  "capillary": {
    "rho": 998, "eta": 0.001, "gamma": 0.0728,
    "theta_deg": 0, "r": 0.0005, "beta_deg": 0,
    "tmax": 5, "n_points": 200
  },
  "measurements": {
    "t": [0.001, 0.002, 0.004, 0.008],
    "x": [0.000599, 0.001186, 0.002328, 0.004502]
  }
}
```
