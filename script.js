// Coupled Waveguide Simulator (Asymmetric 5-Layer Waveguide Solver, Continuous Wave & Auto-Zoom)
const canvas = document.getElementById('sim-canvas');
const ctx = canvas.getContext('2d');
const graphCanvas = document.getElementById('graph-canvas');
const gCtx = graphCanvas.getContext('2d');

// DOM Elements: Sliders & Controls
const sliderGap = document.getElementById('slider-gap');
const sliderWidth = document.getElementById('slider-width');
const sliderLambda = document.getElementById('slider-lambda');
const sliderNcore1 = document.getElementById('slider-ncore1');
const sliderNcore2 = document.getElementById('slider-ncore2');
const sliderNtop = document.getElementById('slider-ntop');
const sliderNmid = document.getElementById('slider-nmid');
const sliderNbot = document.getElementById('slider-nbot');
const sliderSpeed = document.getElementById('slider-speed');

const valSpeed = document.getElementById('val-speed');

const btnPlayPause = document.getElementById('btn-play-pause');
const btnReset = document.getElementById('btn-reset');
const playIcon = btnPlayPause.querySelector('.play-icon');
const pauseIcon = btnPlayPause.querySelector('.pause-icon');

// Simulation Constants & Scale Factors
// Horizontal: x-direction (propagation axis). 800px = 400 um (1 px = 0.5 um, pixelScaleX = 2 px/um)
const pixelScaleX = 2.0; 
let pixelScaleY = 3.0; // Dynamic scale: updated on gap changes for auto-zoom

let isPlaying = true;
let time = 0;
const launchMode = 'wg1'; // Fixed to Waveguide 1 launch

// Physics parameters
let gap = 15.0;       // um (separation d)
let width = 10.0;     // um (core width w)
let lambda = 1.55;    // um (wavelength)
let nCore1 = 1.500;   // Core 1 refractive index
let nCore2 = 1.500;   // Core 2 refractive index
let nTop = 1.450;     // Top cladding index
let nMid = 1.450;     // Middle cladding index
let nBot = 1.450;     // Bottom cladding index
let speed = 1.0;      // Time speed multiplier

// Interactive Profile Probe tracking
let trackerX = 200.0; // Probe position in um along the waveguide

// Register mouse move to slide the vertical probe line and inspect the transverse profile
canvas.addEventListener('mousemove', (e) => {
    const rect = canvas.getBoundingClientRect();
    const x_px = e.clientX - rect.left;
    trackerX = Math.max(0.0, Math.min(400.0, x_px / pixelScaleX));
});

// Solved constants for Asymmetric waveguide layers
let k0, kx1, kx2, gTop, gMid1, gMid2, gBot, phi0, psi0, kappa;
let phi1 = []; // Precomputed transverse profile for WG1
let phi2 = []; // Precomputed transverse profile for WG2

// Calculate the number of guided modes for a 3-layer asymmetric slab waveguide
function getGuidedModeCount(nCore, nCladL, nCladR, width_um, lambda_um) {
    const k0_val = 2 * Math.PI / lambda_um;
    const cladMax = Math.max(nCladL, nCladR);
    const cladMin = Math.min(nCladL, nCladR);
    
    if (nCore <= cladMax) return 0;
    
    // Core transverse wavevector at cutoff (beta = k0 * cladMax)
    const kx_cutoff = k0_val * Math.sqrt(nCore * nCore - cladMax * cladMax);
    
    // Decaying wavevector of the other cladding at cutoff
    const g_other_cutoff = k0_val * Math.sqrt(Math.max(0, cladMax * cladMax - cladMin * cladMin));
    
    // Phase thickness at cutoff
    const phi_cutoff = kx_cutoff * width_um - Math.atan(g_other_cutoff / kx_cutoff);
    
    // Number of modes is 1 + floor(phi_cutoff / pi)
    const modeCount = 1 + Math.floor(phi_cutoff / Math.PI);
    return Math.max(1, modeCount); // At least fundamental mode is guided
}

// Recalculate Physics & Modes
function updatePhysics() {
    // 1. Auto-Zoom Calculation: scale up Y axis when gap is small, scale down when gap is large
    // Ensures details in the gap region are always legible and waveguides don't clip off-screen.
    pixelScaleY = Math.min(28.0, Math.max(2.0, 280 / (width * 2 + gap)));

    // 2. Validation for Waveguide 1 (Core 1 must be denser than adjacent claddings)
    let cladMax1 = Math.max(nTop, nMid);
    if (nCore1 <= cladMax1) {
        nCore1 = cladMax1 + 0.005;
        sliderNcore1.value = nCore1;
        document.getElementById('input-ncore1').value = nCore1.toFixed(3);
    }
    
    // 3. Validation for Waveguide 2 (Core 2 must be denser than adjacent claddings)
    let cladMax2 = Math.max(nMid, nBot);
    if (nCore2 <= cladMax2) {
        nCore2 = cladMax2 + 0.005;
        sliderNcore2.value = nCore2;
        document.getElementById('input-ncore2').value = nCore2.toFixed(3);
    }
    
    k0 = 2 * Math.PI / lambda;
    
    // --- Waveguide 1 Solver (Asymmetric 3-layer Slab: Top, Core1, Mid) ---
    let betaMin1 = k0 * cladMax1 + 0.001;
    let betaMax1 = k0 * nCore1 - 0.001;
    let beta1 = solveBeta(nCore1, nTop, nMid, betaMin1, betaMax1);
    
    kx1 = Math.sqrt(k0*k0*nCore1*nCore1 - beta1*beta1);
    gTop = Math.sqrt(beta1*beta1 - k0*k0*nTop*nTop);
    gMid1 = Math.sqrt(beta1*beta1 - k0*k0*nMid*nMid);
    phi0 = 0.5 * (Math.atan(gTop / kx1) - Math.atan(gMid1 / kx1));
    
    // --- Waveguide 2 Solver (Asymmetric 3-layer Slab: Mid, Core2, Bot) ---
    let betaMin2 = k0 * cladMax2 + 0.001;
    let betaMax2 = k0 * nCore2 - 0.001;
    let beta2 = solveBeta(nCore2, nMid, nBot, betaMin2, betaMax2);
    
    kx2 = Math.sqrt(k0*k0*nCore2*nCore2 - beta2*beta2);
    gMid2 = Math.sqrt(beta2*beta2 - k0*k0*nMid*nMid);
    gBot = Math.sqrt(beta2*beta2 - k0*k0*nBot*nBot);
    psi0 = 0.5 * (Math.atan(gMid2 / kx2) - Math.atan(gBot / kx2));
    
    // Coupled mode propagation terms
    let betaAvg = (beta1 + beta2) / 2;
    let delta = (beta1 - beta2) / 2; // Phase mismatch
    
    // Coupling coefficient kappa
    let gMidAvg = (gMid1 + gMid2) / 2;
    const denominator = Math.sqrt(beta1 * beta2) * (width / 2 + 1 / gTop + 1 / gMid1) * (width / 2 + 1 / gBot + 1 / gMid2);
    kappa = (2 * kx1 * kx2 * gMidAvg * Math.exp(-gMidAvg * gap)) / denominator;
    
    // Store variables globally for use in E-field and Power Graph calculations
    window.g_beta1 = beta1;
    window.g_beta2 = beta2;
    window.g_delta = delta;
    window.g_q = Math.sqrt(kappa * kappa + delta * delta);
    
    // Mode profile normalizations (normalizes core peaks so they render consistently)
    let Ic1 = width/2 + Math.sin(kx1 * width) * Math.cos(2 * phi0) / (2 * kx1);
    let It1 = Math.pow(Math.cos(kx1 * width/2 + phi0), 2) / (2 * gTop);
    let Ib1 = Math.pow(Math.cos(kx1 * width/2 - phi0), 2) / (2 * gMid1);
    let norm1 = 1 / Math.sqrt(Ic1 + It1 + Ib1);

    let Ic2 = width/2 + Math.sin(kx2 * width) * Math.cos(2 * psi0) / (2 * kx2);
    let It2 = Math.pow(Math.cos(kx2 * width/2 + psi0), 2) / (2 * gMid2);
    let Ib2 = Math.pow(Math.cos(kx2 * width/2 - psi0), 2) / (2 * gBot);
    let norm2 = 1 / Math.sqrt(Ic2 + It2 + Ib2);

    // Calculate and display the number of guided modes
    let m1 = getGuidedModeCount(nCore1, nTop, nMid, width, lambda);
    let m2 = getGuidedModeCount(nCore2, nMid, nBot, width, lambda);
    document.getElementById('mode-count-1').textContent = `${m1}개`;
    document.getElementById('mode-count-2').textContent = `${m2}개`;

    precomputeTransverseModes(norm1, norm2);
}

// Bisection search solver for beta (propagation constant)
function solveBeta(nCore, nCladL, nCladR, low, high) {
    let mid = (low + high) / 2;
    for (let i = 0; i < 20; i++) {
        mid = (low + high) / 2;
        let kx = Math.sqrt(k0*k0*nCore*nCore - mid*mid);
        let gL = Math.sqrt(mid*mid - k0*k0*nCladL*nCladL);
        let gR = Math.sqrt(mid*mid - k0*k0*nCladR*nCladR);
        
        let val = kx * width - Math.atan(gL / kx) - Math.atan(gR / kx);
        if (val > 0) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return (low + high) / 2;
}

// Precompute transverse mode profiles along Y axis
function precomputeTransverseModes(norm1, norm2) {
    const cy = canvas.height / 2;
    const gap_px = gap * pixelScaleY;
    const w_px = width * pixelScaleY;
    
    // Core center positions
    const y1 = cy - (gap_px / 2) - (w_px / 2);
    const y2 = cy + (gap_px / 2) + (w_px / 2);
    
    phi1 = new Array(canvas.height);
    phi2 = new Array(canvas.height);
    
    // WG1 mode profile
    for (let y = 0; y < canvas.height; y++) {
        let dist = y - y1;
        let y_um = dist / pixelScaleY;
        if (Math.abs(dist) <= w_px / 2) {
            phi1[y] = norm1 * Math.cos(kx1 * y_um - phi0);
        } else if (dist < -w_px / 2) {
            // top cladding
            let boundary_um = -width / 2;
            let tail_um = -y_um - width / 2;
            phi1[y] = norm1 * Math.cos(kx1 * boundary_um - phi0) * Math.exp(-gTop * tail_um);
        } else {
            // middle cladding
            let boundary_um = width / 2;
            let tail_um = y_um - width / 2;
            phi1[y] = norm1 * Math.cos(kx1 * boundary_um - phi0) * Math.exp(-gMid1 * tail_um);
        }
    }
    
    // WG2 mode profile
    for (let y = 0; y < canvas.height; y++) {
        let dist = y - y2;
        let y_um = dist / pixelScaleY;
        if (Math.abs(dist) <= w_px / 2) {
            phi2[y] = norm2 * Math.cos(kx2 * y_um - psi0);
        } else if (dist < -w_px / 2) {
            // middle cladding
            let boundary_um = -width / 2;
            let tail_um = -y_um - width / 2;
            phi2[y] = norm2 * Math.cos(kx2 * boundary_um - psi0) * Math.exp(-gMid2 * tail_um);
        } else {
            // bottom cladding
            let boundary_um = width / 2;
            let tail_um = y_um - width / 2;
            phi2[y] = norm2 * Math.cos(kx2 * boundary_um - psi0) * Math.exp(-gBot * tail_um);
        }
    }
}

// Render electric field: 2D Continuous Wave color map + Propagating 1D Transverse E(y) profile curve + Stationary borders
function drawElectricField(theta) {
    const width_px = canvas.width;
    const height_px = canvas.height;
    
    const imgData = ctx.createImageData(width_px, height_px);
    const data = imgData.data;
    
    const w_px = width * pixelScaleY;
    const gap_px = gap * pixelScaleY;
    const cy = height_px / 2;
    const y1 = cy - (gap_px / 2) - (w_px / 2);
    const y2 = cy + (gap_px / 2) + (w_px / 2);
    
    const delta = window.g_delta;
    const q = window.g_q;
    const beta1 = window.g_beta1;
    const beta2 = window.g_beta2;
    
    // 1. Draw 2D Phase Continuous Wave (CW) in the background using color density
    for (let x = 0; x < width_px; x += 2) {
        let x_um = x / pixelScaleX;
        
        let cos_qx = Math.cos(q * x_um);
        let sin_qx = Math.sin(q * x_um);
        let cos_b1x = Math.cos(beta1 * x_um - theta);
        let sin_b1x = Math.sin(beta1 * x_um - theta);
        let sin_b2x = Math.sin(beta2 * x_um - theta);
        
        // Continuous wave coupling amplitudes along x
        let c1 = cos_qx * cos_b1x + (delta / q) * sin_qx * sin_b1x;
        let c2 = (kappa / q) * sin_qx * sin_b2x;
        
        for (let y = 0; y < height_px; y += 2) {
            let E = phi1[y] * c1 + phi2[y] * c2;
            
            let r = 8, g = 12, b = 22; // Base dark blue background
            
            if (E > 0.015) {
                let amt = Math.min(1.0, E) * 230;
                r = Math.round(r + amt * 0.1);
                g = Math.round(g + amt * 0.85);
                b = Math.round(b + amt * 1.0);
            } else if (E < -0.015) {
                let amt = Math.min(1.0, -E) * 230;
                r = Math.round(r + amt * 1.0);
                g = Math.round(g + amt * 0.15);
                b = Math.round(b + amt * 0.7);
            }
            
            draw2x2Block(data, x, y, width_px, r, g, b);
        }
    }
    
    ctx.putImageData(imgData, 0, 0);
    
    // 2. Draw Waveguide Core Regions - STRICTLY STATIONARY BORDERS
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    // Core 1 borders
    ctx.moveTo(0, y1 - w_px/2); ctx.lineTo(width_px, y1 - w_px/2);
    ctx.moveTo(0, y1 + w_px/2); ctx.lineTo(width_px, y1 + w_px/2);
    // Core 2 borders
    ctx.moveTo(0, y2 - w_px/2); ctx.lineTo(width_px, y2 - w_px/2);
    ctx.moveTo(0, y2 + w_px/2); ctx.lineTo(width_px, y2 + w_px/2);
    ctx.stroke();
    
    // 3. Draw Propagating 1D Transverse Electric Field Profile Curve E(y) at the Probe Position
    let trackerX_px = trackerX * pixelScaleX;
    
    // Draw vertical reference dashed line (local axis) at the probe position
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.35)';
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(trackerX_px, 0);
    ctx.lineTo(trackerX_px, height_px);
    ctx.stroke();
    ctx.setLineDash([]);
    
    // Draw the 1D transverse profile shape wiggling horizontally around the probe reference line
    ctx.strokeStyle = '#ffffff'; // Solid white wave line
    ctx.lineWidth = 2.5;
    ctx.shadowColor = '#06b6d4'; // Cyan neon glow
    ctx.shadowBlur = 8;
    ctx.beginPath();
    
    let cos_qx_c = Math.cos(q * trackerX);
    let sin_qx_c = Math.sin(q * trackerX);
    let cos_b1x_c = Math.cos(beta1 * trackerX - theta);
    let sin_b1x_c = Math.sin(beta1 * trackerX - theta);
    let sin_b2x_c = Math.sin(beta2 * trackerX - theta);
    
    let c1_prof = cos_qx_c * cos_b1x_c + (delta / q) * sin_qx_c * sin_b1x_c;
    let c2_prof = (kappa / q) * sin_qx_c * sin_b2x_c;
    
    const profileScale = 50.0; // Profile displacement scaling factor
    
    for (let y = 0; y < height_px; y += 2) {
        let E_prof = phi1[y] * c1_prof + phi2[y] * c2_prof;
        let x_disp = trackerX_px + E_prof * profileScale;
        
        if (y === 0) ctx.moveTo(x_disp, y);
        else ctx.lineTo(x_disp, y);
    }
    ctx.stroke();
    ctx.shadowBlur = 0; // Reset shadow glow
    
    // Core & cladding static labels
    ctx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    ctx.font = 'bold 10px Inter';
    ctx.fillText(`CORE 1 (n = ${nCore1.toFixed(3)})`, 15, y1 + 4);
    ctx.fillText(`CORE 2 (n = ${nCore2.toFixed(3)})`, 15, y2 + 4);
    
    ctx.font = '9px Inter';
    ctx.fillText(`상부 클래딩 (n = ${nTop.toFixed(3)})`, 15, y1 - w_px/2 - 8);
    ctx.fillText(`중간 클래딩 (n = ${nMid.toFixed(3)})`, 15, cy + 4);
    ctx.fillText(`하부 클래딩 (n = ${nBot.toFixed(3)})`, 15, y2 + w_px/2 + 14);
    
    // Draw axis arrows & coordinate system details
    drawAxes();
}

// Draw Coordinate Axis Indicator in the corner of the canvas
function drawAxes() {
    const startX = 40;
    const startY = canvas.height - 40;
    
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.75)';
    ctx.fillStyle = 'rgba(255, 255, 255, 0.75)';
    ctx.lineWidth = 2.0;
    
    // X Axis arrow (propagation along waveguide)
    ctx.beginPath();
    ctx.moveTo(startX, startY);
    ctx.lineTo(startX + 50, startY);
    ctx.lineTo(startX + 44, startY - 4);
    ctx.moveTo(startX + 50, startY);
    ctx.lineTo(startX + 44, startY + 4);
    ctx.stroke();
    
    // Y Axis arrow (transverse across waveguide, pointing upward)
    ctx.beginPath();
    ctx.moveTo(startX, startY);
    ctx.lineTo(startX, startY - 50);
    ctx.lineTo(startX - 4, startY - 44);
    ctx.moveTo(startX, startY - 50);
    ctx.lineTo(startX + 4, startY - 44);
    ctx.stroke();
    
    // Labels
    ctx.font = 'bold 11px Inter';
    ctx.fillText('x (진행 방향)', startX + 55, startY + 4);
    ctx.fillText('y (수직 횡방향)', startX - 12, startY - 56);
}

// Helper to write a 2x2 pixel block into the ImageData array
function draw2x2Block(data, x, y, width, r, g, b) {
    const idx1 = (y * width + x) * 4;
    const idx2 = ((y + 1) * width + x) * 4;
    
    data[idx1] = r; data[idx1+1] = g; data[idx1+2] = b; data[idx1+3] = 255;
    data[idx1+4] = r; data[idx1+5] = g; data[idx1+6] = b; data[idx1+7] = 255;
    
    data[idx2] = r; data[idx2+1] = g; data[idx2+2] = b; data[idx2+3] = 255;
    data[idx2+4] = r; data[idx2+5] = g; data[idx2+6] = b; data[idx2+7] = 255;
}

// Draw Longitudinal Power Density Graph along x axis
function drawPowerGraph() {
    const w = graphCanvas.width;
    const h = graphCanvas.height;
    
    gCtx.fillStyle = '#040711';
    gCtx.fillRect(0, 0, w, h);
    
    // Draw Grid Lines
    gCtx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
    gCtx.lineWidth = 1;
    for (let i = 0; i <= 2; i++) {
        let y_px = h * 0.15 + (i * h * 0.35);
        gCtx.beginPath();
        gCtx.moveTo(40, y_px);
        gCtx.lineTo(w - 20, y_px);
        gCtx.stroke();
    }
    
    for (let x_val = 0; x_val <= 400; x_val += 100) {
        let x_px = 40 + (x_val / 400) * (w - 60);
        gCtx.beginPath();
        gCtx.moveTo(x_px, h * 0.15);
        gCtx.lineTo(x_px, h * 0.85);
        gCtx.stroke();
        
        gCtx.fillStyle = 'rgba(255,255,255,0.3)';
        gCtx.font = '9px Fira Code';
        gCtx.fillText(`${x_val}um`, x_px - 12, h * 0.95);
    }
    
    gCtx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    gCtx.font = '9px Fira Code';
    gCtx.fillText('1.0', 15, h * 0.18);
    gCtx.fillText('0.5', 15, h * 0.53);
    gCtx.fillText('0.0', 15, h * 0.88);
    
    // Plot lines
    gCtx.lineWidth = 2.5;
    
    const delta = window.g_delta;
    const q = window.g_q;
    
    // 1. WG1 Power Curve (Continuous spatial power flow along x)
    gCtx.strokeStyle = '#38bdf8'; // Cyan
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let x_um = ((x - 40) / (w - 60)) * 400;
        let p1 = Math.pow(Math.cos(q * x_um), 2) + Math.pow(delta / q, 2) * Math.pow(Math.sin(q * x_um), 2);
        
        let y_px = h * 0.85 - p1 * (h * 0.7);
        if (x === 40) gCtx.moveTo(x, y_px);
        else gCtx.lineTo(x, y_px);
    }
    gCtx.stroke();
    
    // 2. WG2 Power Curve
    gCtx.strokeStyle = '#f43f5e'; // Pink
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let x_um = ((x - 40) / (w - 60)) * 400;
        let p2 = Math.pow(kappa / q, 2) * Math.pow(Math.sin(q * x_um), 2);
        
        let y_px = h * 0.85 - p2 * (h * 0.7);
        if (x === 40) gCtx.moveTo(x, y_px);
        else gCtx.lineTo(x, y_px);
    }
    gCtx.stroke();
    
    // 3. Draw Vertical Probe Position Indicator in the Power Graph
    let trackerX_graph_px = 40 + (trackerX / 400) * (w - 60);
    gCtx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
    gCtx.lineWidth = 1;
    gCtx.setLineDash([2, 2]);
    gCtx.beginPath();
    gCtx.moveTo(trackerX_graph_px, h * 0.15);
    gCtx.lineTo(trackerX_graph_px, h * 0.85);
    gCtx.stroke();
    gCtx.setLineDash([]);
    
    // Legend labels
    gCtx.font = 'bold 9px Inter';
    gCtx.fillStyle = '#38bdf8';
    gCtx.fillText('도파로 1 전력 (WG1 Power)', w - 160, h * 0.25);
    gCtx.fillStyle = '#f43f5e';
    gCtx.fillText('도파로 2 전력 (WG2 Power)', w - 160, h * 0.38);
}

// Main Simulation Loop
function tick() {
    if (isPlaying) {
        time += 0.05 * speed;
    }
    
    drawElectricField(time);
    drawPowerGraph();
    
    requestAnimationFrame(tick);
}

// Helper function to bind sliders with numeric input fields and apply buttons
function bindParam(sliderId, inputId, applyBtnId, minVal, maxVal, isInt, updateVarFn) {
    const slider = document.getElementById(sliderId);
    const input = document.getElementById(inputId);
    const btn = document.getElementById(applyBtnId);
    
    // Sync slider changes to numerical text field instantly during dragging
    slider.addEventListener('input', (e) => {
        let val = parseFloat(e.target.value);
        input.value = isInt ? val.toFixed(0) : val.toFixed(3);
        updateVarFn(val);
        updatePhysics();
    });
    
    // Apply typed value on button click or Enter keypress
    function applyManualValue() {
        let val = parseFloat(input.value);
        if (isNaN(val)) return;
        
        // Guard boundaries
        val = Math.max(minVal, Math.min(maxVal, val));
        input.value = isInt ? val.toFixed(0) : val.toFixed(3);
        
        // Sync slider position if it lies within the slider's track bounds
        let sliderMin = parseFloat(slider.min);
        let sliderMax = parseFloat(slider.max);
        slider.value = Math.max(sliderMin, Math.min(sliderMax, val));
        
        updateVarFn(val);
        updatePhysics();
    }
    
    btn.addEventListener('click', applyManualValue);
    
    input.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') {
            applyManualValue();
        }
    });
}

// Slider event listeners & display updates
function setupSliders() {
    bindParam('slider-gap', 'input-gap', 'btn-apply-gap', 1.0, 100.0, false, (v) => { gap = v; });
    bindParam('slider-width', 'input-width', 'btn-apply-width', 0.1, 100.0, false, (v) => { width = v; });
    bindParam('slider-lambda', 'input-lambda', 'btn-apply-lambda', 0.1, 10.0, false, (v) => { lambda = v; });
    bindParam('slider-ncore1', 'input-ncore1', 'btn-apply-ncore1', 1.0, 4.0, false, (v) => { nCore1 = v; });
    bindParam('slider-ncore2', 'input-ncore2', 'btn-apply-ncore2', 1.0, 4.0, false, (v) => { nCore2 = v; });
    bindParam('slider-ntop', 'input-ntop', 'btn-apply-ntop', 1.0, 4.0, false, (v) => { nTop = v; });
    bindParam('slider-nmid', 'input-nmid', 'btn-apply-nmid', 1.0, 4.0, false, (v) => { nMid = v; });
    bindParam('slider-nbot', 'input-nbot', 'btn-apply-nbot', 1.0, 4.0, false, (v) => { nBot = v; });
    
    // Speed slider
    sliderSpeed.addEventListener('input', (e) => {
        speed = parseFloat(e.target.value);
        valSpeed.textContent = `${speed.toFixed(1)}x`;
    });
}

function setupControls() {
    btnPlayPause.addEventListener('click', () => {
        isPlaying = !isPlaying;
        if (isPlaying) {
            playIcon.style.display = 'none';
            pauseIcon.style.display = 'inline';
        } else {
            playIcon.style.display = 'inline';
            pauseIcon.style.display = 'none';
        }
    });
    
    btnReset.addEventListener('click', () => {
        time = 0;
        sliderGap.value = 15; document.getElementById('input-gap').value = '15.000'; gap = 15;
        sliderWidth.value = 10; document.getElementById('input-width').value = '10.000'; width = 10;
        sliderLambda.value = 1.55; document.getElementById('input-lambda').value = '1.550'; lambda = 1.55;
        sliderNcore1.value = 1.500; document.getElementById('input-ncore1').value = '1.500'; nCore1 = 1.500;
        sliderNcore2.value = 1.500; document.getElementById('input-ncore2').value = '1.500'; nCore2 = 1.500;
        sliderNtop.value = 1.450; document.getElementById('input-ntop').value = '1.450'; nTop = 1.450;
        sliderNmid.value = 1.450; document.getElementById('input-nmid').value = '1.450'; nMid = 1.450;
        sliderNbot.value = 1.450; document.getElementById('input-nbot').value = '1.450'; nBot = 1.450;
        sliderSpeed.value = 1.0; valSpeed.textContent = '1.0x'; speed = 1.0;
        
        isPlaying = true;
        playIcon.style.display = 'none';
        pauseIcon.style.display = 'inline';
        
        trackerX = 200.0;
        
        updatePhysics();
    });
}

// Initialise Application
updatePhysics();
setupSliders();
setupControls();

// Start Sim loop
playIcon.style.display = 'none';
pauseIcon.style.display = 'inline';
tick();
