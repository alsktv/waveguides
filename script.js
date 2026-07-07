// Coupled Waveguide Simulator (Asymmetric 5-Layer Waveguide Solver, 2D Field & Propagating 1D Mode Profile)
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

const valGap = document.getElementById('val-gap');
const valWidth = document.getElementById('val-width');
const valLambda = document.getElementById('val-lambda');
const valNcore1 = document.getElementById('val-ncore1');
const valNcore2 = document.getElementById('val-ncore2');
const valNtop = document.getElementById('val-ntop');
const valNmid = document.getElementById('val-nmid');
const valNbot = document.getElementById('val-nbot');
const valSpeed = document.getElementById('val-speed');

const btnPlayPause = document.getElementById('btn-play-pause');
const btnReset = document.getElementById('btn-reset');
const playIcon = btnPlayPause.querySelector('.play-icon');
const pauseIcon = btnPlayPause.querySelector('.pause-icon');

// Simulation Constants & Scale Factors
// Horizontal: x-direction (propagation axis). 800px = 400 um (1 px = 0.5 um, pixelScaleX = 2 px/um)
// Vertical: y-direction (transverse axis). 360px = 120 um (1 px = 0.33 um, pixelScaleY = 3 px/um)
const pixelScaleX = 2.0; 
const pixelScaleY = 3.0;

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

// Pulse parameters
const sigmaPulse = 25.0; // Pulse width in um

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
    // 1. Validation for Waveguide 1 (Core 1 must be denser than adjacent claddings)
    let cladMax1 = Math.max(nTop, nMid);
    if (nCore1 <= cladMax1) {
        nCore1 = cladMax1 + 0.005;
        sliderNcore1.value = nCore1;
        valNcore1.textContent = nCore1.toFixed(3);
    }
    
    // 2. Validation for Waveguide 2 (Core 2 must be denser than adjacent claddings)
    let cladMax2 = Math.max(nMid, nBot);
    if (nCore2 <= cladMax2) {
        nCore2 = cladMax2 + 0.005;
        sliderNcore2.value = nCore2;
        valNcore2.textContent = nCore2.toFixed(3);
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
    
    // Coupling coefficient kappa (calculated with average cladding decay in the mid region)
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

// Render electric field: 2D Color Density map + Moving 1D Transverse profile curve + Stationary borders
function drawElectricField(theta, pulseCenter) {
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
    
    // 1. Draw 2D Phase Wave Packet in the background using color density
    for (let x = 0; x < width_px; x += 2) {
        let x_um = x / pixelScaleX;
        
        let distToCenter = x_um - pulseCenter;
        let gauss = Math.exp(- (distToCenter * distToCenter) / (2 * sigmaPulse * sigmaPulse));
        
        if (gauss < 0.001) {
            for (let y = 0; y < height_px; y += 2) {
                draw2x2Block(data, x, y, width_px, 8, 12, 22);
            }
            continue;
        }
        
        let cos_qx = Math.cos(q * x_um);
        let sin_qx = Math.sin(q * x_um);
        let cos_b1x = Math.cos(beta1 * x_um - theta);
        let sin_b1x = Math.sin(beta1 * x_um - theta);
        let sin_b2x = Math.sin(beta2 * x_um - theta);
        
        // Asymmetric Coupled Mode theory amplitudes
        let c1 = gauss * (cos_qx * cos_b1x + (delta / q) * sin_qx * sin_b1x);
        let c2 = gauss * (kappa / q) * sin_qx * sin_b2x;
        
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
    
    // 3. Draw Propagating 1D Transverse Electric Field Profile Curve E(y) at the Pulse Center
    let pulseCenter_px = pulseCenter * pixelScaleX;
    
    // Draw vertical reference dashed line (local axis) at the center of the propagating pulse
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(pulseCenter_px, 0);
    ctx.lineTo(pulseCenter_px, height_px);
    ctx.stroke();
    ctx.setLineDash([]);
    
    // Draw the 1D transverse profile shape wiggling horizontally around the reference line
    ctx.strokeStyle = '#ffffff'; // Solid white wave line
    ctx.lineWidth = 2.5;
    ctx.shadowColor = '#06b6d4'; // Cyan neon glow
    ctx.shadowBlur = 8;
    ctx.beginPath();
    
    let cos_qx_c = Math.cos(q * pulseCenter);
    let sin_qx_c = Math.sin(q * pulseCenter);
    let cos_b1x_c = Math.cos(beta1 * pulseCenter - theta);
    let sin_b1x_c = Math.sin(beta1 * pulseCenter - theta);
    let sin_b2x_c = Math.sin(beta2 * pulseCenter - theta);
    
    let c1_prof = cos_qx_c * cos_b1x_c + (delta / q) * sin_qx_c * sin_b1x_c;
    let c2_prof = (kappa / q) * sin_qx_c * sin_b2x_c;
    
    const profileScale = 50.0; // Profile displacement scaling factor
    
    for (let y = 0; y < height_px; y += 2) {
        let E_prof = phi1[y] * c1_prof + phi2[y] * c2_prof;
        let x_disp = pulseCenter_px + E_prof * profileScale;
        
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
    const startX = 35;
    const startY = canvas.height - 35;
    
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)';
    ctx.fillStyle = 'rgba(255, 255, 255, 0.5)';
    ctx.lineWidth = 1.5;
    
    // X Axis arrow (propagation along waveguide)
    ctx.beginPath();
    ctx.moveTo(startX, startY);
    ctx.lineTo(startX + 40, startY);
    ctx.lineTo(startX + 35, startY - 3);
    ctx.moveTo(startX + 40, startY);
    ctx.lineTo(startX + 35, startY + 3);
    ctx.stroke();
    
    // Y Axis arrow (transverse across waveguide, pointing upward)
    ctx.beginPath();
    ctx.moveTo(startX, startY);
    ctx.lineTo(startX, startY - 40);
    ctx.lineTo(startX - 3, startY - 35);
    ctx.moveTo(startX, startY - 40);
    ctx.lineTo(startX + 3, startY - 35);
    ctx.stroke();
    
    // Labels
    ctx.font = 'italic 11px Inter';
    ctx.fillText('x (진행 방향)', startX + 45, startY + 4);
    ctx.fillText('y (수직 횡방향)', startX - 10, startY - 48);
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
function drawPowerGraph(pulseCenter) {
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
    
    // 1. WG1 Power Curve
    gCtx.strokeStyle = '#38bdf8'; // Cyan
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let x_um = ((x - 40) / (w - 60)) * 400;
        let distToCenter = x_um - pulseCenter;
        let gauss = Math.exp(- (distToCenter * distToCenter) / (2 * sigmaPulse * sigmaPulse));
        let p1 = gauss * gauss * (Math.pow(Math.cos(q * x_um), 2) + Math.pow(delta / q, 2) * Math.pow(Math.sin(q * x_um), 2));
        
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
        let distToCenter = x_um - pulseCenter;
        let gauss = Math.exp(- (distToCenter * distToCenter) / (2 * sigmaPulse * sigmaPulse));
        let p2 = gauss * gauss * Math.pow(kappa / q, 2) * Math.pow(Math.sin(q * x_um), 2);
        
        let y_px = h * 0.85 - p2 * (h * 0.7);
        if (x === 40) gCtx.moveTo(x, y_px);
        else gCtx.lineTo(x, y_px);
    }
    gCtx.stroke();
    
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
    
    // Group velocity movement for Gaussian pulse center (loops along x propagation direction)
    let pulseCenter = (time * 12) % (400 + 160) - 80;
    
    drawElectricField(time, pulseCenter);
    drawPowerGraph(pulseCenter);
    
    requestAnimationFrame(tick);
}

// Slider event listeners & display updates
function setupSliders() {
    sliderGap.addEventListener('input', (e) => {
        gap = parseFloat(e.target.value);
        valGap.textContent = `${gap} μm`;
        updatePhysics();
    });
    
    sliderWidth.addEventListener('input', (e) => {
        width = parseFloat(e.target.value);
        valWidth.textContent = `${width} μm`;
        updatePhysics();
    });
    
    sliderLambda.addEventListener('input', (e) => {
        lambda = parseFloat(e.target.value);
        valLambda.textContent = `${lambda.toFixed(2)} μm`;
        updatePhysics();
    });
    
    sliderNcore1.addEventListener('input', (e) => {
        nCore1 = parseFloat(e.target.value);
        valNcore1.textContent = nCore1.toFixed(3);
        updatePhysics();
    });

    sliderNcore2.addEventListener('input', (e) => {
        nCore2 = parseFloat(e.target.value);
        valNcore2.textContent = nCore2.toFixed(3);
        updatePhysics();
    });

    sliderNtop.addEventListener('input', (e) => {
        nTop = parseFloat(e.target.value);
        valNtop.textContent = nTop.toFixed(3);
        updatePhysics();
    });

    sliderNmid.addEventListener('input', (e) => {
        nMid = parseFloat(e.target.value);
        valNmid.textContent = nMid.toFixed(3);
        updatePhysics();
    });

    sliderNbot.addEventListener('input', (e) => {
        nBot = parseFloat(e.target.value);
        valNbot.textContent = nBot.toFixed(3);
        updatePhysics();
    });
    
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
        sliderGap.value = 15; valGap.textContent = '15 μm'; gap = 15;
        sliderWidth.value = 10; valWidth.textContent = '10 μm'; width = 10;
        sliderLambda.value = 1.55; valLambda.textContent = '1.55 μm'; lambda = 1.55;
        sliderNcore1.value = 1.500; valNcore1.textContent = '1.500'; nCore1 = 1.500;
        sliderNcore2.value = 1.500; valNcore2.textContent = '1.500'; nCore2 = 1.500;
        sliderNtop.value = 1.450; valNtop.textContent = '1.450'; nTop = 1.450;
        sliderNmid.value = 1.450; valNmid.textContent = '1.450'; nMid = 1.450;
        sliderNbot.value = 1.450; valNbot.textContent = '1.450'; nBot = 1.450;
        sliderSpeed.value = 1.0; valSpeed.textContent = '1.0x'; speed = 1.0;
        
        isPlaying = true;
        playIcon.style.display = 'none';
        pauseIcon.style.display = 'inline';
        
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
