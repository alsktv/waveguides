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
function drawElectricField(theta, probeX) {
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
    
    // 3. Draw Propagating 1D Transverse Electric Field Profile Curve E(y) at the moving probe position
    let probeX_px = probeX * pixelScaleX;
    
    // Draw the 1D transverse profile shape wiggling horizontally around the probe position
    ctx.strokeStyle = '#ffffff'; // Solid white wave line
    ctx.lineWidth = 2.5;
    ctx.shadowColor = '#06b6d4'; // Cyan neon glow
    ctx.shadowBlur = 8;
    ctx.beginPath();
    
    let cos_qx_c = Math.cos(q * probeX);
    let sin_qx_c = Math.sin(q * probeX);
    let cos_b1x_c = Math.cos(beta1 * probeX - theta);
    let sin_b1x_c = Math.sin(beta1 * probeX - theta);
    let sin_b2x_c = Math.sin(beta2 * probeX - theta);
    
    let c1_prof = cos_qx_c * cos_b1x_c + (delta / q) * sin_qx_c * sin_b1x_c;
    let c2_prof = (kappa / q) * sin_qx_c * sin_b2x_c;
    
    const profileScale = 50.0; // Profile displacement scaling factor
    
    for (let y = 0; y < height_px; y += 2) {
        let E_prof = phi1[y] * c1_prof + phi2[y] * c2_prof;
        let x_disp = probeX_px + E_prof * profileScale;
        
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

// Draw Coordinate Axes Ticks and Labels (linked to physical scale and user parameters)
function drawAxes() {
    // 1. Y Axis Ticks & Labels (Transverse axis)
    // Draw vertical axis line at x = 40px
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    ctx.lineWidth = 1.0;
    ctx.beginPath();
    ctx.moveTo(40, 10);
    ctx.lineTo(40, 350);
    ctx.stroke();
    
    // Draw tick marks and labels
    const yTicks = [20, 100, 180, 260, 340];
    ctx.fillStyle = 'rgba(255, 255, 255, 0.6)';
    ctx.font = '9px Fira Code, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    
    yTicks.forEach(y => {
        // Draw tick mark line (horizontal, 6px wide)
        ctx.beginPath();
        ctx.moveTo(34, y);
        ctx.lineTo(40, y);
        ctx.stroke();
        
        // Calculate physical Y coordinate in um: (y - 180) / pixelScaleY
        let y_um = (y - 180) / pixelScaleY;
        // Format with sign and 1 decimal place (except 0)
        let labelText = Math.abs(y_um) < 0.001 ? '0.0' : (y_um > 0 ? '+' : '') + y_um.toFixed(1);
        ctx.fillText(labelText + ' μm', 30, y);
    });
    
    // Label for Y Axis at the top
    ctx.textAlign = 'left';
    ctx.font = 'bold 10px Inter';
    ctx.fillStyle = '#ffffff';
    ctx.fillText('y (수직 횡방향)', 10, 15);
    
    // 2. X Axis Ticks & Labels (Propagation axis)
    // Draw horizontal axis line at y = 340px
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    ctx.beginPath();
    ctx.moveTo(40, 340);
    ctx.lineTo(780, 340);
    ctx.stroke();
    
    // Draw tick marks and labels
    const xTicks = [40, 220, 400, 580, 760];
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillStyle = 'rgba(255, 255, 255, 0.6)';
    ctx.font = '9px Fira Code, monospace';
    
    xTicks.forEach(x => {
        // Draw tick mark line (vertical, 6px tall)
        ctx.beginPath();
        ctx.moveTo(x, 340);
        ctx.lineTo(x, 346);
        ctx.stroke();
        
        // Calculate physical X coordinate in um: (x - 40) / pixelScaleX
        let x_um = (x - 40) / pixelScaleX;
        ctx.fillText(x_um.toFixed(0) + ' μm', x, 348);
    });
    
    // Label for X Axis at the right end
    ctx.textAlign = 'right';
    ctx.font = 'bold 10px Inter';
    ctx.fillStyle = '#ffffff';
    ctx.fillText('x (진행 방향) →', 780, 325);
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
function drawPowerGraph(probeX) {
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
    
    // --- Calculate dynamic marker points along the waveguide length ---
    let peakPoints = [];
    let returnPoints = [];
    let crossPoints = [];
    let maxPowerWG2 = (kappa * kappa) / (q * q);
    
    if (q > 0) {
        // A. Peak energy transfer points to WG2
        for (let n = 0; n < 8; n++) {
            let x_val = (2 * n + 1) * Math.PI / (2 * q);
            if (x_val > 400.0) break;
            peakPoints.push(x_val);
        }
        
        // B. Complete power return points to WG1 (P1 = 1.0)
        for (let n = 1; n < 8; n++) {
            let x_val = n * Math.PI / q;
            if (x_val > 400.0) break;
            returnPoints.push(x_val);
        }
        
        // C. 50:50 power crossover points (P1 = P2 = 0.5)
        if (kappa > 0 && maxPowerWG2 >= 0.5) {
            let ratio = q / (Math.sqrt(2) * kappa);
            if (ratio <= 1.0) {
                let base_asin = Math.asin(ratio);
                for (let n = 0; n < 8; n++) {
                    let x1 = (base_asin + n * Math.PI) / q;
                    if (x1 <= 400.0) crossPoints.push(x1);
                    
                    let x2 = (Math.PI - base_asin + n * Math.PI) / q;
                    if (x2 <= 400.0) crossPoints.push(x2);
                }
                crossPoints.sort((a, b) => a - b);
            }
        }
    }
    
    // --- Draw Crossover (50:50) Point Indicators ---
    gCtx.lineWidth = 1;
    if (crossPoints.length > 0) {
        crossPoints.slice(0, 3).forEach(x_val => {
            let x_px = 40 + (x_val / 400) * (w - 60);
            
            // Vertical dotted indicator
            gCtx.strokeStyle = 'rgba(168, 85, 247, 0.5)'; // Purple
            gCtx.setLineDash([2, 2]);
            gCtx.beginPath();
            gCtx.moveTo(x_px, h * 0.15);
            gCtx.lineTo(x_px, h * 0.85);
            gCtx.stroke();
            gCtx.setLineDash([]);
            
            // Intersection marker circle (at y corresponding to 0.5 power)
            gCtx.fillStyle = '#ffffff';
            gCtx.beginPath();
            gCtx.arc(x_px, h * 0.5, 3.5, 0, 2 * Math.PI);
            gCtx.fill();
            gCtx.strokeStyle = '#a855f7';
            gCtx.lineWidth = 1;
            gCtx.stroke();
            
            // Text Label
            gCtx.fillStyle = '#e9d5ff';
            gCtx.font = 'bold 8px Fira Code, monospace';
            gCtx.textAlign = 'center';
            gCtx.fillText(`${x_val.toFixed(1)}μm(50:50)`, x_px, h * 0.11);
        });
    }
    
    // --- Draw Peak Energy Transfer Points (P2 Max) ---
    peakPoints.slice(0, 2).forEach(x_val => {
        let x_px = 40 + (x_val / 400) * (w - 60);
        let y_px = h * 0.85 - maxPowerWG2 * (h * 0.7);
        
        // Vertical line
        gCtx.strokeStyle = 'rgba(244, 63, 94, 0.4)'; // Pink
        gCtx.setLineDash([2, 3]);
        gCtx.beginPath();
        gCtx.moveTo(x_px, h * 0.15);
        gCtx.lineTo(x_px, h * 0.85);
        gCtx.stroke();
        gCtx.setLineDash([]);
        
        // Marker circle at peak
        gCtx.fillStyle = '#f43f5e';
        gCtx.beginPath();
        gCtx.arc(x_px, y_px, 3.5, 0, 2 * Math.PI);
        gCtx.fill();
        
        // Text Label
        gCtx.fillStyle = '#fecdd3';
        gCtx.font = 'bold 8px Fira Code, monospace';
        gCtx.textAlign = 'center';
        let label = maxPowerWG2 > 0.999 ? 'P₂=1.0' : `Max P₂(${maxPowerWG2.toFixed(2)})`;
        gCtx.fillText(`${x_val.toFixed(1)}μm(${label})`, x_px, h * 0.9);
    });
    
    // --- Draw Return Points (P1 = 1.0) ---
    returnPoints.slice(0, 2).forEach(x_val => {
        let x_px = 40 + (x_val / 400) * (w - 60);
        
        // Vertical line
        gCtx.strokeStyle = 'rgba(56, 189, 248, 0.4)'; // Cyan
        gCtx.setLineDash([2, 3]);
        gCtx.beginPath();
        gCtx.moveTo(x_px, h * 0.15);
        gCtx.lineTo(x_px, h * 0.85);
        gCtx.stroke();
        gCtx.setLineDash([]);
        
        // Marker circle at P1 peak (y corresponding to 1.0)
        gCtx.fillStyle = '#38bdf8';
        gCtx.beginPath();
        gCtx.arc(x_px, h * 0.15, 3.5, 0, 2 * Math.PI);
        gCtx.fill();
        
        // Text Label
        gCtx.fillStyle = '#bae6fd';
        gCtx.font = 'bold 8px Fira Code, monospace';
        gCtx.textAlign = 'center';
        gCtx.fillText(`${x_val.toFixed(1)}μm(P₁=1.0)`, x_px, h * 0.06);
    });
    
    // 3. Draw Vertical Probe Position Indicator in the Power Graph
    let probeX_graph_px = 40 + (probeX / 400) * (w - 60);
    gCtx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
    gCtx.lineWidth = 1;
    gCtx.setLineDash([2, 2]);
    gCtx.beginPath();
    gCtx.moveTo(probeX_graph_px, h * 0.15);
    gCtx.lineTo(probeX_graph_px, h * 0.85);
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
    
    // Auto-propagate probeX from left (0 um) to right (400 um) over time
    let probeX = (time * 10) % 400;
    
    drawElectricField(time, probeX);
    drawPowerGraph(probeX);
    
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
