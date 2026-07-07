// Coupled Waveguide Simulator (Asymmetric 5-Layer Waveguide Solver & String Waveform Render)
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
const pixelScaleZ = 2.0; 
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

// Render electric field as horizontal wiggling strings showing cladding/evanescent leakage
function drawElectricField(theta, pulseCenter) {
    const width_px = canvas.width;
    const height_px = canvas.height;
    
    // Clear canvas
    ctx.fillStyle = '#040711';
    ctx.fillRect(0, 0, width_px, height_px);
    
    const w_px = width * pixelScaleY;
    const gap_px = gap * pixelScaleY;
    const cy = height_px / 2;
    const y1 = cy - (gap_px / 2) - (w_px / 2);
    const y2 = cy + (gap_px / 2) + (w_px / 2);
    
    // 1. Draw Waveguide Core Regions in the background
    // Core 1
    ctx.fillStyle = 'rgba(56, 189, 248, 0.04)';
    ctx.fillRect(0, y1 - w_px/2, width_px, w_px);
    // Core 2
    ctx.fillStyle = 'rgba(244, 63, 94, 0.04)';
    ctx.fillRect(0, y2 - w_px/2, width_px, w_px);
    
    // Draw Core borders
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.07)';
    ctx.lineWidth = 1;
    ctx.setLineDash([5, 5]);
    ctx.beginPath();
    ctx.moveTo(0, y1 - w_px/2); ctx.lineTo(width_px, y1 - w_px/2);
    ctx.moveTo(0, y1 + w_px/2); ctx.lineTo(width_px, y1 + w_px/2);
    ctx.moveTo(0, y2 - w_px/2); ctx.lineTo(width_px, y2 - w_px/2);
    ctx.moveTo(0, y2 + w_px/2); ctx.lineTo(width_px, y2 + w_px/2);
    ctx.stroke();
    ctx.setLineDash([]); // Reset dash
    
    // 2. Draw propagating transverse field curves (Y-slices)
    const stepY = 8; // Spacing between horizontal string lines
    const wiggleScale = 22.0; // Waveform wiggle amplitude
    
    const delta = window.g_delta;
    const q = window.g_q;
    const beta1 = window.g_beta1;
    const beta2 = window.g_beta2;
    
    for (let y = 12; y < height_px - 10; y += stepY) {
        // Color transition depending on layer position
        let color = '';
        let opacity = 0.25;
        
        if (y < y1 - w_px/2) {
            // Top cladding
            color = `rgba(139, 92, 246, ${opacity})`; // Purple
        } else if (Math.abs(y - y1) <= w_px/2) {
            // Core 1
            color = `rgba(56, 189, 248, ${opacity + 0.15})`; // Cyan
        } else if (y > y1 + w_px/2 && y < y2 - w_px/2) {
            // Middle cladding (evanescent gap)
            color = `rgba(168, 85, 247, ${opacity})`; // Violet
        } else if (Math.abs(y - y2) <= w_px/2) {
            // Core 2
            color = `rgba(244, 63, 94, ${opacity + 0.15})`; // Pink
        } else {
            // Bottom cladding
            color = `rgba(249, 115, 22, ${opacity})`; // Orange
        }
        
        ctx.beginPath();
        ctx.strokeStyle = color;
        ctx.lineWidth = (Math.abs(y - y1) <= w_px/2 || Math.abs(y - y2) <= w_px/2) ? 1.5 : 1.0;
        
        for (let z = 0; z < width_px; z += 4) {
            let z_um = z / pixelScaleZ;
            
            // Gaussian envelope
            let distToCenter = z_um - pulseCenter;
            let gauss = Math.exp(- (distToCenter * distToCenter) / (2 * sigmaPulse * sigmaPulse));
            
            let E = 0;
            if (gauss > 0.001) {
                // Wave propagation constants (asymmetric Coupled Mode equations)
                let cos_qz = Math.cos(q * z_um);
                let sin_qz = Math.sin(q * z_um);
                let cos_b1z = Math.cos(beta1 * z_um - theta);
                let sin_b1z = Math.sin(beta1 * z_um - theta);
                let sin_b2z = Math.sin(beta2 * z_um - theta);
                
                // WG1 input contribution
                let c1 = gauss * (cos_qz * cos_b1z + (delta / q) * sin_qz * sin_b1z);
                // WG2 coupling contribution (evanescent field transfer)
                let c2 = gauss * (kappa / q) * sin_qz * sin_b2z;
                
                E = phi1[y] * c1 + phi2[y] * c2;
            }
            
            // Offset coordinates horizontally/vertically (displacement waveform)
            let y_disp = y + E * wiggleScale;
            
            if (z === 0) ctx.moveTo(z, y_disp);
            else ctx.lineTo(z, y_disp);
        }
        ctx.stroke();
    }
    
    // Core label markings
    ctx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    ctx.font = 'bold 10px Inter';
    ctx.fillText(`CORE 1 (n = ${nCore1.toFixed(3)})`, 15, y1 + 4);
    ctx.fillText(`CORE 2 (n = ${nCore2.toFixed(3)})`, 15, y2 + 4);
    
    ctx.font = '9px Inter';
    ctx.fillText(`상부 클래딩 (n = ${nTop.toFixed(3)})`, 15, y1 - w_px/2 - 8);
    ctx.fillText(`중간 클래딩 (n = ${nMid.toFixed(3)})`, 15, cy + 4);
    ctx.fillText(`하부 클래딩 (n = ${nBot.toFixed(3)})`, 15, y2 + w_px/2 + 14);
}

// Draw Longitudinal Power Density Graph
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
    
    for (let z_val = 0; z_val <= 400; z_val += 100) {
        let x_px = 40 + (z_val / 400) * (w - 60);
        gCtx.beginPath();
        gCtx.moveTo(x_px, h * 0.15);
        gCtx.lineTo(x_px, h * 0.85);
        gCtx.stroke();
        
        gCtx.fillStyle = 'rgba(255,255,255,0.3)';
        gCtx.font = '9px Fira Code';
        gCtx.fillText(`${z_val}um`, x_px - 12, h * 0.95);
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
        let z_um = ((x - 40) / (w - 60)) * 400;
        let distToCenter = z_um - pulseCenter;
        let gauss = Math.exp(- (distToCenter * distToCenter) / (2 * sigmaPulse * sigmaPulse));
        let p1 = gauss * gauss * (Math.pow(Math.cos(q * z_um), 2) + Math.pow(delta / q, 2) * Math.pow(Math.sin(q * z_um), 2));
        
        let y_px = h * 0.85 - p1 * (h * 0.7);
        if (x === 40) gCtx.moveTo(x, y_px);
        else gCtx.lineTo(x, y_px);
    }
    gCtx.stroke();
    
    // 2. WG2 Power Curve
    gCtx.strokeStyle = '#f43f5e'; // Pink
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let z_um = ((x - 40) / (w - 60)) * 400;
        let distToCenter = z_um - pulseCenter;
        let gauss = Math.exp(- (distToCenter * distToCenter) / (2 * sigmaPulse * sigmaPulse));
        let p2 = gauss * gauss * Math.pow(kappa / q, 2) * Math.pow(Math.sin(q * z_um), 2);
        
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
    
    // Group velocity movement for Gaussian pulse center (loops along z)
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
