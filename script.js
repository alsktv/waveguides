// Coupled Waveguide Simulator
const canvas = document.getElementById('sim-canvas');
const ctx = canvas.getContext('2d');
const graphCanvas = document.getElementById('graph-canvas');
const gCtx = graphCanvas.getContext('2d');

// DOM Elements: Sliders & Controls
const sliderGap = document.getElementById('slider-gap');
const sliderWidth = document.getElementById('slider-width');
const sliderLambda = document.getElementById('slider-lambda');
const sliderDn = document.getElementById('slider-dn');
const sliderSpeed = document.getElementById('slider-speed');

const valGap = document.getElementById('val-gap');
const valWidth = document.getElementById('val-width');
const valLambda = document.getElementById('val-lambda');
const valDn = document.getElementById('val-dn');
const valSpeed = document.getElementById('val-speed');

const btnPlayPause = document.getElementById('btn-play-pause');
const btnReset = document.getElementById('btn-reset');
const playIcon = btnPlayPause.querySelector('.play-icon');
const pauseIcon = btnPlayPause.querySelector('.pause-icon');

const modeButtons = document.querySelectorAll('.mode-btn');

// Simulation Constants & Scale Factors
// Canvas dimensions: 800 x 360
// Horizontal: z-direction. 800px = 400 um (1 px = 0.5 um, pixelScaleZ = 2 px/um)
// Vertical: x-direction. 360px = 120 um (1 px = 0.33 um, pixelScaleY = 3 px/um)
const pixelScaleZ = 2.0; 
const pixelScaleY = 3.0;

let isPlaying = true;
let time = 0;
let launchMode = 'wg1'; // 'wg1', 'wg2', 'symmetric', 'antisymmetric'

// Physics parameters (recalculated when sliders change)
let gap = 15.0;       // um (separation d)
let width = 10.0;     // um (core width w)
let lambda = 1.55;    // um (wavelength)
let dn = 0.03;        // index difference (delta n)
let speed = 1.0;      // time evolution speed multiplier

// Solved constants
const nCore = 1.5;
let nClad, k0, kDelta, V, kx, gamma, beta, kappa;
let phi1 = []; // Precomputed transverse profile for WG1
let phi2 = []; // Precomputed transverse profile for WG2

// Recalculate Physics & Modes
function updatePhysics() {
    nClad = nCore - dn;
    k0 = 2 * Math.PI / lambda;
    kDelta = k0 * Math.sqrt(nCore * nCore - nClad * nClad);
    
    // V parameter
    V = kDelta * (width / 2);
    
    // Solve transcendental equation for single-slab mode: u * tan(u) = sqrt(V^2 - u^2)
    // where u = kx * w/2
    let u = solveSlabMode(V);
    kx = (2 * u) / width;
    
    // Cladding decay constant (gamma)
    gamma = Math.sqrt(Math.max(0.001, kDelta * kDelta - kx * kx));
    
    // Propagation constant (beta)
    beta = Math.sqrt(Math.max(0.001, k0 * k0 * nCore * nCore - kx * kx));
    
    // Coupling coefficient (kappa)
    // Formula for identical slab waveguides:
    // kappa = (gamma^2 * kx^2 * exp(-gamma * d)) / (beta * (w/2 + 1/gamma) * (kx^2 + gamma^2))
    const denominator = beta * (width / 2 + 1 / gamma) * (kx * kx + gamma * gamma);
    kappa = (gamma * gamma * kx * kx * Math.exp(-gamma * gap)) / denominator;
    
    precomputeTransverseModes();
}

// Bisection solver for u * tan(u) = sqrt(V^2 - u^2)
function solveSlabMode(V) {
    let low = 0.001;
    let high = Math.min(V, Math.PI / 2 - 0.001);
    
    // If V is too small to support guided modes, clamp to safe minimum
    if (V < 0.05) return V * 0.99;

    for (let i = 0; i < 15; i++) {
        let mid = (low + high) / 2;
        let lhs = mid * Math.tan(mid);
        let rhs = Math.sqrt(V * V - mid * mid);
        if (lhs > rhs) {
            high = mid;
        } else {
            low = mid;
        }
    }
    return (low + high) / 2;
}

// Precompute transverse mode profiles along the Y axis
function precomputeTransverseModes() {
    const cy = canvas.height / 2; // Y-center = 180px
    const gap_px = gap * pixelScaleY;
    const w_px = width * pixelScaleY;
    
    // Core centers
    const y1 = cy - (gap_px / 2) - (w_px / 2);
    const y2 = cy + (gap_px / 2) + (w_px / 2);
    
    phi1 = new Array(canvas.height);
    phi2 = new Array(canvas.height);
    
    // Find mode amplitude normalization factor A
    // (A * cos(kx * w/2) accounts for continuity at boundary)
    const norm = 1.0; 
    
    for (let y = 0; y < canvas.height; y++) {
        // Mode 1 (WG1)
        let dist1 = y - y1;
        let val1 = 0;
        if (Math.abs(dist1) <= w_px / 2) {
            // inside core
            val1 = norm * Math.cos(kx * (dist1 / pixelScaleY));
        } else {
            // outside core (evanescent tail)
            let boundary = (dist1 > 0 ? w_px / 2 : -w_px / 2);
            let tail = Math.abs(dist1) - w_px / 2;
            val1 = norm * Math.cos(kx * (boundary / pixelScaleY)) * Math.exp(-gamma * (tail / pixelScaleY));
        }
        phi1[y] = val1;
        
        // Mode 2 (WG2)
        let dist2 = y - y2;
        let val2 = 0;
        if (Math.abs(dist2) <= w_px / 2) {
            val2 = norm * Math.cos(kx * (dist2 / pixelScaleY));
        } else {
            let boundary = (dist2 > 0 ? w_px / 2 : -w_px / 2);
            let tail = Math.abs(dist2) - w_px / 2;
            val2 = norm * Math.cos(kx * (boundary / pixelScaleY)) * Math.exp(-gamma * (tail / pixelScaleY));
        }
        phi2[y] = val2;
    }
}

// Optimised Pixel Block Rendering
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
    
    // Draw columns (Z) and rows (Y) in steps of 2 for performance
    for (let z = 0; z < width_px; z += 2) {
        let z_um = z / pixelScaleZ;
        
        // Phase values for oscillation
        let argZ = beta * z_um - theta;
        let cos_bz = Math.cos(argZ);
        let sin_bz = Math.sin(argZ);
        
        // Evanescent coupling envelopes
        let cos_kz = Math.cos(kappa * z_um);
        let sin_kz = Math.sin(kappa * z_um);
        
        let c1 = 0, c2 = 0;
        
        // Envelope coefficients based on selected launch mode
        if (launchMode === 'wg1') {
            c1 = cos_kz * cos_bz;
            c2 = -sin_kz * sin_bz; 
        } else if (launchMode === 'wg2') {
            c1 = -sin_kz * sin_bz; 
            c2 = cos_kz * cos_bz;
        } else if (launchMode === 'symmetric') {
            let phaseS = (beta + kappa) * z_um - theta;
            let c_sym = Math.cos(phaseS) * Math.SQRT1_2;
            c1 = c_sym;
            c2 = c_sym;
        } else if (launchMode === 'antisymmetric') {
            let phaseA = (beta - kappa) * z_um - theta;
            let c_anti = Math.cos(phaseA) * Math.SQRT1_2;
            c1 = c_anti;
            c2 = -c_anti;
        }
        
        for (let y = 0; y < height_px; y += 2) {
            // Overlay field amplitude
            let E = phi1[y] * c1 + phi2[y] * c2;
            
            // Map positive field to blue/cyan, negative to pink/magenta
            let r = 8, g = 12, b = 22; // Base dark blue background
            
            if (E > 0.02) {
                let amt = Math.min(1.0, E) * 230;
                r = Math.round(r + amt * 0.1);
                g = Math.round(g + amt * 0.8);
                b = Math.round(b + amt * 1.0);
            } else if (E < -0.02) {
                let amt = Math.min(1.0, -E) * 230;
                r = Math.round(r + amt * 1.0);
                g = Math.round(g + amt * 0.1);
                b = Math.round(b + amt * 0.7);
            }
            
            // Add visual highlight overlay for the waveguide core boundaries
            let inCore = (Math.abs(y - y1) <= w_px / 2) || (Math.abs(y - y2) <= w_px / 2);
            if (inCore) {
                r = Math.min(255, r + 15);
                g = Math.min(255, g + 15);
                b = Math.min(255, b + 25);
            }
            
            // Fast ImageData block fill
            draw2x2Block(data, z, y, width_px, r, g, b);
        }
    }
    
    ctx.putImageData(imgData, 0, 0);
    
    // Draw waveguide border guide lines on top of canvas
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.18)';
    ctx.lineWidth = 1;
    
    // Waveguide 1 borders
    ctx.beginPath();
    ctx.moveTo(0, y1 - w_px/2); ctx.lineTo(width_px, y1 - w_px/2);
    ctx.moveTo(0, y1 + w_px/2); ctx.lineTo(width_px, y1 + w_px/2);
    // Waveguide 2 borders
    ctx.moveTo(0, y2 - w_px/2); ctx.lineTo(width_px, y2 - w_px/2);
    ctx.moveTo(0, y2 + w_px/2); ctx.lineTo(width_px, y2 + w_px/2);
    ctx.stroke();
    
    // Core label markings
    ctx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    ctx.font = 'bold 10px Inter';
    ctx.fillText('CORE 1 (WG1)', 15, y1 + 4);
    ctx.fillText('CORE 2 (WG2)', 15, y2 + 4);
}

// Inline helper for filling 2x2 blocks inside ImageData array
function draw2x2Block(data, z, y, width, r, g, b) {
    const idx1 = (y * width + z) * 4;
    const idx2 = ((y + 1) * width + z) * 4;
    
    data[idx1] = r; data[idx1+1] = g; data[idx1+2] = b; data[idx1+3] = 255;
    data[idx1+4] = r; data[idx1+5] = g; data[idx1+6] = b; data[idx1+7] = 255;
    
    data[idx2] = r; data[idx2+1] = g; data[idx2+2] = b; data[idx2+3] = 255;
    data[idx2+4] = r; data[idx2+5] = g; data[idx2+6] = b; data[idx2+7] = 255;
}

// Draw Longitudinal Power Density Graph
function drawPowerGraph() {
    const w = graphCanvas.width;
    const h = graphCanvas.height;
    
    gCtx.fillStyle = '#040711';
    gCtx.fillRect(0, 0, w, h);
    
    // Draw Grid Lines
    gCtx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
    gCtx.lineWidth = 1;
    // Y grid: 0.0, 0.5, 1.0
    for (let i = 0; i <= 2; i++) {
        let y_px = h * 0.15 + (i * h * 0.35);
        gCtx.beginPath();
        gCtx.moveTo(40, y_px);
        gCtx.lineTo(w - 20, y_px);
        gCtx.stroke();
    }
    
    // X grid: distance marks
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
    
    // Draw Y axis labels
    gCtx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    gCtx.font = '9px Fira Code';
    gCtx.fillText('1.0', 15, h * 0.18);
    gCtx.fillText('0.5', 15, h * 0.53);
    gCtx.fillText('0.0', 15, h * 0.88);
    
    // Plot lines
    gCtx.lineWidth = 2;
    
    // WG1 Power curve
    gCtx.strokeStyle = '#38bdf8'; // Cyan
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let z_um = ((x - 40) / (w - 60)) * 400;
        let p1 = 0;
        
        if (launchMode === 'wg1') {
            p1 = Math.pow(Math.cos(kappa * z_um), 2);
        } else if (launchMode === 'wg2') {
            p1 = Math.pow(Math.sin(kappa * z_um), 2);
        } else {
            p1 = 0.5; // Symmetric/Antisymmetric mode power remains 50:50
        }
        
        let y_px = h * 0.85 - p1 * (h * 0.7);
        if (x === 40) gCtx.moveTo(x, y_px);
        else gCtx.lineTo(x, y_px);
    }
    gCtx.stroke();
    
    // WG2 Power curve
    gCtx.strokeStyle = '#f43f5e'; // Pink
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let z_um = ((x - 40) / (w - 60)) * 400;
        let p2 = 0;
        
        if (launchMode === 'wg1') {
            p2 = Math.pow(Math.sin(kappa * z_um), 2);
        } else if (launchMode === 'wg2') {
            p2 = Math.pow(Math.cos(kappa * z_um), 2);
        } else {
            p2 = 0.5;
        }
        
        let y_px = h * 0.85 - p2 * (h * 0.7);
        if (x === 40) gCtx.moveTo(x, y_px);
        else gCtx.lineTo(x, y_px);
    }
    gCtx.stroke();
    
    // Legend labels
    gCtx.fillStyle = '#38bdf8';
    gCtx.fillText('Power 1 (WG1)', w - 110, h * 0.25);
    gCtx.fillStyle = '#f43f5e';
    gCtx.fillText('Power 2 (WG2)', w - 110, h * 0.38);
}

// Main Simulation Loop
function tick() {
    if (isPlaying) {
        // Advance time phase
        time += 0.05 * speed;
        drawElectricField(time);
    }
    requestAnimationFrame(tick);
}

// Slider event listeners & display updates
function setupSliders() {
    sliderGap.addEventListener('input', (e) => {
        gap = parseFloat(e.target.value);
        valGap.textContent = `${gap} μm`;
        updatePhysics();
        drawPowerGraph();
        if (!isPlaying) drawElectricField(time);
    });
    
    sliderWidth.addEventListener('input', (e) => {
        width = parseFloat(e.target.value);
        valWidth.textContent = `${width} μm`;
        updatePhysics();
        drawPowerGraph();
        if (!isPlaying) drawElectricField(time);
    });
    
    sliderLambda.addEventListener('input', (e) => {
        lambda = parseFloat(e.target.value);
        valLambda.textContent = `${lambda.toFixed(2)} μm`;
        updatePhysics();
        drawPowerGraph();
        if (!isPlaying) drawElectricField(time);
    });
    
    sliderDn.addEventListener('input', (e) => {
        dn = parseFloat(e.target.value);
        valDn.textContent = `${dn.toFixed(3)}`;
        updatePhysics();
        drawPowerGraph();
        if (!isPlaying) drawElectricField(time);
    });
    
    sliderSpeed.addEventListener('input', (e) => {
        speed = parseFloat(e.target.value);
        valSpeed.textContent = `${speed.toFixed(1)}x`;
    });
}

// Mode selection triggers
function setupControls() {
    modeButtons.forEach(btn => {
        btn.addEventListener('click', () => {
            modeButtons.forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            launchMode = btn.dataset.launch;
            drawPowerGraph();
            if (!isPlaying) drawElectricField(time);
        });
    });
    
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
        sliderDn.value = 0.03; valDn.textContent = '0.030'; dn = 0.03;
        sliderSpeed.value = 1.0; valSpeed.textContent = '1.0x'; speed = 1.0;
        
        modeButtons.forEach(b => b.classList.remove('active'));
        document.querySelector('[data-launch="wg1"]').classList.add('active');
        launchMode = 'wg1';
        
        isPlaying = true;
        playIcon.style.display = 'none';
        pauseIcon.style.display = 'inline';
        
        updatePhysics();
        drawPowerGraph();
    });
}

// Initialise Application
updatePhysics();
drawPowerGraph();
setupSliders();
setupControls();

// Start Sim loop
playIcon.style.display = 'none';
pauseIcon.style.display = 'inline';
tick();
