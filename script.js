// Coupled Waveguide Simulator (Gaussian Pulse Propagation)
const canvas = document.getElementById('sim-canvas');
const ctx = canvas.getContext('2d');
const graphCanvas = document.getElementById('graph-canvas');
const gCtx = graphCanvas.getContext('2d');

// DOM Elements: Sliders & Controls
const sliderGap = document.getElementById('slider-gap');
const sliderWidth = document.getElementById('slider-width');
const sliderLambda = document.getElementById('slider-lambda');
const sliderNcore = document.getElementById('slider-ncore');
const sliderNclad = document.getElementById('slider-nclad');
const sliderSpeed = document.getElementById('slider-speed');

const valGap = document.getElementById('val-gap');
const valWidth = document.getElementById('val-width');
const valLambda = document.getElementById('val-lambda');
const valNcore = document.getElementById('val-ncore');
const valNclad = document.getElementById('val-nclad');
const valSpeed = document.getElementById('val-speed');

const btnPlayPause = document.getElementById('btn-play-pause');
const btnReset = document.getElementById('btn-reset');
const playIcon = btnPlayPause.querySelector('.play-icon');
const pauseIcon = btnPlayPause.querySelector('.pause-icon');

// Simulation Constants & Scale Factors
// Canvas dimensions: 800 x 360
// Horizontal: z-direction. 800px = 400 um (1 px = 0.5 um, pixelScaleZ = 2 px/um)
// Vertical: x-direction. 360px = 120 um (1 px = 0.33 um, pixelScaleY = 3 px/um)
const pixelScaleZ = 2.0; 
const pixelScaleY = 3.0;

let isPlaying = true;
let time = 0;
const launchMode = 'wg1'; // Fixed to Waveguide 1 launch as requested

// Physics parameters (recalculated when sliders change)
let gap = 15.0;       // um (separation d)
let width = 10.0;     // um (core width w)
let lambda = 1.55;    // um (wavelength)
let nCore = 1.50;     // Core refractive index
let nClad = 1.45;     // Cladding refractive index
let speed = 1.0;      // Time speed multiplier

// Pulse parameters
const sigmaPulse = 25.0; // Pulse width in um

// Solved constants
let k0, kDelta, V, kx, gamma, beta, kappa;
let phi1 = []; // Precomputed transverse profile for WG1
let phi2 = []; // Precomputed transverse profile for WG2

// Recalculate Physics & Modes
function updatePhysics() {
    // 1. Validation: Core index must be strictly greater than cladding index for guiding
    if (nClad >= nCore) {
        nClad = nCore - 0.005;
        sliderNclad.value = nClad;
        valNclad.textContent = nClad.toFixed(3);
    }
    
    k0 = 2 * Math.PI / lambda;
    kDelta = k0 * Math.sqrt(nCore * nCore - nClad * nClad);
    
    // V parameter
    V = kDelta * (width / 2);
    
    // Solve transcendental equation for single-slab mode
    let u = solveSlabMode(V);
    kx = (2 * u) / width;
    
    // Cladding decay constant (gamma)
    gamma = Math.sqrt(Math.max(0.001, kDelta * kDelta - kx * kx));
    
    // Propagation constant (beta)
    beta = Math.sqrt(Math.max(0.001, k0 * k0 * nCore * nCore - kx * kx));
    
    // Coupling coefficient (kappa)
    const denominator = beta * (width / 2 + 1 / gamma) * (kx * kx + gamma * gamma);
    kappa = (gamma * gamma * kx * kx * Math.exp(-gamma * gap)) / denominator;
    
    precomputeTransverseModes();
}

// Bisection solver for slab waveguide dispersion relation
function solveSlabMode(V) {
    let low = 0.001;
    let high = Math.min(V, Math.PI / 2 - 0.001);
    
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

// Precompute transverse mode profiles along Y axis
function precomputeTransverseModes() {
    const cy = canvas.height / 2; // Y-center = 180px
    const gap_px = gap * pixelScaleY;
    const w_px = width * pixelScaleY;
    
    // Core centers
    const y1 = cy - (gap_px / 2) - (w_px / 2);
    const y2 = cy + (gap_px / 2) + (w_px / 2);
    
    phi1 = new Array(canvas.height);
    phi2 = new Array(canvas.height);
    
    const norm = 1.0; 
    
    for (let y = 0; y < canvas.height; y++) {
        // Mode 1 (WG1)
        let dist1 = y - y1;
        let val1 = 0;
        if (Math.abs(dist1) <= w_px / 2) {
            val1 = norm * Math.cos(kx * (dist1 / pixelScaleY));
        } else {
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

// Render electric field intensity of propagating Gaussian Pulse
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
    
    // Step by 2 pixels for optimal rendering performance
    for (let z = 0; z < width_px; z += 2) {
        let z_um = z / pixelScaleZ;
        
        // Gaussian envelope for the pulse
        let distToCenter = z_um - pulseCenter;
        let gauss = Math.exp(- (distToCenter * distToCenter) / (2 * sigmaPulse * sigmaPulse));
        
        // Skip calculations if envelope is extremely small
        if (gauss < 0.001) {
            for (let y = 0; y < height_px; y += 2) {
                draw2x2Block(data, z, y, width_px, 8, 12, 22);
            }
            continue;
        }

        // Phase values for oscillation
        let argZ = beta * z_um - theta;
        let cos_bz = Math.cos(argZ);
        let sin_bz = Math.sin(argZ);
        
        // Evanescent coupling envelopes
        let cos_kz = Math.cos(kappa * z_um);
        let sin_kz = Math.sin(kappa * z_um);
        
        // Calculate coefficients for Waveguide 1 launch
        let c1 = gauss * cos_kz * cos_bz;
        let c2 = -gauss * sin_kz * sin_bz;
        
        for (let y = 0; y < height_px; y += 2) {
            let E = phi1[y] * c1 + phi2[y] * c2;
            
            // Map field intensity
            let r = 8, g = 12, b = 22; // Dark base color
            
            if (E > 0.015) {
                let amt = Math.min(1.0, E) * 235;
                r = Math.round(r + amt * 0.1);
                g = Math.round(g + amt * 0.85);
                b = Math.round(b + amt * 1.0);
            } else if (E < -0.015) {
                let amt = Math.min(1.0, -E) * 235;
                r = Math.round(r + amt * 1.0);
                g = Math.round(g + amt * 0.1);
                b = Math.round(b + amt * 0.7);
            }
            
            // Core highlight overlay
            let inCore = (Math.abs(y - y1) <= w_px / 2) || (Math.abs(y - y2) <= w_px / 2);
            if (inCore) {
                r = Math.min(255, r + 15);
                g = Math.min(255, g + 15);
                b = Math.min(255, b + 25);
            }
            
            draw2x2Block(data, z, y, width_px, r, g, b);
        }
    }
    
    ctx.putImageData(imgData, 0, 0);
    
    // Draw guide lines
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.18)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, y1 - w_px/2); ctx.lineTo(width_px, y1 - w_px/2);
    ctx.moveTo(0, y1 + w_px/2); ctx.lineTo(width_px, y1 + w_px/2);
    ctx.moveTo(0, y2 - w_px/2); ctx.lineTo(width_px, y2 - w_px/2);
    ctx.moveTo(0, y2 + w_px/2); ctx.lineTo(width_px, y2 + w_px/2);
    ctx.stroke();
    
    ctx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    ctx.font = 'bold 10px Inter';
    ctx.fillText('CORE 1 (입사 채널)', 15, y1 + 4);
    ctx.fillText('CORE 2 (결합 채널)', 15, y2 + 4);
}

// Helper to write a 2x2 pixel block into the ImageData array
function draw2x2Block(data, z, y, width, r, g, b) {
    const idx1 = (y * width + z) * 4;
    const idx2 = ((y + 1) * width + z) * 4;
    
    data[idx1] = r; data[idx1+1] = g; data[idx1+2] = b; data[idx1+3] = 255;
    data[idx1+4] = r; data[idx1+5] = g; data[idx1+6] = b; data[idx1+7] = 255;
    
    data[idx2] = r; data[idx2+1] = g; data[idx2+2] = b; data[idx2+3] = 255;
    data[idx2+4] = r; data[idx2+5] = g; data[idx2+6] = b; data[idx2+7] = 255;
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
    
    // 1. WG1 Power Curve
    gCtx.strokeStyle = '#38bdf8'; // Cyan
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let z_um = ((x - 40) / (w - 60)) * 400;
        let distToCenter = z_um - pulseCenter;
        let gauss = Math.exp(- (distToCenter * distToCenter) / (2 * sigmaPulse * sigmaPulse));
        let p1 = gauss * gauss * Math.pow(Math.cos(kappa * z_um), 2);
        
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
        let p2 = gauss * gauss * Math.pow(Math.sin(kappa * z_um), 2);
        
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
        // Advance time phase
        time += 0.05 * speed;
    }
    
    // Group velocity movement for Gaussian pulse center
    // Loops from -80 um to 480 um (canvas visible range is 0 to 400 um)
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
    
    sliderNcore.addEventListener('input', (e) => {
        nCore = parseFloat(e.target.value);
        valNcore.textContent = nCore.toFixed(3);
        updatePhysics();
    });

    sliderNclad.addEventListener('input', (e) => {
        nClad = parseFloat(e.target.value);
        valNclad.textContent = nClad.toFixed(3);
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
        sliderNcore.value = 1.50; valNcore.textContent = '1.500'; nCore = 1.50;
        sliderNclad.value = 1.45; valNclad.textContent = '1.450'; nClad = 1.45;
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
