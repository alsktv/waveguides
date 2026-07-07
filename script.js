// Combined Waveguide Simulator (2D Coupled Propagation & 3D Slab/Rib Mode Solver using EIM)
const canvas = document.getElementById('sim-canvas');
const ctx = canvas.getContext('2d');
const graphCanvas = document.getElementById('graph-canvas');
const gCtx = graphCanvas.getContext('2d');

// Tabs management
let activeTab = 'coupled';

function switchTab(tabId) {
    activeTab = tabId;
    
    // Toggle active classes
    document.getElementById('btn-tab-coupled').classList.toggle('active', tabId === 'coupled');
    document.getElementById('btn-tab-3d').classList.toggle('active', tabId === '3d');
    
    document.getElementById('tab-coupled').classList.toggle('active-tab', tabId === 'coupled');
    document.getElementById('tab-3d').classList.toggle('active-tab', tabId === '3d');
    
    if (tabId === '3d') {
        updatePhysics3D();
    } else {
        updatePhysics();
    }
}

document.getElementById('btn-tab-coupled').addEventListener('click', () => switchTab('coupled'));
document.getElementById('btn-tab-3d').addEventListener('click', () => switchTab('3d'));

// ----------------------------------------------------
// PART 1: 2D Planar Coupled Waveguide Simulator
// ----------------------------------------------------
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

const pixelScaleX = 2.0; 
let pixelScaleY = 3.0;

let isPlaying = true;
let time = 0;

let gap = 15.0;       
let width = 10.0;     
let lambda = 1.55;    
let nCore1 = 1.500;   
let nCore2 = 1.500;   
let nTop = 1.450;     
let nMid = 1.450;     
let nBot = 1.450;     
let speed = 1.0;      

let probeX = 200.0;

let k0, kx1, kx2, gTop, gMid1, gMid2, gBot, phi0, psi0, kappa;
let phi1 = []; 
let phi2 = []; 

function getGuidedModeCount(nCore, nCladL, nCladR, width_um, lambda_um) {
    const k0_val = 2 * Math.PI / lambda_um;
    const cladMax = Math.max(nCladL, nCladR);
    const cladMin = Math.min(nCladL, nCladR);
    if (nCore <= cladMax) return 0;
    
    const kx_cutoff = k0_val * Math.sqrt(nCore * nCore - cladMax * cladMax);
    const g_other_cutoff = k0_val * Math.sqrt(Math.max(0, cladMax * cladMax - cladMin * cladMin));
    const phi_cutoff = kx_cutoff * width_um - Math.atan(g_other_cutoff / kx_cutoff);
    const modeCount = 1 + Math.floor(phi_cutoff / Math.PI);
    return Math.max(1, modeCount);
}

function updatePhysics() {
    pixelScaleY = Math.min(28.0, Math.max(2.0, 280 / (width * 2 + gap)));

    let cladMax1 = Math.max(nTop, nMid);
    if (nCore1 <= cladMax1) {
        nCore1 = cladMax1 + 0.005;
        sliderNcore1.value = nCore1;
        document.getElementById('input-ncore1').value = nCore1.toFixed(3);
    }
    
    let cladMax2 = Math.max(nMid, nBot);
    if (nCore2 <= cladMax2) {
        nCore2 = cladMax2 + 0.005;
        sliderNcore2.value = nCore2;
        document.getElementById('input-ncore2').value = nCore2.toFixed(3);
    }
    
    k0 = 2 * Math.PI / lambda;
    
    let betaMin1 = k0 * cladMax1 + 0.001;
    let betaMax1 = k0 * nCore1 - 0.001;
    let beta1 = solveBeta(nCore1, nTop, nMid, betaMin1, betaMax1);
    
    kx1 = Math.sqrt(k0*k0*nCore1*nCore1 - beta1*beta1);
    gTop = Math.sqrt(beta1*beta1 - k0*k0*nTop*nTop);
    gMid1 = Math.sqrt(beta1*beta1 - k0*k0*nMid*nMid);
    phi0 = 0.5 * (Math.atan(gTop / kx1) - Math.atan(gMid1 / kx1));
    
    let betaMin2 = k0 * cladMax2 + 0.001;
    let betaMax2 = k0 * nCore2 - 0.001;
    let beta2 = solveBeta(nCore2, nMid, nBot, betaMin2, betaMax2);
    
    kx2 = Math.sqrt(k0*k0*nCore2*nCore2 - beta2*beta2);
    gMid2 = Math.sqrt(beta2*beta2 - k0*k0*nMid*nMid);
    gBot = Math.sqrt(beta2*beta2 - k0*k0*nBot*nBot);
    psi0 = 0.5 * (Math.atan(gMid2 / kx2) - Math.atan(gBot / kx2));
    
    let betaAvg = (beta1 + beta2) / 2;
    let delta = (beta1 - beta2) / 2;
    
    let gMidAvg = (gMid1 + gMid2) / 2;
    const denominator = Math.sqrt(beta1 * beta2) * (width / 2 + 1 / gTop + 1 / gMid1) * (width / 2 + 1 / gBot + 1 / gMid2);
    kappa = (2 * kx1 * kx2 * gMidAvg * Math.exp(-gMidAvg * gap)) / denominator;
    
    window.g_beta1 = beta1;
    window.g_beta2 = beta2;
    window.g_delta = delta;
    window.g_q = Math.sqrt(kappa * kappa + delta * delta);
    
    let Ic1 = width/2 + Math.sin(kx1 * width) * Math.cos(2 * phi0) / (2 * kx1);
    let It1 = Math.pow(Math.cos(kx1 * width/2 + phi0), 2) / (2 * gTop);
    let Ib1 = Math.pow(Math.cos(kx1 * width/2 - phi0), 2) / (2 * gMid1);
    let norm1 = 1 / Math.sqrt(Ic1 + It1 + Ib1);

    let Ic2 = width/2 + Math.sin(kx2 * width) * Math.cos(2 * psi0) / (2 * kx2);
    let It2 = Math.pow(Math.cos(kx2 * width/2 + psi0), 2) / (2 * gMid2);
    let Ib2 = Math.pow(Math.cos(kx2 * width/2 - psi0), 2) / (2 * gBot);
    let norm2 = 1 / Math.sqrt(Ic2 + It2 + Ib2);

    let m1 = getGuidedModeCount(nCore1, nTop, nMid, width, lambda);
    let m2 = getGuidedModeCount(nCore2, nMid, nBot, width, lambda);
    document.getElementById('mode-count-1').textContent = `${m1}개`;
    document.getElementById('mode-count-2').textContent = `${m2}개`;

    precomputeTransverseModes(norm1, norm2);
}

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

function precomputeTransverseModes(norm1, norm2) {
    const cy = canvas.height / 2;
    const gap_px = gap * pixelScaleY;
    const w_px = width * pixelScaleY;
    
    const y1 = cy - (gap_px / 2) - (w_px / 2);
    const y2 = cy + (gap_px / 2) + (w_px / 2);
    
    phi1 = new Array(canvas.height);
    phi2 = new Array(canvas.height);
    
    for (let y = 0; y < canvas.height; y++) {
        let dist = y - y1;
        let y_um = dist / pixelScaleY;
        if (Math.abs(dist) <= w_px / 2) {
            phi1[y] = norm1 * Math.cos(kx1 * y_um - phi0);
        } else if (dist < -w_px / 2) {
            let boundary_um = -width / 2;
            let tail_um = -y_um - width / 2;
            phi1[y] = norm1 * Math.cos(kx1 * boundary_um - phi0) * Math.exp(-gTop * tail_um);
        } else {
            let boundary_um = width / 2;
            let tail_um = y_um - width / 2;
            phi1[y] = norm1 * Math.cos(kx1 * boundary_um - phi0) * Math.exp(-gMid1 * tail_um);
        }
    }
    
    for (let y = 0; y < canvas.height; y++) {
        let dist = y - y2;
        let y_um = dist / pixelScaleY;
        if (Math.abs(dist) <= w_px / 2) {
            phi2[y] = norm2 * Math.cos(kx2 * y_um - psi0);
        } else if (dist < -w_px / 2) {
            let boundary_um = -width / 2;
            let tail_um = -y_um - width / 2;
            phi2[y] = norm2 * Math.cos(kx2 * boundary_um - psi0) * Math.exp(-gMid2 * tail_um);
        } else {
            let boundary_um = width / 2;
            let tail_um = y_um - width / 2;
            phi2[y] = norm2 * Math.cos(kx2 * boundary_um - psi0) * Math.exp(-gBot * tail_um);
        }
    }
}

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
    
    for (let x = 0; x < width_px; x += 2) {
        let x_um = x / pixelScaleX;
        
        let cos_qx = Math.cos(q * x_um);
        let sin_qx = Math.sin(q * x_um);
        let cos_b1x = Math.cos(beta1 * x_um - theta);
        let sin_b1x = Math.sin(beta1 * x_um - theta);
        let sin_b2x = Math.sin(beta2 * x_um - theta);
        
        let c1 = cos_qx * cos_b1x + (delta / q) * sin_qx * sin_b1x;
        let c2 = (kappa / q) * sin_qx * sin_b2x;
        
        for (let y = 0; y < height_px; y += 2) {
            let E = phi1[y] * c1 + phi2[y] * c2;
            
            let r = 8, g = 12, b = 22; 
            
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
    
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(0, y1 - w_px/2); ctx.lineTo(width_px, y1 - w_px/2);
    ctx.moveTo(0, y1 + w_px/2); ctx.lineTo(width_px, y1 + w_px/2);
    ctx.moveTo(0, y2 - w_px/2); ctx.lineTo(width_px, y2 - w_px/2);
    ctx.moveTo(0, y2 + w_px/2); ctx.lineTo(width_px, y2 + w_px/2);
    ctx.stroke();
    
    let probeX_px = probeX * pixelScaleX;
    
    ctx.strokeStyle = '#ffffff'; 
    ctx.lineWidth = 2.5;
    ctx.shadowColor = '#06b6d4'; 
    ctx.shadowBlur = 8;
    ctx.beginPath();
    
    let cos_qx_c = Math.cos(q * probeX);
    let sin_qx_c = Math.sin(q * probeX);
    let cos_b1x_c = Math.cos(beta1 * probeX - theta);
    let sin_b1x_c = Math.sin(beta1 * probeX - theta);
    let sin_b2x_c = Math.sin(beta2 * probeX - theta);
    
    let c1_prof = cos_qx_c * cos_b1x_c + (delta / q) * sin_qx_c * sin_b1x_c;
    let c2_prof = (kappa / q) * sin_qx_c * sin_b2x_c;
    
    const profileScale = 50.0; 
    
    for (let y = 0; y < height_px; y += 2) {
        let E_prof = phi1[y] * c1_prof + phi2[y] * c2_prof;
        let x_disp = probeX_px + E_prof * profileScale;
        
        if (y === 0) ctx.moveTo(x_disp, y);
        else ctx.lineTo(x_disp, y);
    }
    ctx.stroke();
    ctx.shadowBlur = 0; 
    
    ctx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    ctx.font = 'bold 10px Inter';
    ctx.fillText(`CORE 1 (n = ${nCore1.toFixed(3)})`, 15, y1 + 4);
    ctx.fillText(`CORE 2 (n = ${nCore2.toFixed(3)})`, 15, y2 + 4);
    
    ctx.font = '9px Inter';
    ctx.fillText(`상부 클래딩 (n = ${nTop.toFixed(3)})`, 15, y1 - w_px/2 - 8);
    ctx.fillText(`중간 클래딩 (n = ${nMid.toFixed(3)})`, 15, cy + 4);
    ctx.fillText(`하부 클래딩 (n = ${nBot.toFixed(3)})`, 15, y2 + w_px/2 + 14);
    
    drawAxes();
}

function drawAxes() {
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    ctx.lineWidth = 1.0;
    ctx.beginPath();
    ctx.moveTo(40, 10);
    ctx.lineTo(40, 350);
    ctx.stroke();
    
    const yTicks = [20, 100, 180, 260, 340];
    ctx.fillStyle = 'rgba(255, 255, 255, 0.6)';
    ctx.font = '9px Fira Code, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    
    yTicks.forEach(y => {
        ctx.beginPath();
        ctx.moveTo(34, y);
        ctx.lineTo(40, y);
        ctx.stroke();
        
        let y_um = (y - 180) / pixelScaleY;
        let labelText = Math.abs(y_um) < 0.001 ? '0.0' : (y_um > 0 ? '+' : '') + y_um.toFixed(1);
        ctx.fillText(labelText + ' μm', 30, y);
    });
    
    ctx.textAlign = 'left';
    ctx.font = 'bold 10px Inter';
    ctx.fillStyle = '#ffffff';
    ctx.fillText('y (수직 횡방향)', 10, 15);
    
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    ctx.beginPath();
    ctx.moveTo(40, 340);
    ctx.lineTo(780, 340);
    ctx.stroke();
    
    const xTicks = [40, 220, 400, 580, 760];
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillStyle = 'rgba(255, 255, 255, 0.6)';
    ctx.font = '9px Fira Code, monospace';
    
    xTicks.forEach(x => {
        ctx.beginPath();
        ctx.moveTo(x, 340);
        ctx.lineTo(x, 346);
        ctx.stroke();
        
        let x_um = (x - 40) / pixelScaleX;
        ctx.fillText(x_um.toFixed(0) + ' μm', x, 348);
    });
    
    ctx.textAlign = 'right';
    ctx.font = 'bold 10px Inter';
    ctx.fillStyle = '#ffffff';
    ctx.fillText('x (진행 방향) →', 780, 325);
}

function draw2x2Block(data, x, y, width, r, g, b) {
    const idx1 = (y * width + x) * 4;
    const idx2 = ((y + 1) * width + x) * 4;
    data[idx1] = r; data[idx1+1] = g; data[idx1+2] = b; data[idx1+3] = 255;
    data[idx1+4] = r; data[idx1+5] = g; data[idx1+6] = b; data[idx1+7] = 255;
    data[idx2] = r; data[idx2+1] = g; data[idx2+2] = b; data[idx2+3] = 255;
    data[idx2+4] = r; data[idx2+5] = g; data[idx2+6] = b; data[idx2+7] = 255;
}

function drawPowerGraph(probeX) {
    const w = graphCanvas.width;
    const h = graphCanvas.height;
    
    gCtx.fillStyle = '#040711';
    gCtx.fillRect(0, 0, w, h);
    
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
    
    gCtx.lineWidth = 2.5;
    const delta = window.g_delta;
    const q = window.g_q;
    
    gCtx.strokeStyle = '#38bdf8'; 
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let x_um = ((x - 40) / (w - 60)) * 400;
        let p1 = Math.pow(Math.cos(q * x_um), 2) + Math.pow(delta / q, 2) * Math.pow(Math.sin(q * x_um), 2);
        let y_px = h * 0.85 - p1 * (h * 0.7);
        if (x === 40) gCtx.moveTo(x, y_px);
        else gCtx.lineTo(x, y_px);
    }
    gCtx.stroke();
    
    gCtx.strokeStyle = '#f43f5e'; 
    gCtx.beginPath();
    for (let x = 40; x < w - 20; x++) {
        let x_um = ((x - 40) / (w - 60)) * 400;
        let p2 = Math.pow(kappa / q, 2) * Math.pow(Math.sin(q * x_um), 2);
        let y_px = h * 0.85 - p2 * (h * 0.7);
        if (x === 40) gCtx.moveTo(x, y_px);
        else gCtx.lineTo(x, y_px);
    }
    gCtx.stroke();
    
    let peakPoints = [];
    let returnPoints = [];
    let crossPoints = [];
    let maxPowerWG2 = (kappa * kappa) / (q * q);
    
    if (q > 0) {
        for (let n = 0; n < 8; n++) {
            let x_val = (2 * n + 1) * Math.PI / (2 * q);
            if (x_val > 400.0) break;
            peakPoints.push(x_val);
        }
        for (let n = 1; n < 8; n++) {
            let x_val = n * Math.PI / q;
            if (x_val > 400.0) break;
            returnPoints.push(x_val);
        }
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
    
    gCtx.lineWidth = 1;
    if (crossPoints.length > 0) {
        crossPoints.slice(0, 3).forEach(x_val => {
            let x_px = 40 + (x_val / 400) * (w - 60);
            gCtx.strokeStyle = 'rgba(168, 85, 247, 0.5)';
            gCtx.setLineDash([2, 2]);
            gCtx.beginPath();
            gCtx.moveTo(x_px, h * 0.15);
            gCtx.lineTo(x_px, h * 0.85);
            gCtx.stroke();
            gCtx.setLineDash([]);
            
            gCtx.fillStyle = '#ffffff';
            gCtx.beginPath();
            gCtx.arc(x_px, h * 0.5, 3.5, 0, 2 * Math.PI);
            gCtx.fill();
            gCtx.strokeStyle = '#a855f7';
            gCtx.stroke();
            
            gCtx.fillStyle = '#e9d5ff';
            gCtx.font = 'bold 8px Fira Code, monospace';
            gCtx.textAlign = 'center';
            gCtx.fillText(`${x_val.toFixed(1)}μm(50:50)`, x_px, h * 0.11);
        });
    }
    
    peakPoints.slice(0, 2).forEach(x_val => {
        let x_px = 40 + (x_val / 400) * (w - 60);
        let y_px = h * 0.85 - maxPowerWG2 * (h * 0.7);
        gCtx.strokeStyle = 'rgba(244, 63, 94, 0.4)';
        gCtx.setLineDash([2, 3]);
        gCtx.beginPath();
        gCtx.moveTo(x_px, h * 0.15);
        gCtx.lineTo(x_px, h * 0.85);
        gCtx.stroke();
        gCtx.setLineDash([]);
        
        gCtx.fillStyle = '#f43f5e';
        gCtx.beginPath();
        gCtx.arc(x_px, y_px, 3.5, 0, 2 * Math.PI);
        gCtx.fill();
        
        gCtx.fillStyle = '#fecdd3';
        gCtx.font = 'bold 8px Fira Code, monospace';
        gCtx.textAlign = 'center';
        let label = maxPowerWG2 > 0.999 ? 'P₂=1.0' : `Max P₂(${maxPowerWG2.toFixed(2)})`;
        gCtx.fillText(`${x_val.toFixed(1)}μm(${label})`, x_px, h * 0.9);
    });
    
    returnPoints.slice(0, 2).forEach(x_val => {
        let x_px = 40 + (x_val / 400) * (w - 60);
        gCtx.strokeStyle = 'rgba(56, 189, 248, 0.4)';
        gCtx.setLineDash([2, 3]);
        gCtx.beginPath();
        gCtx.moveTo(x_px, h * 0.15);
        gCtx.lineTo(x_px, h * 0.85);
        gCtx.stroke();
        gCtx.setLineDash([]);
        
        gCtx.fillStyle = '#38bdf8';
        gCtx.beginPath();
        gCtx.arc(x_px, h * 0.15, 3.5, 0, 2 * Math.PI);
        gCtx.fill();
        
        gCtx.fillStyle = '#bae6fd';
        gCtx.font = 'bold 8px Fira Code, monospace';
        gCtx.textAlign = 'center';
        gCtx.fillText(`${x_val.toFixed(1)}μm(P₁=1.0)`, x_px, h * 0.06);
    });
    
    let probeX_graph_px = 40 + (probeX / 400) * (w - 60);
    gCtx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
    gCtx.lineWidth = 1;
    gCtx.setLineDash([2, 2]);
    gCtx.beginPath();
    gCtx.moveTo(probeX_graph_px, h * 0.15);
    gCtx.lineTo(probeX_graph_px, h * 0.85);
    gCtx.stroke();
    gCtx.setLineDash([]);
    
    if (!isPlaying) {
        let p1_val = Math.pow(Math.cos(q * probeX), 2) + Math.pow(delta / q, 2) * Math.pow(Math.sin(q * probeX), 2);
        let p2_val = Math.pow(kappa / q, 2) * Math.pow(Math.sin(q * probeX), 2);
        let y_p1_px = h * 0.85 - p1_val * (h * 0.7);
        let y_p2_px = h * 0.85 - p2_val * (h * 0.7);
        
        gCtx.fillStyle = '#38bdf8';
        gCtx.beginPath();
        gCtx.arc(probeX_graph_px, y_p1_px, 5.5, 0, 2 * Math.PI);
        gCtx.fill();
        gCtx.fillStyle = '#ffffff';
        gCtx.beginPath();
        gCtx.arc(probeX_graph_px, y_p1_px, 2.0, 0, 2 * Math.PI);
        gCtx.fill();
        
        gCtx.fillStyle = '#f43f5e';
        gCtx.beginPath();
        gCtx.arc(probeX_graph_px, y_p2_px, 5.5, 0, 2 * Math.PI);
        gCtx.fill();
        gCtx.fillStyle = '#ffffff';
        gCtx.beginPath();
        gCtx.arc(probeX_graph_px, y_p2_px, 2.0, 0, 2 * Math.PI);
        gCtx.fill();
        
        gCtx.font = 'bold 9px Fira Code, monospace';
        gCtx.textAlign = 'left';
        let offsetP1 = (y_p1_px < y_p2_px) ? -6 : 12;
        let offsetP2 = (y_p2_px < y_p1_px) ? -6 : 12;
        gCtx.fillStyle = '#bae6fd';
        gCtx.fillText(`P₁:${p1_val.toFixed(3)}`, probeX_graph_px + 8, y_p1_px + offsetP1);
        gCtx.fillStyle = '#fecdd3';
        gCtx.fillText(`P₂:${p2_val.toFixed(3)}`, probeX_graph_px + 8, y_p2_px + offsetP2);
    }
    
    gCtx.font = 'bold 9px Inter';
    gCtx.fillStyle = '#38bdf8';
    gCtx.fillText('도파로 1 전력 (WG1 Power)', w - 160, h * 0.25);
    gCtx.fillStyle = '#f43f5e';
    gCtx.fillText('도파로 2 전력 (WG2 Power)', w - 160, h * 0.38);
}

function bindParam(sliderId, inputId, applyBtnId, minVal, maxVal, isInt, updateVarFn) {
    const slider = document.getElementById(sliderId);
    const input = document.getElementById(inputId);
    const btn = document.getElementById(applyBtnId);
    
    slider.addEventListener('input', (e) => {
        let val = parseFloat(e.target.value);
        input.value = isInt ? val.toFixed(0) : val.toFixed(3);
        updateVarFn(val);
        updatePhysics();
    });
    
    function applyManualValue() {
        let val = parseFloat(input.value);
        if (isNaN(val)) return;
        
        val = Math.max(minVal, Math.min(maxVal, val));
        input.value = isInt ? val.toFixed(0) : val.toFixed(3);
        slider.value = Math.max(parseFloat(slider.min), Math.min(parseFloat(slider.max), val));
        
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

function setupSliders() {
    bindParam('slider-gap', 'input-gap', 'btn-apply-gap', 1.0, 100.0, false, (v) => { gap = v; });
    bindParam('slider-width', 'input-width', 'btn-apply-width', 0.1, 100.0, false, (v) => { width = v; });
    bindParam('slider-lambda', 'input-lambda', 'btn-apply-lambda', 0.1, 10.0, false, (v) => { lambda = v; });
    bindParam('slider-ncore1', 'input-ncore1', 'btn-apply-ncore1', 1.0, 4.0, false, (v) => { nCore1 = v; });
    bindParam('slider-ncore2', 'input-ncore2', 'btn-apply-ncore2', 1.0, 4.0, false, (v) => { nCore2 = v; });
    bindParam('slider-ntop', 'input-ntop', 'btn-apply-ntop', 1.0, 4.0, false, (v) => { nTop = v; });
    bindParam('slider-nmid', 'input-nmid', 'btn-apply-nmid', 1.0, 4.0, false, (v) => { nMid = v; });
    bindParam('slider-nbot', 'input-nbot', 'btn-apply-nbot', 1.0, 4.0, false, (v) => { nBot = v; });
    
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

// ----------------------------------------------------
// PART 2: 3D Slab/Rib EIM Waveguide Solver
// ----------------------------------------------------
const canvas3D = document.getElementById('canvas-3d');
const ctx3D = canvas3D.getContext('2d');
const cutlinesCanvas = document.getElementById('canvas-3d-cutlines');
const cutCtx = cutlinesCanvas.getContext('2d');

let wgType3D = 'rib'; 
let width3D = 3.0;     
let height3D = 2.0;    
let slab3D = 0.8;      
let lambda3D = 1.55;   
let nCore3D = 1.500;   
let nCover3D = 1.000;  
let nSub3D = 1.450;    

let nEff3D = 0.0;
let cutoff3D = false;

// 2D grids precalculated profiles
let G_I = [];  
let G_II = []; 
let F_y = [];  
let kz1_3d, kz2_3d, ky_3d;
let gCover1_3d, gSub1_3d, gCover2_3d, gSub2_3d, gY_3d;
let phi01_3d, phi02_3d;

// Bisection slab solver helper (with mode order m support)
function solveSlab3D(nCore, nCladL, nCladR, thickness, wavelength, m = 0) {
    let k0_val = 2 * Math.PI / wavelength;
    let cladMax = Math.max(nCladL, nCladR);
    if (nCore <= cladMax) return null;
    
    let low = k0_val * cladMax + 0.0001;
    let high = k0_val * nCore - 0.0001;
    let beta = null;
    
    for (let i = 0; i < 24; i++) {
        let mid = (low + high) / 2;
        let kx = Math.sqrt(k0_val*k0_val*nCore*nCore - mid*mid);
        let gL = Math.sqrt(mid*mid - k0_val*k0_val*nCladL*nCladL);
        let gR = Math.sqrt(mid*mid - k0_val*k0_val*nCladR*nCladR);
        
        let val = kx * thickness - Math.atan(gL / kx) - Math.atan(gR / kx) - m * Math.PI;
        if (val > 0) {
            low = mid;
            beta = mid;
        } else {
            high = mid;
        }
    }
    return beta;
}

// Jet-based thermal color map function (optimized dark navy background for v=0)
function getJetColor(v) {
    v = Math.max(0, Math.min(1, v));
    let r = 0, g = 0, b = 0;
    
    if (v < 0.25) {
        // Dark Blue (10, 15, 60) to Cyan (0, 255, 255)
        let t = v / 0.25;
        r = Math.round(10 * (1 - t) + 0 * t);
        g = Math.round(15 * (1 - t) + 255 * t);
        b = Math.round(60 * (1 - t) + 255 * t);
    } else if (v < 0.5) {
        // Cyan (0, 255, 255) to Green (0, 255, 0)
        let t = (v - 0.25) / 0.25;
        r = 0;
        g = 255;
        b = Math.round(255 * (1 - t));
    } else if (v < 0.75) {
        // Green (0, 255, 0) to Yellow (255, 255, 0)
        let t = (v - 0.5) / 0.25;
        r = Math.round(255 * t);
        g = 255;
        b = 0;
    } else {
        // Yellow (255, 255, 0) to Dark Red (180, 0, 0)
        let t = (v - 0.75) / 0.25;
        r = Math.round(255 * (1 - t) + 180 * t);
        g = Math.round(255 * (1 - t) + 0 * t);
        b = 0;
    }
    return { r, g, b };
}

function updatePhysics3D() {
    // Validations
    let cladMax = Math.max(nCover3D, nSub3D);
    if (nCore3D <= cladMax) {
        nCore3D = cladMax + 0.005;
        document.getElementById('slider-3d-ncore').value = nCore3D;
        document.getElementById('input-3d-ncore').value = nCore3D.toFixed(3);
    }
    
    if (slab3D > height3D) {
        slab3D = height3D;
        document.getElementById('slider-3d-slab').value = slab3D;
        document.getElementById('input-3d-slab').value = slab3D.toFixed(1);
    }
    
    const k0_3d = 2 * Math.PI / lambda3D;
    
    // Read the selected mode index (0 = fundamental, 1 = first-order)
    let selectedMode3D = parseInt(document.getElementById('select-3d-mode').value);
    
    // Step 1: Vertical Slab I (Center thickness H, always fundamental m_z = 0)
    let beta_I = solveSlab3D(nCore3D, nCover3D, nSub3D, height3D, lambda3D, 0);
    let nEff_I = beta_I ? beta_I / k0_3d : nSub3D;
    
    // Step 2: Vertical Slab II (Slab wing height h, always fundamental m_z = 0)
    let nEff_II = nSub3D;
    let beta_II = null;
    
    if (wgType3D === 'rib' && slab3D > 0.02) {
        beta_II = solveSlab3D(nCore3D, nCover3D, nSub3D, slab3D, lambda3D, 0);
        if (beta_II) {
            nEff_II = beta_II / k0_3d;
        }
    }
    
    // Step 3: Horizontal Slab (Width W, core nEff_I, cladding nEff_II, mode m_y = selectedMode3D)
    let beta_final = null;
    if (beta_I && nEff_I > nEff_II) {
        beta_final = solveSlab3D(nEff_I, nEff_II, nEff_II, width3D, lambda3D, selectedMode3D);
    }
    
    // Update badge description text
    let geomText = wgType3D === 'rib' ? 'Rib Waveguide' : 'Strip Waveguide';
    let modeText = selectedMode3D === 0 ? 'm = 0 (Fundamental Mode)' : 'm = 1 (1st-order Mode)';
    document.getElementById('badge-mode-type').textContent = `${geomText} (${modeText})`;
    
    if (beta_final) {
        cutoff3D = false;
        nEff3D = beta_final / k0_3d;
        document.getElementById('val-3d-neff').textContent = nEff3D.toFixed(4);
        document.getElementById('val-3d-status').textContent = '도파 모드 존재';
        document.getElementById('val-3d-status').style.color = '#22c55e'; // Green
        
        // Solve structural waves parameters for profiles
        kz1_3d = Math.sqrt(k0_3d*k0_3d*nCore3D*nCore3D - beta_I*beta_I);
        gCover1_3d = Math.sqrt(beta_I*beta_I - k0_3d*k0_3d*nCover3D*nCover3D);
        gSub1_3d = Math.sqrt(beta_I*beta_I - k0_3d*k0_3d*nSub3D*nSub3D);
        phi01_3d = 0.5 * (Math.atan(gCover1_3d / kz1_3d) - Math.atan(gSub1_3d / kz1_3d));

        if (beta_II) {
            kz2_3d = Math.sqrt(k0_3d*k0_3d*nCore3D*nCore3D - beta_II*beta_II);
            gCover2_3d = Math.sqrt(beta_II*beta_II - k0_3d*k0_3d*nCover3D*nCover3D);
            gSub2_3d = Math.sqrt(beta_II*beta_II - k0_3d*k0_3d*nSub3D*nSub3D);
            phi02_3d = 0.5 * (Math.atan(gCover2_3d / kz2_3d) - Math.atan(gSub2_3d / kz2_3d));
        }
        
        ky_3d = Math.sqrt(k0_3d*k0_3d*nEff_I*nEff_I - beta_final*beta_final);
        gY_3d = Math.sqrt(beta_final*beta_final - k0_3d*k0_3d*nEff_II*nEff_II);
        
        precompute3DProfiles(selectedMode3D);
    } else {
        cutoff3D = true;
        document.getElementById('val-3d-neff').textContent = 'N/A';
        document.getElementById('val-3d-status').textContent = '차단됨 (Cutoff)';
        document.getElementById('val-3d-status').style.color = '#ef4444'; // Red
    }
    
    render3D();
}

function precompute3DProfiles(selectedMode) {
    const h_px = canvas3D.height;
    const w_px = canvas3D.width;
    
    G_I = new Array(h_px);
    G_II = new Array(h_px);
    F_y = new Array(w_px);
    
    // 1D vertical Region I profile
    for (let z_px = 0; z_px < h_px; z_px++) {
        let z_um = (260 - z_px) / 60.0;
        if (z_um >= 0 && z_um <= height3D) {
            G_I[z_px] = Math.cos(kz1_3d * z_um - phi01_3d);
        } else if (z_um < 0) {
            G_I[z_px] = Math.cos(phi01_3d) * Math.exp(gSub1_3d * z_um);
        } else {
            G_I[z_px] = Math.cos(kz1_3d * height3D - phi01_3d) * Math.exp(-gCover1_3d * (z_um - height3D));
        }
    }
    
    // 1D vertical Region II profile
    for (let z_px = 0; z_px < h_px; z_px++) {
        let z_um = (260 - z_px) / 60.0;
        if (wgType3D === 'rib' && kz2_3d) {
            if (z_um >= 0 && z_um <= slab3D) {
                G_II[z_px] = Math.cos(kz2_3d * z_um - phi02_3d);
            } else if (z_um < 0) {
                G_II[z_px] = Math.cos(phi02_3d) * Math.exp(gSub2_3d * z_um);
            } else {
                G_II[z_px] = Math.cos(kz2_3d * slab3D - phi02_3d) * Math.exp(-gCover2_3d * (z_um - slab3D));
            }
        } else {
            // For ridge/strip with no slab, exponential decay into air/cladding from substrate z=0
            if (z_um < 0) {
                G_II[z_px] = Math.exp(gSub1_3d * z_um);
            } else {
                G_II[z_px] = Math.exp(-gCover1_3d * z_um);
            }
        }
    }
    
    // 1D horizontal profile F(y) (which acts as F(x) in E(x, y) coordinates)
    for (let y_px = 0; y_px < w_px; y_px++) {
        let y_um = (y_px - 400) / 60.0;
        if (selectedMode === 1) {
            // Asymmetric first-order mode (two lobes of opposite sign)
            if (Math.abs(y_um) <= width3D / 2) {
                F_y[y_px] = Math.sin(ky_3d * y_um);
            } else {
                let edgeVal = Math.sin(ky_3d * width3D / 2);
                if (y_um > 0) {
                    F_y[y_px] = edgeVal * Math.exp(-gY_3d * (y_um - width3D / 2));
                } else {
                    F_y[y_px] = -edgeVal * Math.exp(-gY_3d * (-y_um - width3D / 2));
                }
            }
        } else {
            // Symmetric fundamental mode (single lobe)
            if (Math.abs(y_um) <= width3D / 2) {
                F_y[y_px] = Math.cos(ky_3d * y_um);
            } else {
                F_y[y_px] = Math.cos(ky_3d * width3D / 2) * Math.exp(-gY_3d * (Math.abs(y_um) - width3D / 2));
            }
        }
    }
}

function render3D() {
    const w = canvas3D.width;
    const h = canvas3D.height;
    
    if (cutoff3D) {
        ctx3D.fillStyle = '#040711';
        ctx3D.fillRect(0, 0, w, h);
        ctx3D.fillStyle = '#ef4444';
        ctx3D.font = 'bold 16px Inter';
        ctx3D.textAlign = 'center';
        ctx3D.fillText('도파로 차단 조건 (No Guided Mode)', w / 2, h / 2);
        
        // Draw blank cutlines graph
        cutCtx.fillStyle = '#040711';
        cutCtx.fillRect(0, 0, cutlinesCanvas.width, cutlinesCanvas.height);
        return;
    }
    
    const imgData = ctx3D.createImageData(w, h);
    const data = imgData.data;
    
    // 1. Draw 2D EIM mode profile heatmap (mapped from blue to red)
    for (let x = 0; x < w; x += 2) {
        let y_um = (x - 400) / 60.0;
        let F = F_y[x];
        
        for (let y = 0; y < h; y += 2) {
            let G = (Math.abs(y_um) <= width3D / 2) ? G_I[y] : G_II[y];
            let E = F * G;
            
            // Map absolute amplitude field to Jet colormap (max amplitude is 1.0, zero field is dark blue)
            let intensity = Math.abs(E);
            let col = getJetColor(intensity);
            
            draw2x2Block(data, x, y, w, col.r, col.g, col.b);
        }
    }
    ctx3D.putImageData(imgData, 0, 0);
    
    // 2. Draw Waveguide geometry boundary overlays (Yellow neon lines)
    ctx3D.strokeStyle = '#eab308';
    ctx3D.lineWidth = 2.0;
    ctx3D.beginPath();
    
    let x_left = 400 - (width3D / 2) * 60;
    let x_right = 400 + (width3D / 2) * 60;
    let y_sub = 260; // substrate y=0
    let y_top = 260 - height3D * 60;
    let y_wing = 260 - slab3D * 60;
    
    if (wgType3D === 'rib') {
        // Substrate horizontal line
        ctx3D.moveTo(0, y_sub); ctx3D.lineTo(w, y_sub);
        // Left wing top
        ctx3D.moveTo(0, y_wing); ctx3D.lineTo(x_left, y_wing);
        // Left rib vertical wall
        ctx3D.lineTo(x_left, y_top);
        // Rib center top
        ctx3D.lineTo(x_right, y_top);
        // Right rib vertical wall
        ctx3D.lineTo(x_right, y_wing);
        // Right wing top
        ctx3D.lineTo(w, y_wing);
    } else {
        // Ridge/Strip waveguide
        ctx3D.moveTo(0, y_sub); ctx3D.lineTo(x_left, y_sub);
        ctx3D.lineTo(x_left, y_top);
        ctx3D.lineTo(x_right, y_top);
        ctx3D.lineTo(x_right, y_sub);
        ctx3D.lineTo(w, y_sub);
    }
    ctx3D.stroke();
    
    // Draw 3D coordinate system ticks and units
    draw3DCoords(x_left, x_right, y_sub, y_top, y_wing);
    
    // 3. Render 1D Cutlines Graph
    render3DCutlines();
}

function draw3DCoords(x_left, x_right, y_sub, y_top, y_wing) {
    ctx3D.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    ctx3D.fillStyle = 'rgba(255, 255, 255, 0.6)';
    ctx3D.lineWidth = 1.0;
    ctx3D.font = '9px Fira Code, monospace';
    
    // X-axis (transverse horizontal axis x ticks at bottom)
    ctx3D.beginPath();
    ctx3D.moveTo(40, 340);
    ctx3D.lineTo(760, 340);
    ctx3D.stroke();
    
    const ticksX = [-5, -3, -1, 0, 1, 3, 5];
    ticksX.forEach(val => {
        let x_px = 400 + val * 60;
        ctx3D.beginPath();
        ctx3D.moveTo(x_px, 340);
        ctx3D.lineTo(x_px, 345);
        ctx3D.stroke();
        
        ctx3D.textAlign = 'center';
        ctx3D.fillText(val + 'μm', x_px, 355);
    });
    
    // Y-axis (transverse vertical axis y ticks at left)
    ctx3D.beginPath();
    ctx3D.moveTo(40, 40);
    ctx3D.lineTo(40, 340);
    ctx3D.stroke();
    
    const ticksY = [-1.0, 0, 1.0, 2.0, 3.0];
    ticksY.forEach(val => {
        let y_px = 260 - val * 60;
        ctx3D.beginPath();
        ctx3D.moveTo(35, y_px);
        ctx3D.lineTo(40, y_px);
        ctx3D.stroke();
        
        ctx3D.textAlign = 'right';
        ctx3D.textBaseline = 'middle';
        ctx3D.fillText(val.toFixed(1) + 'μm', 30, y_px);
    });
    
    // Text labels
    ctx3D.fillStyle = '#ffffff';
    ctx3D.font = 'bold 10px Inter';
    ctx3D.textAlign = 'left';
    ctx3D.fillText('y (수직 횡방향 축)', 10, 30);
    ctx3D.textAlign = 'right';
    ctx3D.fillText('x (수평 횡방향 축) →', 780, 325);
    
    // Overlay labels for Rib dimensions
    ctx3D.fillStyle = '#eab308';
    ctx3D.font = '9px Inter';
    ctx3D.textAlign = 'center';
    ctx3D.fillText(`W = ${width3D.toFixed(1)} μm`, 400, y_top - 6);
    
    ctx3D.textAlign = 'left';
    ctx3D.fillText(`H = ${height3D.toFixed(1)} μm`, x_right + 6, y_top + (y_sub - y_top)/2);
    if (wgType3D === 'rib') {
        ctx3D.fillText(`h = ${slab3D.toFixed(1)} μm`, x_right + 35, y_wing + (y_sub - y_wing)/2);
    }
}

function render3DCutlines() {
    const w = cutlinesCanvas.width;
    const h = cutlinesCanvas.height;
    
    cutCtx.fillStyle = '#040711';
    cutCtx.fillRect(0, 0, w, h);
    
    // Draw Grid Lines
    cutCtx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
    cutCtx.lineWidth = 1;
    for (let i = 0; i <= 2; i++) {
        let y_px = h * 0.15 + (i * h * 0.35);
        cutCtx.beginPath();
        cutCtx.moveTo(40, y_px);
        cutCtx.lineTo(w - 20, y_px);
        cutCtx.stroke();
    }
    
    // Draw horizontal ticks (-5 to 5 um)
    for (let val = -5; val <= 5; val += 2) {
        let x_px = 40 + ((val + 5) / 10) * (w - 60);
        cutCtx.beginPath();
        cutCtx.moveTo(x_px, h * 0.15);
        cutCtx.lineTo(x_px, h * 0.85);
        cutCtx.stroke();
        
        cutCtx.fillStyle = 'rgba(255,255,255,0.3)';
        cutCtx.font = '8px Fira Code';
        cutCtx.textAlign = 'center';
        cutCtx.fillText(`${val}um`, x_px, h * 0.94);
    }
    
    cutCtx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    gCtx.textAlign = 'right';
    cutCtx.font = '8px Fira Code';
    cutCtx.fillText('1.0', 25, h * 0.18);
    cutCtx.fillText('0.5', 25, h * 0.53);
    cutCtx.fillText('0.0', 25, h * 0.88);
    
    // Plot Cutlines
    // 1. Horizontal Cutline: |E(x, y = H/2)| -> plotted along graph in cyan
    cutCtx.strokeStyle = '#38bdf8'; // Cyan
    cutCtx.lineWidth = 2;
    cutCtx.beginPath();
    
    const center_y_px = Math.round(260 - (height3D / 2) * 60);
    const G_val = G_I[center_y_px];
    
    for (let x = 40; x < w - 20; x++) {
        let x_canvas = 400 + ((x - 40) / (w - 60) * 10 - 5) * 60; // maps graph coordinate to 3d canvas x pixel
        x_canvas = Math.max(0, Math.min(canvas3D.width - 1, Math.round(x_canvas)));
        
        let F = F_y[x_canvas];
        let E_val = F * G_val;
        
        let y_px = h * 0.85 - Math.abs(E_val) * (h * 0.7); // absolute field profile
        if (x === 40) cutCtx.moveTo(x, y_px);
        else cutCtx.lineTo(x, y_px);
    }
    cutCtx.stroke();
    
    // 2. Vertical Cutline: |E(x = 0, y)| -> plotted along graph in pink
    cutCtx.strokeStyle = '#f43f5e'; // Pink
    cutCtx.beginPath();
    
    const F_center = F_y[400]; // at x = 0
    
    for (let x = 40; x < w - 20; x++) {
        let y_um = ((x - 40) / (w - 60)) * 4.5 - 1.5; // maps graph coordinate to y [-1.5, 3.0] um
        let y_px = Math.round(260 - y_um * 60);
        y_px = Math.max(0, Math.min(canvas3D.height - 1, y_px));
        
        let G = G_I[y_px];
        let E_val = F_center * G;
        
        let y_px = h * 0.85 - Math.abs(E_val) * (h * 0.7);
        if (x === 40) cutCtx.moveTo(x, y_px);
        else cutCtx.lineTo(x, y_px);
    }
    cutCtx.stroke();
    
    // Legend
    cutCtx.font = 'bold 9px Inter';
    cutCtx.fillStyle = '#38bdf8';
    cutCtx.textAlign = 'left';
    cutCtx.fillText('수평 단면 프로파일 E(x, y=H/2)', w - 200, h * 0.25);
    cutCtx.fillStyle = '#f43f5e';
    cutCtx.fillText('수직 단면 프로파일 E(x=0, y)', w - 200, h * 0.38);
}

function bindParam3D(sliderId, inputId, applyBtnId, minVal, maxVal, isInt, updateVarFn) {
    const slider = document.getElementById(sliderId);
    const input = document.getElementById(inputId);
    const btn = document.getElementById(applyBtnId);
    
    slider.addEventListener('input', (e) => {
        let val = parseFloat(e.target.value);
        input.value = isInt ? val.toFixed(0) : val.toFixed(3);
        updateVarFn(val);
        updatePhysics3D();
    });
    
    function applyManualValue() {
        let val = parseFloat(input.value);
        if (isNaN(val)) return;
        
        val = Math.max(minVal, Math.min(maxVal, val));
        input.value = isInt ? val.toFixed(0) : val.toFixed(3);
        slider.value = Math.max(parseFloat(slider.min), Math.min(parseFloat(slider.max), val));
        
        updateVarFn(val);
        updatePhysics3D();
    }
    
    btn.addEventListener('click', applyManualValue);
    input.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') {
            applyManualValue();
        }
    });
}

function setupSliders3D() {
    bindParam3D('slider-3d-width', 'input-3d-width', 'btn-apply-3d-width', 0.5, 10.0, false, (v) => { width3D = v; });
    bindParam3D('slider-3d-height', 'input-3d-height', 'btn-apply-3d-height', 0.5, 5.0, false, (v) => { height3D = v; });
    bindParam3D('slider-3d-slab', 'input-3d-slab', 'btn-apply-3d-slab', 0.0, 4.0, false, (v) => { slab3D = v; });
    bindParam3D('slider-3d-lambda', 'input-3d-lambda', 'btn-apply-3d-lambda', 0.5, 2.5, false, (v) => { lambda3D = v; });
    bindParam3D('slider-3d-ncore', 'input-3d-ncore', 'btn-apply-3d-ncore', 1.0, 4.0, false, (v) => { nCore3D = v; });
    bindParam3D('slider-3d-ncover', 'input-3d-ncover', 'btn-apply-3d-ncover', 1.0, 4.0, false, (v) => { nCover3D = v; });
    bindParam3D('slider-3d-nsub', 'input-3d-nsub', 'btn-apply-3d-nsub', 1.0, 4.0, false, (v) => { nSub3D = v; });
    
    // Waveguide Type selection change listener
    document.getElementById('select-wg-type').addEventListener('change', (e) => {
        wgType3D = e.target.value;
        const slabGroup = document.getElementById('group-3d-slab');
        if (wgType3D === 'ridge') {
            slabGroup.style.display = 'none';
            slab3D = 0.0;
        } else {
            slabGroup.style.display = 'block';
            slab3D = parseFloat(document.getElementById('slider-3d-slab').value);
        }
        updatePhysics3D();
    });

    // Mode selection change listener
    document.getElementById('select-3d-mode').addEventListener('change', () => {
        updatePhysics3D();
    });
}

// ----------------------------------------------------
// PART 3: Simulation Loop Integration
// ----------------------------------------------------
function tick() {
    if (activeTab === 'coupled') {
        if (isPlaying) {
            time += 0.05 * speed;
        }
        let currentProbeX = (time * 10) % 400;
        drawElectricField(time, currentProbeX);
        drawPowerGraph(currentProbeX);
    }
    
    requestAnimationFrame(tick);
}

// Initialise application components
updatePhysics();
setupSliders();
setupControls();

setupSliders3D();
updatePhysics3D();

// Run simulator loop
playIcon.style.display = 'none';
pauseIcon.style.display = 'inline';
tick();
