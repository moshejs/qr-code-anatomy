/* ============================================================
   qrcode.js — a from-scratch QR encoder + Reed–Solomon decoder.
   No dependencies. Byte mode, versions 1–10, ECC levels L/M/Q/H.

   The encoder returns a fully annotated matrix: every module knows
   its role (finder / timing / alignment / format / data / ecc …) and,
   for data/ecc modules, which codeword bit it carries and what mask
   bit was XORed over it. That lets the UI read the codewords back
   out of a *damaged* matrix and run real error correction on them.

   RS decode = syndromes → Berlekamp–Massey → Chien search → Forney.
   GF(256), primitive element α=2, prim. polynomial 0x11D.
   ============================================================ */
(function (global) {
  "use strict";

  /* ---------------- GF(256) arithmetic ---------------- */
  const EXP = new Uint8Array(512);
  const LOG = new Uint8Array(256);
  (function initTables() {
    let x = 1;
    for (let i = 0; i < 255; i++) {
      EXP[i] = x;
      LOG[x] = i;
      x <<= 1;
      if (x & 0x100) x ^= 0x11d;
    }
    for (let i = 255; i < 512; i++) EXP[i] = EXP[i - 255];
  })();
  const gmul = (a, b) => (a === 0 || b === 0 ? 0 : EXP[LOG[a] + LOG[b]]);
  const gpow = (a, n) => (a === 0 ? 0 : EXP[(LOG[a] * n) % 255]);
  const ginv = (a) => EXP[255 - LOG[a]];

  // polynomials are big-endian: p[0] is the highest-degree coefficient.
  function polyMul(p, q) {
    const r = new Uint8Array(p.length + q.length - 1);
    for (let i = 0; i < p.length; i++)
      if (p[i]) for (let j = 0; j < q.length; j++) r[i + j] ^= gmul(p[i], q[j]);
    return r;
  }
  function polyEval(p, x) {
    let y = p[0];
    for (let i = 1; i < p.length; i++) y = gmul(y, x) ^ p[i];
    return y;
  }
  function polyScale(p, s) {
    const r = new Uint8Array(p.length);
    for (let i = 0; i < p.length; i++) r[i] = gmul(p[i], s);
    return r;
  }
  function polyAdd(a, b) {
    const r = new Uint8Array(Math.max(a.length, b.length));
    for (let i = 0; i < a.length; i++) r[i + r.length - a.length] = a[i];
    for (let i = 0; i < b.length; i++) r[i + r.length - b.length] ^= b[i];
    return r;
  }

  function rsGenPoly(nsym) {
    let g = Uint8Array.of(1);
    for (let i = 0; i < nsym; i++) g = polyMul(g, Uint8Array.of(1, EXP[i]));
    return g;
  }
  function rsEncode(data, nsym) {
    const gen = rsGenPoly(nsym);
    const res = new Uint8Array(data.length + nsym);
    res.set(data);
    for (let i = 0; i < data.length; i++) {
      const coef = res[i];
      if (coef !== 0)
        for (let j = 0; j < gen.length; j++) res[i + j] ^= gmul(gen[j], coef);
    }
    return res.slice(data.length);
  }

  /* ---- RS decode: recover a corrupted codeword block ---- */
  function rsSyndromes(msg, nsym) {
    const s = new Uint8Array(nsym);
    for (let i = 0; i < nsym; i++) s[i] = polyEval(msg, EXP[i]);
    return s;
  }
  function rsErrorLocator(synd, nsym) {
    let errLoc = Uint8Array.of(1);
    let oldLoc = Uint8Array.of(1);
    for (let i = 0; i < nsym; i++) {
      let delta = synd[i];
      for (let j = 1; j < errLoc.length; j++)
        delta ^= gmul(errLoc[errLoc.length - 1 - j], synd[i - j]);
      oldLoc = concat(oldLoc, Uint8Array.of(0));
      if (delta !== 0) {
        if (oldLoc.length > errLoc.length) {
          const newLoc = polyScale(oldLoc, delta);
          oldLoc = polyScale(errLoc, ginv(delta));
          errLoc = newLoc;
        }
        errLoc = polyAdd(errLoc, polyScale(oldLoc, delta));
      }
    }
    errLoc = stripLeadingZeros(errLoc);
    return errLoc;
  }
  function rsErrorPositions(errLoc, n) {
    // Λ(x) = ∏(1 - Y_p·x) with Y_p = α^(n-1-p); a position p is an error
    // when Λ(Y_p⁻¹) = 0. (Chien search, matched to our syndrome convention.)
    const errs = errLoc.length - 1;
    const pos = [];
    for (let p = 0; p < n; p++) {
      const Yp = gpow(2, n - 1 - p);
      if (polyEval(errLoc, ginv(Yp)) === 0) pos.push(p);
    }
    return pos.length === errs ? pos : null;
  }
  // Forney: given syndromes + error positions, compute & apply magnitudes.
  function rsCorrect(msg, synd, errPos) {
    const n = msg.length;
    // error locator from positions
    let q = Uint8Array.of(1);
    const X = [];
    for (const p of errPos) {
      const Xi = gpow(2, n - 1 - p);
      X.push(Xi);
      q = polyMul(q, Uint8Array.of(Xi, 1)); // (Xi*x + 1)
    }
    // syndrome polynomial (big-endian, high degree first)
    const sPoly = new Uint8Array(synd.length);
    for (let i = 0; i < synd.length; i++) sPoly[i] = synd[synd.length - 1 - i];
    // error evaluator Ω(x) = S(x)·Λ(x) mod x^nsym
    let omega = polyMul(sPoly, q);
    omega = omega.slice(Math.max(0, omega.length - synd.length));
    // Forney per position
    const out = msg.slice();
    for (let k = 0; k < X.length; k++) {
      const Xi = X[k];
      const XiInv = ginv(Xi);
      // formal derivative of Λ evaluated at Xi^-1 → product form
      let denom = 1;
      for (let j = 0; j < X.length; j++) {
        if (j === k) continue;
        denom = gmul(denom, 1 ^ gmul(XiInv, X[j]));
      }
      const num = polyEval(omega, XiInv);
      // e = Ω(Xᵢ⁻¹) / Λ'(Xᵢ⁻¹); with fcr=0 the Xᵢ factor cancels against Λ'.
      const magnitude = gmul(num, ginv(denom));
      out[errPos[k]] ^= magnitude;
    }
    return out;
  }
  // Returns {data, ecc, corrected, ok}. data/ecc are corrected codewords.
  function rsDecodeBlock(codewords, nsym) {
    const msg = Uint8Array.from(codewords);
    const synd = rsSyndromes(msg, nsym);
    let clean = true;
    for (const s of synd) if (s !== 0) { clean = false; break; }
    if (clean) return { msg, corrected: 0, ok: true };
    const errLoc = rsErrorLocator(synd, nsym);
    const pos = rsErrorPositions(errLoc, msg.length);
    if (!pos) return { msg, corrected: -1, ok: false }; // beyond capacity
    const fixed = rsCorrect(msg, synd, pos);
    // verify
    const s2 = rsSyndromes(fixed, nsym);
    for (const s of s2) if (s !== 0) return { msg, corrected: -1, ok: false };
    return { msg: fixed, corrected: pos.length, ok: true };
  }

  function concat(a, b) {
    const r = new Uint8Array(a.length + b.length);
    r.set(a); r.set(b, a.length); return r;
  }
  function stripLeadingZeros(p) {
    let i = 0;
    while (i < p.length - 1 && p[i] === 0) i++;
    return p.slice(i);
  }

  /* ---------------- QR error-correction tables (v1–10) ---------------- */
  // ECC[version][level] = [ecPerBlock, [[count, dataPer], ...]]  level L,M,Q,H
  const ECC = {
    1:  { L:[7,[[1,19]]],  M:[10,[[1,16]]],  Q:[13,[[1,13]]],  H:[17,[[1,9]]] },
    2:  { L:[10,[[1,34]]], M:[16,[[1,28]]],  Q:[22,[[1,22]]],  H:[28,[[1,16]]] },
    3:  { L:[15,[[1,55]]], M:[26,[[1,44]]],  Q:[18,[[2,17]]],  H:[22,[[2,13]]] },
    4:  { L:[20,[[1,80]]], M:[18,[[2,32]]],  Q:[26,[[2,24]]],  H:[16,[[4,9]]] },
    5:  { L:[26,[[1,108]]],M:[24,[[2,43]]],  Q:[18,[[2,15],[2,16]]], H:[22,[[2,11],[2,12]]] },
    6:  { L:[18,[[2,68]]], M:[16,[[4,27]]],  Q:[24,[[4,19]]],  H:[28,[[4,15]]] },
    7:  { L:[20,[[2,78]]], M:[18,[[4,31]]],  Q:[18,[[2,14],[4,15]]], H:[26,[[4,13],[1,14]]] },
    8:  { L:[24,[[2,97]]], M:[22,[[2,38],[2,39]]], Q:[22,[[4,18],[2,19]]], H:[26,[[4,14],[2,15]]] },
    9:  { L:[30,[[2,116]]],M:[22,[[3,36],[2,37]]], Q:[20,[[4,16],[4,17]]], H:[24,[[4,12],[4,13]]] },
    10: { L:[18,[[2,68],[2,69]]], M:[26,[[4,43],[1,44]]], Q:[24,[[6,19],[2,20]]], H:[28,[[6,15],[2,16]]] },
  };
  const ALIGN = {
    1:[], 2:[6,18], 3:[6,22], 4:[6,26], 5:[6,30],
    6:[6,34], 7:[6,22,38], 8:[6,24,42], 9:[6,26,46], 10:[6,28,50],
  };
  const ECL_BITS = { L:0b01, M:0b00, Q:0b11, H:0b10 };

  function blockStructure(version, level) {
    const [ecPerBlock, groups] = ECC[version][level];
    const blocks = [];
    let totalData = 0;
    for (const [count, dataPer] of groups)
      for (let i = 0; i < count; i++) { blocks.push(dataPer); totalData += dataPer; }
    return { ecPerBlock, blocks, totalData };
  }
  function dataCapacityBytes(version, level) {
    return blockStructure(version, level).totalData;
  }
  function charCountBits(version) { return version >= 10 ? 16 : 8; }

  function chooseVersion(byteLen, level) {
    for (let v = 1; v <= 10; v++) {
      const cap = dataCapacityBytes(v, level);
      const need = Math.ceil((4 + charCountBits(v) + 8 * byteLen) / 8);
      if (need <= cap) return v;
    }
    return null; // too big for v1–10
  }

  /* ---------------- bit buffer ---------------- */
  function BitBuffer() { this.bits = []; }
  BitBuffer.prototype.put = function (val, len) {
    for (let i = len - 1; i >= 0; i--) this.bits.push((val >>> i) & 1);
  };

  /* ---------------- encode: text → annotated matrix ---------------- */
  function encode(text, level, forcedVersion) {
    level = level || "M";
    const bytes = utf8Bytes(text);
    const version = forcedVersion || chooseVersion(bytes.length, level);
    if (!version) throw new Error("Message too long for versions 1–10 at level " + level);
    const { ecPerBlock, blocks, totalData } = blockStructure(version, level);

    // --- build data bit stream ---
    const bb = new BitBuffer();
    bb.put(0b0100, 4); // byte mode
    bb.put(bytes.length, charCountBits(version));
    for (const b of bytes) bb.put(b, 8);
    const capacityBits = totalData * 8;
    const term = Math.min(4, capacityBits - bb.bits.length);
    bb.put(0, term);
    while (bb.bits.length % 8 !== 0) bb.bits.push(0);
    const padBytes = [0xec, 0x11];
    let pi = 0;
    while (bb.bits.length < capacityBits) { bb.put(padBytes[pi & 1], 8); pi++; }

    // bits → data codewords
    const dataCw = [];
    for (let i = 0; i < bb.bits.length; i += 8) {
      let v = 0;
      for (let j = 0; j < 8; j++) v = (v << 1) | bb.bits[i + j];
      dataCw.push(v);
    }

    // split into blocks, compute ECC per block
    const dataBlocks = [], eccBlocks = [];
    let p = 0;
    for (const n of blocks) {
      const d = dataCw.slice(p, p + n); p += n;
      dataBlocks.push(d);
      eccBlocks.push(Array.from(rsEncode(Uint8Array.from(d), ecPerBlock)));
    }

    // interleave. Track, for each final codeword, which block it came from
    // and whether it is data or ecc — used for colouring + read-back.
    const finalCw = [];       // byte values
    const cwMeta = [];        // {block, kind:'data'|'ecc'}
    const maxData = Math.max(...blocks);
    for (let i = 0; i < maxData; i++)
      for (let b = 0; b < dataBlocks.length; b++)
        if (i < dataBlocks[b].length) { finalCw.push(dataBlocks[b][i]); cwMeta.push({ block: b, kind: "data" }); }
    for (let i = 0; i < ecPerBlock; i++)
      for (let b = 0; b < eccBlocks.length; b++) { finalCw.push(eccBlocks[b][i]); cwMeta.push({ block: b, kind: "ecc" }); }

    // final bit stream (MSB first)
    const dataBits = [];
    for (let k = 0; k < finalCw.length; k++)
      for (let j = 7; j >= 0; j--) dataBits.push({ bit: (finalCw[k] >> j) & 1, cw: k });

    // --- build the module matrix ---
    const size = version * 4 + 17;
    const grid = makeGrid(size, version);

    // reserve format & version areas so data placement skips them
    reserveFormat(grid, size);
    if (version >= 7) reserveVersion(grid, size);

    // place data bits in zig-zag, recording (r,c) → data-bit index
    const dataModules = placeData(grid, size, dataBits);

    // choose the best mask
    let best = null;
    for (let m = 0; m < 8; m++) {
      const trial = applyMask(grid, size, dataModules, m);
      const pen = penalty(trial, size);
      if (!best || pen < best.pen) best = { m, pen };
    }
    const mask = best.m;
    const masked = applyMask(grid, size, dataModules, mask);

    placeFormat(masked, size, level, mask);
    if (version >= 7) placeVersion(masked, size, version);

    return {
      version, level, mask, size,
      grid: masked,
      dataModules,       // [{r,c,cw,bit,maskBit}]
      cwMeta,            // per-codeword {block,kind}
      cleanCodewords: finalCw, // pristine codeword values (for a damage meter)
      blocks, ecPerBlock,
      text,
    };
  }

  /* ---- grid + function patterns ---- */
  function cell(type) { return { dark: false, type, fn: true }; }
  function makeGrid(size, version) {
    const g = [];
    for (let r = 0; r < size; r++) {
      const row = [];
      for (let c = 0; c < size; c++) row.push(null);
      g.push(row);
    }
    // finder patterns + separators
    const finders = [[0, 0], [size - 7, 0], [0, size - 7]];
    for (const [fr, fc] of finders) placeFinder(g, fr, fc, size);
    // timing patterns
    for (let i = 8; i < size - 8; i++) {
      if (g[6][i] === null) { g[6][i] = cell("timing"); g[6][i].dark = i % 2 === 0; }
      if (g[i][6] === null) { g[i][6] = cell("timing"); g[i][6].dark = i % 2 === 0; }
    }
    // alignment patterns
    const coords = ALIGN[version];
    for (const ar of coords)
      for (const ac of coords) {
        if ((ar === 6 && ac === 6) ||
            (ar === 6 && ac === size - 7) ||
            (ar === size - 7 && ac === 6)) continue;
        placeAlignment(g, ar, ac);
      }
    // dark module
    g[size - 8][8] = cell("dark"); g[size - 8][8].dark = true;
    return g;
  }
  function placeFinder(g, r, c, size) {
    for (let i = -1; i <= 7; i++)
      for (let j = -1; j <= 7; j++) {
        const rr = r + i, cc = c + j;
        if (rr < 0 || rr >= size || cc < 0 || cc >= size) continue;
        const inRing =
          (i >= 0 && i <= 6 && (j === 0 || j === 6)) ||
          (j >= 0 && j <= 6 && (i === 0 || i === 6));
        const inCore = i >= 2 && i <= 4 && j >= 2 && j <= 4;
        if (i >= 0 && i <= 6 && j >= 0 && j <= 6) {
          g[rr][cc] = cell("finder");
          g[rr][cc].dark = inRing || inCore;
        } else {
          g[rr][cc] = cell("separator");
          g[rr][cc].dark = false;
        }
      }
  }
  function placeAlignment(g, r, c) {
    for (let i = -2; i <= 2; i++)
      for (let j = -2; j <= 2; j++) {
        g[r + i][c + j] = cell("alignment");
        const ring = Math.max(Math.abs(i), Math.abs(j));
        g[r + i][c + j].dark = ring !== 1;
      }
  }
  function reserveFormat(g, size) {
    for (let i = 0; i < 9; i++) {
      if (g[8][i] === null) g[8][i] = cell("format");
      if (g[i][8] === null) g[i][8] = cell("format");
    }
    for (let i = 0; i < 8; i++) {
      if (g[8][size - 1 - i] === null) g[8][size - 1 - i] = cell("format");
      if (g[size - 1 - i][8] === null) g[size - 1 - i][8] = cell("format");
    }
  }
  function reserveVersion(g, size) {
    for (let i = 0; i < 6; i++)
      for (let j = 0; j < 3; j++) {
        g[i][size - 11 + j] = cell("version");
        g[size - 11 + j][i] = cell("version");
      }
  }

  /* ---- data placement (zig-zag) ---- */
  function placeData(g, size, dataBits) {
    const modules = [];
    let bitIdx = 0;
    let upward = true;
    for (let col = size - 1; col > 0; col -= 2) {
      if (col === 6) col--; // skip vertical timing column
      for (let i = 0; i < size; i++) {
        const row = upward ? size - 1 - i : i;
        for (let k = 0; k < 2; k++) {
          const c = col - k;
          if (g[row][c] !== null) continue; // function module
          const db = dataBits[bitIdx];
          const bit = db ? db.bit : 0; // remainder bits = 0
          g[row][c] = { dark: bit === 1, type: db ? "data" : "remainder", fn: false };
          modules.push({ r: row, c, cw: db ? db.cw : -1, bit, maskBit: 0 });
          bitIdx++;
        }
      }
      upward = !upward;
    }
    return modules;
  }

  /* ---- masking ---- */
  const MASK_FN = [
    (r, c) => (r + c) % 2 === 0,
    (r, c) => r % 2 === 0,
    (r, c) => c % 3 === 0,
    (r, c) => (r + c) % 3 === 0,
    (r, c) => (Math.floor(r / 2) + Math.floor(c / 3)) % 2 === 0,
    (r, c) => ((r * c) % 2) + ((r * c) % 3) === 0,
    (r, c) => (((r * c) % 2) + ((r * c) % 3)) % 2 === 0,
    (r, c) => (((r + c) % 2) + ((r * c) % 3)) % 2 === 0,
  ];
  // Returns a deep-ish clone of the grid with the mask applied to data modules.
  function applyMask(g, size, modules, m) {
    const out = g.map((row) => row.map((cell) => (cell ? { ...cell } : null)));
    const fn = MASK_FN[m];
    for (const mod of modules) {
      const flip = fn(mod.r, mod.c);
      const maskBit = flip ? 1 : 0;
      out[mod.r][mod.c].dark = (mod.bit ^ maskBit) === 1;
      // record maskBit on the module for later read-back (mutates shared list;
      // fine because we re-derive per final chosen mask)
      mod.maskBit = maskBit;
    }
    return out;
  }

  function penalty(g, size) {
    let pen = 0;
    // rule 1: runs of 5+ same colour
    for (let r = 0; r < size; r++) {
      let run = 1;
      for (let c = 1; c < size; c++) {
        if (g[r][c].dark === g[r][c - 1].dark) { run++; if (run === 5) pen += 3; else if (run > 5) pen++; }
        else run = 1;
      }
    }
    for (let c = 0; c < size; c++) {
      let run = 1;
      for (let r = 1; r < size; r++) {
        if (g[r][c].dark === g[r - 1][c].dark) { run++; if (run === 5) pen += 3; else if (run > 5) pen++; }
        else run = 1;
      }
    }
    // rule 2: 2x2 blocks
    for (let r = 0; r < size - 1; r++)
      for (let c = 0; c < size - 1; c++) {
        const v = g[r][c].dark;
        if (v === g[r][c + 1].dark && v === g[r + 1][c].dark && v === g[r + 1][c + 1].dark) pen += 3;
      }
    // rule 3: finder-like patterns
    const pat1 = [1,0,1,1,1,0,1,0,0,0,0];
    const pat2 = [0,0,0,0,1,0,1,1,1,0,1];
    const matches = (arr) => {
      let a = true, b = true;
      for (let i = 0; i < 11; i++) { if ((arr[i] ? 1 : 0) !== pat1[i]) a = false; if ((arr[i] ? 1 : 0) !== pat2[i]) b = false; }
      return a || b;
    };
    for (let r = 0; r < size; r++)
      for (let c = 0; c < size - 10; c++) {
        const seg = []; for (let i = 0; i < 11; i++) seg.push(g[r][c + i].dark);
        if (matches(seg)) pen += 40;
      }
    for (let c = 0; c < size; c++)
      for (let r = 0; r < size - 10; r++) {
        const seg = []; for (let i = 0; i < 11; i++) seg.push(g[r + i][c].dark);
        if (matches(seg)) pen += 40;
      }
    // rule 4: dark-module ratio
    let dark = 0;
    for (let r = 0; r < size; r++) for (let c = 0; c < size; c++) if (g[r][c].dark) dark++;
    const pct = (dark * 100) / (size * size);
    pen += Math.floor(Math.abs(pct - 50) / 5) * 10;
    return pen;
  }

  /* ---- format & version info ---- */
  function bchFormat(data) {
    let d = data << 10;
    const g = 0b10100110111;
    for (let i = 4; i >= 0; i--) if (d & (1 << (i + 10))) d ^= g << i;
    return ((data << 10) | d) ^ 0b101010000010010;
  }
  function bchVersion(version) {
    let d = version << 12;
    const g = 0b1111100100101;
    for (let i = 5; i >= 0; i--) if (d & (1 << (i + 12))) d ^= g << i;
    return (version << 12) | d;
  }
  function placeFormat(g, size, level, mask) {
    const data = (ECL_BITS[level] << 3) | mask;
    const bits = bchFormat(data);
    const arr = [];
    for (let i = 14; i >= 0; i--) arr.push((bits >> i) & 1);
    // around top-left
    const map1 = [[8,0],[8,1],[8,2],[8,3],[8,4],[8,5],[8,7],[8,8],[7,8],[5,8],[4,8],[3,8],[2,8],[1,8],[0,8]];
    for (let i = 0; i < 15; i++) { const [r, c] = map1[i]; g[r][c].dark = arr[i] === 1; g[r][c].type = "format"; }
    // split across top-right & bottom-left
    const map2 = [[size-1,8],[size-2,8],[size-3,8],[size-4,8],[size-5,8],[size-6,8],[size-7,8],
                  [8,size-8],[8,size-7],[8,size-6],[8,size-5],[8,size-4],[8,size-3],[8,size-2],[8,size-1]];
    for (let i = 0; i < 15; i++) { const [r, c] = map2[i]; g[r][c].dark = arr[i] === 1; g[r][c].type = "format"; }
  }
  function placeVersion(g, size, version) {
    const bits = bchVersion(version);
    for (let i = 0; i < 18; i++) {
      const bit = (bits >> i) & 1;
      const r = Math.floor(i / 3), c = i % 3;
      g[r][size - 11 + c].dark = bit === 1; g[r][size - 11 + c].type = "version";
      g[size - 11 + c][r].dark = bit === 1; g[size - 11 + c][r].type = "version";
    }
  }

  /* ---- read codewords back out of a (possibly damaged) grid ---- */
  function readCodewords(qr) {
    const { grid, dataModules } = qr;
    // reconstruct data-bit stream from module current state, unmasking
    const bits = new Array(dataModules.length);
    for (let i = 0; i < dataModules.length; i++) {
      const m = dataModules[i];
      const cur = grid[m.r][m.c].dark ? 1 : 0;
      bits[i] = cur ^ m.maskBit;
    }
    // group into codewords using cw index
    const totalCw = qr.cwMeta.length;
    const cw = new Array(totalCw).fill(0);
    const cwBitCount = new Array(totalCw).fill(0);
    for (let i = 0; i < dataModules.length; i++) {
      const m = dataModules[i];
      if (m.cw < 0) continue; // remainder bit
      cw[m.cw] = (cw[m.cw] << 1) | bits[i];
      cwBitCount[m.cw]++;
    }
    return cw;
  }

  // Decode a grid: read codewords, de-interleave into blocks, RS-correct each,
  // reassemble the message. Returns rich info for the UI.
  function decode(qr) {
    const cw = readCodewords(qr);
    const meta = qr.cwMeta;
    const nBlocks = qr.blocks.length;
    const ecPerBlock = qr.ecPerBlock;
    // regroup interleaved codewords back into per-block [data..., ecc...]
    const blocksData = qr.blocks.map(() => []);
    const blocksEcc = qr.blocks.map(() => []);
    for (let i = 0; i < cw.length; i++) {
      if (meta[i].kind === "data") blocksData[meta[i].block].push(cw[i]);
      else blocksEcc[meta[i].block].push(cw[i]);
    }
    let totalCorrected = 0, ok = true;
    const recoveredData = [];
    const perBlock = [];
    for (let b = 0; b < nBlocks; b++) {
      const full = blocksData[b].concat(blocksEcc[b]);
      const res = rsDecodeBlock(full, ecPerBlock);
      perBlock.push({ corrected: res.corrected, ok: res.ok });
      if (!res.ok) { ok = false; }
      else if (res.corrected > 0) totalCorrected += res.corrected;
      for (let k = 0; k < blocksData[b].length; k++) recoveredData.push(res.msg[k]);
    }
    let text = null;
    if (ok) { try { text = decodeBytes(recoveredData, qr.version); } catch (e) { text = null; } }
    return { ok, text, totalCorrected, perBlock };
  }

  // strip mode/count/padding from a data-codeword stream → utf8 text
  function decodeBytes(dataCw, version) {
    const bits = [];
    for (const b of dataCw) for (let j = 7; j >= 0; j--) bits.push((b >> j) & 1);
    let p = 0;
    const read = (n) => { let v = 0; for (let i = 0; i < n; i++) v = (v << 1) | bits[p++]; return v; };
    const mode = read(4);
    if (mode !== 0b0100) throw new Error("not byte mode");
    const ccBits = charCountBits(version); // 8 bits for v1–9, 16 for v10+
    const len = read(ccBits);
    const out = new Uint8Array(len);
    for (let i = 0; i < len; i++) out[i] = read(8);
    return utf8Decode(out);
  }

  /* ---- utf-8 helpers ---- */
  function utf8Bytes(str) {
    if (typeof TextEncoder !== "undefined") return Array.from(new TextEncoder().encode(str));
    return Array.from(unescape(encodeURIComponent(str)), (c) => c.charCodeAt(0));
  }
  function utf8Decode(bytes) {
    if (typeof TextDecoder !== "undefined") return new TextDecoder().decode(Uint8Array.from(bytes));
    return decodeURIComponent(escape(String.fromCharCode.apply(null, bytes)));
  }

  const QR = {
    encode, decode, readCodewords,
    // low-level, exposed for tests
    _gf: { gmul, gpow, ginv }, _rs: { rsEncode, rsDecodeBlock, rsGenPoly, rsSyndromes, rsErrorLocator, rsErrorPositions, rsCorrect },
    _tables: { ECC, ALIGN }, dataCapacityBytes, chooseVersion,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = QR;
  global.QR = QR;
})(typeof window !== "undefined" ? window : globalThis);
