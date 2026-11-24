/* ============================================================
   qrcode.js — a from-scratch QR encoder + Reed–Solomon decoder.
   No dependencies. Byte mode, versions 1–10, ECC levels L/M/Q/H.

   The encoder returns a fully annotated matrix: every module knows
   its role (finder / timing / alignment / format / data / ecc …) and,
   for data/ecc modules, which codeword bit it carries and what mask
   bit was XORed over it. That lets the UI read the codewords back
   out of a *damaged* matrix and run real error correction on them.

   RS decode = syndromes → key equation by extended Euclid (Sugiyama)
   → Chien search → Forney.
   GF(256), primitive element α=2, prim. polynomial 0x11D.
   ============================================================ */
(function (global) {
  "use strict";

  /* ---------------- GF(2^8) and Reed–Solomon ----------------
     Independent implementation (2026-11-24) of the standard mathematics, as
     presented in:
       - ISO/IEC 18004 (QR Code), error-correction section: field GF(2^8) with
         primitive polynomial x^8+x^4+x^3+x^2+1 (0x11D), primitive element α = 2,
         generator g(x) = (x − α^0)(x − α^1)…(x − α^(n−k−1)), check symbols =
         remainder of m(x)·x^(n−k) divided by g(x).
       - S. Lin & D. J. Costello, "Error Control Coding" (2nd ed.), ch. 7:
         syndromes, key equation, Euclidean-algorithm decoding, Chien search,
         Forney's formula.
       - R. E. Blahut, "Algebraic Codes for Data Transmission", ch. 7.
       - Y. Sugiyama, M. Kasahara, S. Hirasawa, T. Namekawa, "A method for
         solving key equation for decoding Goppa codes", Information and
         Control 27 (1975) — the extended-Euclid solution of the key equation.

     Internal convention: a polynomial is a plain array of coefficients in
     ASCENDING degree (p[i] multiplies x^i). The QR byte stream is the
     opposite (first codeword = highest power), so the exported functions
     below translate at the boundary.
     Decoding: syndromes → key equation by extended Euclid (Sugiyama) →
     Chien search → Forney with the formal derivative. Errors only. */

  const GF = (() => {
    const POLY = 0x11d;
    // carry-less "shift and add" product, reduced modulo POLY as it goes
    const slowMul = (a, b) => {
      let acc = 0;
      while (b) {
        if (b & 1) acc ^= a;
        b >>= 1;
        a <<= 1;
        if (a & 0x100) a ^= POLY;
      }
      return acc;
    };
    const antilog = new Array(255); // antilog[i] = α^i
    const logOf = new Array(256).fill(-1); // logOf[α^i] = i, logOf[0] = -1
    for (let i = 0, v = 1; i < 255; i++, v = slowMul(v, 2)) {
      antilog[i] = v;
      logOf[v] = i;
    }
    const mod255 = (e) => ((e % 255) + 255) % 255;
    const alpha = (e) => antilog[mod255(e)];
    const mul = (a, b) => (a && b ? antilog[(logOf[a] + logOf[b]) % 255] : 0);
    const inv = (a) => {
      if (!a) throw new RangeError("GF(256): 0 has no inverse");
      return antilog[mod255(-logOf[a])];
    };
    const div = (a, b) => mul(a, inv(b));
    const pow = (a, k) => (a ? alpha(logOf[a] * k) : k === 0 ? 1 : 0);
    return { alpha, mul, inv, div, pow, slowMul };
  })();

  // ---- ascending-order polynomial helpers ----
  const degreeOf = (p) => {
    let d = p.length - 1;
    while (d >= 0 && p[d] === 0) d--;
    return d; // −1 for the zero polynomial
  };
  const trimmed = (p) => p.slice(0, degreeOf(p) + 1);
  function polySum(p, q) {
    const out = new Array(Math.max(p.length, q.length)).fill(0);
    p.forEach((c, i) => { out[i] ^= c; });
    q.forEach((c, i) => { out[i] ^= c; });
    return out;
  }
  function polyProduct(p, q) {
    if (!p.length || !q.length) return [];
    const out = new Array(p.length + q.length - 1).fill(0);
    for (let i = 0; i < p.length; i++) {
      if (!p[i]) continue;
      for (let j = 0; j < q.length; j++) out[i + j] ^= GF.mul(p[i], q[j]);
    }
    return out;
  }
  // value of p at x (ascending order: walk from the top coefficient down)
  function polyAt(p, x) {
    let acc = 0;
    for (let i = p.length - 1; i >= 0; i--) acc = GF.mul(acc, x) ^ p[i];
    return acc;
  }
  // long division a = quotient·b + remainder, b ≠ 0
  function polyDivMod(a, b) {
    const db = degreeOf(b);
    const rem = trimmed(a);
    const lead = GF.inv(b[db]);
    const quot = new Array(Math.max(0, rem.length - db)).fill(0);
    for (let k = rem.length - 1 - db; k >= 0; k--) {
      const f = GF.mul(rem[k + db], lead);
      if (!f) continue;
      quot[k] = f;
      for (let j = 0; j <= db; j++) rem[k + j] ^= GF.mul(f, b[j]);
    }
    return { quot, rem: trimmed(rem) };
  }
  // formal derivative; over GF(2^m) the even-power terms vanish
  function polyDerivative(p) {
    const out = new Array(Math.max(0, p.length - 1)).fill(0);
    for (let i = 1; i < p.length; i += 2) out[i - 1] = p[i];
    return out;
  }

  // g(x) = ∏_{i=0}^{nsym−1} (x + α^i), ascending order
  function generatorAscending(nsym) {
    let g = [1];
    for (let i = 0; i < nsym; i++) g = polyProduct(g, [GF.alpha(i), 1]);
    return g;
  }

  // Exported: generator polynomial, highest power first (leading 1).
  function rsGenPoly(nsym) {
    return Uint8Array.from(generatorAscending(nsym).reverse());
  }

  // Exported: the nsym check codewords for `data` (highest power first).
  // Computed with a feedback shift register: each incoming byte plus the
  // register's top cell is the feedback term that is multiplied by g(x)
  // and folded into the register (division by g, keeping the remainder).
  function rsEncode(data, nsym) {
    const g = generatorAscending(nsym); // g[nsym] = 1
    const reg = new Array(nsym).fill(0); // reg[j] = coefficient of x^j of the remainder
    for (const byte of data) {
      const fb = byte ^ reg[nsym - 1];
      for (let j = nsym - 1; j > 0; j--) reg[j] = reg[j - 1] ^ GF.mul(fb, g[j]);
      reg[0] = GF.mul(fb, g[0]);
    }
    return Uint8Array.from(reg.reverse());
  }

  // Exported: S_j = r(α^j) for j = 0 … nsym−1, where r is the received
  // block read as a polynomial with its first byte as the highest power.
  function rsSyndromes(block, nsym) {
    const s = new Uint8Array(nsym);
    for (let j = 0; j < nsym; j++) {
      const x = GF.alpha(j);
      let acc = 0;
      for (let i = 0; i < block.length; i++) acc = GF.mul(acc, x) ^ block[i];
      s[j] = acc;
    }
    return s;
  }

  // Key equation Λ(x)·S(x) ≡ Ω(x) (mod x^nsym) solved by running Euclid's
  // algorithm on (x^nsym, S(x)) and stopping once the remainder's degree
  // drops below nsym/2 (Sugiyama et al.). Returns Λ normalised to Λ(0) = 1
  // and the matching Ω, or null when no valid locator exists.
  function solveKeyEquation(synd, nsym) {
    let rPrev = new Array(nsym + 1).fill(0);
    rPrev[nsym] = 1; // x^nsym
    let rCur = trimmed(Array.from(synd));
    let tPrev = [];
    let tCur = [1];
    while (2 * degreeOf(rCur) >= nsym) {
      const { quot, rem } = polyDivMod(rPrev, rCur);
      const tNext = polySum(tPrev, polyProduct(quot, tCur));
      rPrev = rCur; rCur = rem;
      tPrev = tCur; tCur = trimmed(tNext);
    }
    const c0 = tCur[0] || 0;
    if (!c0) return null;
    const scale = GF.inv(c0);
    return {
      locator: tCur.map((c) => GF.mul(c, scale)),
      evaluator: rCur.map((c) => GF.mul(c, scale)),
    };
  }

  // Chien search: test Λ(α^−e) for every exponent e of the block (0 … n−1).
  // Term j of Λ is carried along as λ_j·α^(−j·e) and stepped by α^(−j).
  function chienSearch(locator, n) {
    const terms = locator.slice();
    const steps = locator.map((_, j) => GF.alpha(-j));
    const hits = [];
    for (let e = 0; e < n; e++) {
      let sum = 0;
      for (let j = 0; j < terms.length; j++) sum ^= terms[j];
      if (sum === 0) hits.push(e);
      for (let j = 1; j < terms.length; j++) terms[j] = GF.mul(terms[j], steps[j]);
    }
    return hits;
  }

  // Exported: correct one block [data…, check…]. Returns
  // { msg: Uint8Array, corrected: number of bytes fixed (−1 on failure), ok }.
  function rsDecodeBlock(codewords, nsym) {
    const received = Uint8Array.from(codewords);
    const n = received.length;
    const synd = rsSyndromes(received, nsym);
    if (synd.every((s) => s === 0)) return { msg: received, corrected: 0, ok: true };
    const fail = { msg: received, corrected: -1, ok: false };

    const key = solveKeyEquation(synd, nsym);
    if (!key) return fail;
    const { locator, evaluator } = key;
    const nu = degreeOf(locator);
    if (nu < 1 || 2 * nu > nsym) return fail;
    const exps = chienSearch(locator, n);
    if (exps.length !== nu) return fail; // roots outside the block or repeated

    // Forney (first consecutive root α^0): Y = X·Ω(X⁻¹) / Λ'(X⁻¹), X = α^e
    const dLocator = polyDerivative(locator);
    const fixed = received.slice();
    for (const e of exps) {
      const xInv = GF.alpha(-e);
      const denom = polyAt(dLocator, xInv);
      if (!denom) return fail;
      const magnitude = GF.mul(GF.alpha(e), GF.div(polyAt(evaluator, xInv), denom));
      fixed[n - 1 - e] ^= magnitude;
    }
    if (rsSyndromes(fixed, nsym).some((s) => s !== 0)) return fail;
    return { msg: fixed, corrected: nu, ok: true };
  }

  // names kept for the low-level test export
  const gmul = GF.mul;
  const gpow = GF.pow;
  const ginv = GF.inv;

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
    _gf: { gmul, gpow, ginv }, _rs: { rsEncode, rsDecodeBlock, rsGenPoly, rsSyndromes, solveKeyEquation, chienSearch },
    _tables: { ECC, ALIGN }, dataCapacityBytes, chooseVersion,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = QR;
  global.QR = QR;
})(typeof window !== "undefined" ? window : globalThis);
