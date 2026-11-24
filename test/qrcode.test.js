// Tests for qrcode.js. No dependencies: run with `node --test`.
//
// golden-vectors.json was captured on 2026-11-24 by calling the PREVIOUS
// implementation (commit 57d03fe) as a black box: generator polynomials,
// check codewords, syndromes, block decoding with 0…nsym injected errors,
// and whole symbols (codewords, mask, module grid) for several texts,
// levels and versions. The rewrite must reproduce every one of them.
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const path = require("node:path");
const fs = require("node:fs");

const QR = require(path.join(__dirname, "..", "qrcode.js"));
const { rsEncode, rsDecodeBlock, rsGenPoly, rsSyndromes } = QR._rs;
const { gmul, gpow, ginv } = QR._gf;
const golden = JSON.parse(fs.readFileSync(path.join(__dirname, "golden-vectors.json"), "utf8"));

const hex = (a) => Buffer.from(Uint8Array.from(a)).toString("hex");
const unhex = (s) => Array.from(Buffer.from(s, "hex"));

// deterministic PRNG (mulberry32) so failures are reproducible
function prng(seed) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* ---------------- golden vectors (behaviour of the old code) ---------------- */

test("golden: generator polynomials", () => {
  for (const [n, want] of Object.entries(golden.genPoly)) assert.equal(hex(rsGenPoly(+n)), want, `nsym=${n}`);
});

test("golden: check codewords", () => {
  for (const v of golden.encode) assert.equal(hex(rsEncode(Uint8Array.from(unhex(v.data)), v.nsym)), v.ecc, `k=${v.k} nsym=${v.nsym}`);
});

test("golden: syndromes", () => {
  for (const v of golden.syndromes) assert.equal(hex(rsSyndromes(Uint8Array.from(unhex(v.msg)), v.nsym)), v.synd);
});

test("golden: block decoding with 0…nsym injected errors", () => {
  assert.ok(golden.decode.length > 700);
  for (const v of golden.decode) {
    const res = rsDecodeBlock(unhex(v.received), v.nsym);
    const label = `nsym=${v.nsym} errors=${v.errors}`;
    assert.equal(res.ok, v.ok, label);
    assert.equal(res.corrected, v.corrected, label);
    assert.equal(hex(res.msg), v.msg, label);
  }
});

test("golden: whole symbols are bit-identical", () => {
  for (const v of golden.symbols) {
    const q = QR.encode(v.text, v.level, v.forced || undefined);
    const label = `${JSON.stringify(v.text).slice(0, 24)} ${v.level} v${v.version}`;
    assert.equal(q.version, v.version, label);
    assert.equal(q.mask, v.mask, label);
    assert.equal(hex(q.cleanCodewords), v.codewords, label);
    assert.equal(q.grid.map((row) => row.map((m) => (m.dark ? "1" : "0")).join("")).join("/"), v.grid, label);
    const dec = QR.decode(q);
    assert.equal(dec.ok, true, label);
    assert.equal(dec.text, v.text, label);
  }
});

/* ---------------- independent checks (do not rely on the old code) ---------------- */

test("field: table multiply agrees with shift-and-add multiply for all pairs", () => {
  const slow = (a, b) => {
    let r = 0;
    for (; b; b >>= 1) {
      if (b & 1) r ^= a;
      a <<= 1;
      if (a & 0x100) a ^= 0x11d;
    }
    return r;
  };
  for (let a = 0; a < 256; a++) for (let b = 0; b < 256; b++) assert.equal(gmul(a, b), slow(a, b));
  for (let a = 1; a < 256; a++) assert.equal(gmul(a, ginv(a)), 1);
  assert.equal(gpow(2, 255), 1);
  assert.equal(gpow(2, 8), 0x1d);
  // α = 2 is primitive: its powers visit all 255 non-zero elements
  const seen = new Set();
  for (let i = 0; i < 255; i++) seen.add(gpow(2, i));
  assert.equal(seen.size, 255);
});

test("generator: g(α^i) = 0 for i < nsym, and α^nsym is not a root", () => {
  for (const n of [2, 7, 10, 13, 17, 22, 30]) {
    const g = Array.from(rsGenPoly(n)); // highest power first
    const at = (x) => g.reduce((acc, c) => gmul(acc, x) ^ c, 0);
    assert.equal(g.length, n + 1);
    assert.equal(g[0], 1);
    for (let i = 0; i < n; i++) assert.equal(at(gpow(2, i)), 0, `n=${n} i=${i}`);
    assert.notEqual(at(gpow(2, n)), 0);
  }
});

test("known answer: HELLO WORLD, version 1-M, alphanumeric mode", () => {
  // Data codewords built here from the alphanumeric-mode rules (mode 0010,
  // 9-bit count, pairs as 45a+b in 11 bits, terminator, 0xEC/0x11 pad).
  const table = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";
  const text = "HELLO WORLD";
  const bits = [];
  const put = (v, n) => { for (let i = n - 1; i >= 0; i--) bits.push((v >> i) & 1); };
  put(0b0010, 4);
  put(text.length, 9);
  for (let i = 0; i < text.length; i += 2) {
    if (i + 1 < text.length) put(45 * table.indexOf(text[i]) + table.indexOf(text[i + 1]), 11);
    else put(table.indexOf(text[i]), 6);
  }
  put(0, 4);
  while (bits.length % 8) bits.push(0);
  const data = [];
  for (let i = 0; i < bits.length; i += 8) data.push(parseInt(bits.slice(i, i + 8).join(""), 2));
  for (let p = 0; data.length < 16; p++) data.push(p % 2 ? 0x11 : 0xec);
  assert.deepEqual(data, [32, 91, 11, 120, 209, 114, 220, 77, 67, 64, 236, 17, 236, 17, 236, 17]);
  // Check codewords for this widely used worked example; the same values come
  // out of the previous implementation (black box) and of this one.
  assert.deepEqual(Array.from(rsEncode(data, 10)), [196, 35, 39, 119, 235, 215, 231, 226, 93, 23]);
});

test("property: encoded blocks have zero syndromes", () => {
  const R = prng(1);
  for (let trial = 0; trial < 300; trial++) {
    const n = 2 + Math.floor(R() * 29);
    const k = 1 + Math.floor(R() * (255 - n - 1));
    const data = Array.from({ length: k }, () => Math.floor(R() * 256));
    const block = data.concat(Array.from(rsEncode(data, n)));
    assert.ok(rsSyndromes(block, n).every((s) => s === 0));
  }
});

function corrupt(block, count, R) {
  const bad = block.slice();
  const pos = new Set();
  while (pos.size < count) pos.add(Math.floor(R() * block.length));
  for (const p of pos) bad[p] ^= 1 + Math.floor(R() * 255);
  return bad;
}

test("property: up to floor(nsym/2) random errors are always corrected", () => {
  const R = prng(2);
  for (let trial = 0; trial < 2000; trial++) {
    const n = [7, 10, 13, 15, 16, 17, 18, 20, 22, 24, 26, 28, 30][trial % 13];
    const k = 1 + Math.floor(R() * 120);
    const data = Array.from({ length: k }, () => Math.floor(R() * 256));
    const block = data.concat(Array.from(rsEncode(data, n)));
    const errs = Math.floor(R() * (Math.floor(n / 2) + 1));
    const res = rsDecodeBlock(corrupt(block, errs, R), n);
    assert.equal(res.ok, true, `n=${n} k=${k} errs=${errs}`);
    assert.equal(res.corrected, errs);
    assert.deepEqual(Array.from(res.msg), block);
  }
});

test("property: beyond capacity, the decoder fails or lands on a valid codeword — never returns garbage", () => {
  const R = prng(3);
  let failures = 0;
  for (let trial = 0; trial < 2000; trial++) {
    const n = [7, 10, 13, 16, 22, 26, 30][trial % 7];
    const k = 1 + Math.floor(R() * 80);
    const data = Array.from({ length: k }, () => Math.floor(R() * 256));
    const block = data.concat(Array.from(rsEncode(data, n)));
    const t = Math.floor(n / 2);
    const errs = Math.min(block.length, t + 1 + Math.floor(R() * t));
    const bad = corrupt(block, errs, R);
    const res = rsDecodeBlock(bad, n);
    if (!res.ok) {
      failures++;
      assert.equal(res.corrected, -1);
      assert.deepEqual(Array.from(res.msg), bad, "on failure the input is returned untouched");
      continue;
    }
    // A "success" here can only be a miscorrection to a different codeword
    // within distance t of what was received (undetectable by any decoder).
    assert.ok(rsSyndromes(res.msg, n).every((s) => s === 0));
    const dist = Array.from(res.msg).reduce((d, b, i) => d + (b !== bad[i]), 0);
    assert.equal(dist, res.corrected);
    assert.ok(res.corrected <= t);
  }
  assert.ok(failures > 1900, `expected nearly all to be reported as failures, got ${failures}/2000`);
});

test("symbol: damaged matrix is repaired by QR.decode up to the block capacity", () => {
  const R = prng(4);
  for (const level of "LMQH") {
    const q = QR.encode("https://moshejs.github.io/qr-code-anatomy/", level);
    const t = Math.floor(q.ecPerBlock / 2);
    for (let b = 0; b < q.blocks.length; b++) {
      const cws = q.cwMeta.map((m, i) => (m.block === b ? i : -1)).filter((i) => i >= 0);
      const hit = new Set();
      while (hit.size < t) hit.add(cws[Math.floor(R() * cws.length)]);
      for (const m of q.dataModules) {
        if (hit.has(m.cw) && m.bit === 0) q.grid[m.r][m.c].dark = !q.grid[m.r][m.c].dark;
      }
    }
    const dec = QR.decode(q);
    assert.equal(dec.ok, true, level);
    assert.equal(dec.text, q.text, level);
    assert.equal(dec.totalCorrected, t * q.blocks.length, level);
  }
});
