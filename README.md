# Anatomy of a QR code

Build a QR code cell by cell, then take a knife to it and watch it refuse to die.

**[▶ Live demo](https://moshejs.github.io/qr-code-anatomy/)**

A complete QR **encoder** and a real **Reed–Solomon decoder**, running entirely in
your browser — no libraries, no server. Type some text and watch the code assemble
module by module (finder patterns, the timing ruler, data, and parity). Then switch
to Damage mode, scribble across it, and watch the error correction detect and repair
the corrupted bytes live — until you finally push past the safety margin.

## Run it

Static files. Any server works:

```bash
python3 -m http.server 8000   # then open http://localhost:8000
```

## How it's built

Vanilla HTML/CSS/JS, **zero dependencies**. `qrcode.js` (~600 lines) is a from-scratch
codec: GF(256) arithmetic, Reed–Solomon encode *and* decode (syndromes →
Berlekamp–Massey → Chien search → Forney), byte-mode encoding, matrix construction,
masking, and BCH format/version info. The "still scans after damage" demo reads the
codewords back out of the damaged matrix and runs the real decoder on them.

Part of [The Lab](https://github.com/moshejs) — handmade interactive explainers.

## Credits and sources

- **[Reed–Solomon codes for coders](https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders)**
  by Wikiversity contributors, licensed under
  [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). The GF(256) and
  Reed–Solomon routines in `qrcode.js` (log/antilog tables, polynomial add/evaluate,
  generator polynomial, encoding by synthetic division, syndromes, the Berlekamp–Massey
  error locator and the Forney magnitude step) are adapted from the tutorial's Python
  reference code, translated to JavaScript and restricted to errors only (no erasures).
- **QR Code format:** the structure the encoder follows (error-correction block
  tables, alignment positions, mask patterns, format and version information) is
  defined by the ISO/IEC 18004 standard. The tables in `qrcode.js` hold values from
  that specification; this repository does not record which copy they were taken from.

## License

MIT
