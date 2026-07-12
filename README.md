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

## License

MIT
