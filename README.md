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
key equation by the extended Euclidean algorithm → Chien search → Forney),
byte-mode encoding, matrix construction,
masking, and BCH format/version info. The "still scans after damage" demo reads the
codewords back out of the damaged matrix and runs the real decoder on them.

## Tests

No dependencies, Node 18 or newer:

```bash
node --test
```

`test/qrcode.test.js` checks the codec against golden vectors in
`test/golden-vectors.json` (generator polynomials, check codewords, syndromes,
block decoding with 0 to nsym injected errors, and whole symbols for several texts,
levels and versions), plus field identities, a known-answer vector and randomized
error-correction properties.

Part of [The Lab](https://github.com/moshejs) — handmade interactive explainers.

## Credits and sources

- **Reed–Solomon routines:** earlier versions of `qrcode.js` (up to commit
  `57d03fe`) contained GF(256) and Reed–Solomon routines adapted from the Wikiversity
  tutorial [Reed–Solomon codes for coders](https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders)
  ([CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)). They were replaced
  on 2026-11-24 by an independent implementation based on the standard mathematics:
  the error-correction scheme of ISO/IEC 18004 (field GF(2^8) with primitive polynomial 0x11D,
  generator polynomial as the product of (x − α^i)) and the textbook treatment of
  Reed–Solomon decoding found in Lin & Costello, *Error Control Coding*, and Blahut,
  *Algebraic Codes for Data Transmission*: syndromes, the key equation solved by the
  extended Euclidean algorithm (Sugiyama, Kasahara, Hirasawa and Namekawa, 1975),
  Chien search, and Forney's formula with the formal derivative.
- **QR Code format:** the structure the encoder follows (error-correction block
  tables, alignment positions, mask patterns, format and version information) is
  defined by the ISO/IEC 18004 standard. The tables in `qrcode.js` hold values from
  that specification; this repository does not record which copy they were taken from.

## License

MIT
