/* qr.js — the interactive: build animation, role colouring, and the
   destructible Reed–Solomon demo. Depends on qrcode.js (window.QR). */
(function () {
  "use strict";
  const $ = (s) => document.querySelector(s);

  const canvas = $("#qr-canvas");
  const ctx = canvas.getContext("2d");
  const tip = $("#qr-tip");

  const ROLE = {
    finder:    { css: "--c-finder", label: "Finder pattern",   note: "corner locator" },
    separator: { css: "--c-finder", label: "Separator",        note: "quiet border" },
    timing:    { css: "--c-timing", label: "Timing pattern",   note: "the ruler" },
    alignment: { css: "--c-align",  label: "Alignment pattern",note: "distortion fix" },
    dark:      { css: "--c-dark",   label: "Dark module",      note: "always on" },
    format:    { css: "--c-format", label: "Format info",      note: "level + mask" },
    version:   { css: "--c-format", label: "Version info",     note: "size" },
    data:      { css: "--c-data",   label: "Data",             note: "your message" },
    ecc:       { css: "--c-ecc",    label: "Parity",           note: "error correction" },
    remainder: { css: "--c-remain", label: "Remainder",        note: "padding bits" },
  };
  function roleColor(role) {
    return getComputedStyle(document.documentElement).getPropertyValue(ROLE[role].css).trim() || "#888";
  }

  let state = {
    text: "https://moshemalka.com/lab",
    ecl: "M",
    mode: "roles",   // roles | scan
    tool: "inspect", // inspect | damage
    qr: null,
    reveal: 1,       // 0..1 build-animation progress
    cwRole: null,    // per-codeword role map for colouring data vs ecc
    moduleRole: null,// 2D role lookup for hover
    damaged: new Set(), // "r,c" of flipped modules
    anim: null,
  };

  const QUIET = 4;

  /* ---- build the code + derived lookups ---- */
  function rebuild(animate) {
    let qr;
    try {
      qr = QR.encode(state.text || " ", state.ecl);
    } catch (e) {
      // too long — bump to a level/size that fits, or truncate gracefully
      setVerdict(false, "Too long", "Trim the text — this exceeds versions 1–10 at level " + state.ecl + ".");
      return;
    }
    state.qr = qr;
    state.damaged.clear();

    // per-module role grid (data/ecc need codeword kind)
    const size = qr.size;
    const moduleRole = Array.from({ length: size }, () => new Array(size).fill(null));
    for (let r = 0; r < size; r++)
      for (let c = 0; c < size; c++) {
        const cell = qr.grid[r][c];
        moduleRole[r][c] = cell.type;
      }
    // upgrade data/remainder cells to data|ecc by codeword kind
    for (const m of qr.dataModules) {
      if (m.cw < 0) { moduleRole[m.r][m.c] = "remainder"; continue; }
      moduleRole[m.r][m.c] = qr.cwMeta[m.cw].kind; // "data" | "ecc"
    }
    state.moduleRole = moduleRole;

    // module lookup for damage/hover
    state.modByPos = {};
    for (const m of qr.dataModules) state.modByPos[m.r + "," + m.c] = m;

    updateFacts();
    buildLegend();
    updateVerdict();

    if (animate) startReveal();
    else { state.reveal = 1; draw(); }
  }

  function startReveal() {
    cancelAnimationFrame(state.anim);
    const start = performance.now();
    const dur = 1100;
    const step = (t) => {
      const p = Math.min(1, (t - start) / dur);
      state.reveal = easeOutCubic(p);
      draw();
      if (p < 1) state.anim = requestAnimationFrame(step);
    };
    state.reveal = 0;
    state.anim = requestAnimationFrame(step);
  }
  const easeOutCubic = (x) => 1 - Math.pow(1 - x, 3);

  /* ---- drawing ---- */
  function draw() {
    const qr = state.qr; if (!qr) return;
    const size = qr.size;
    const total = state.mode === "scan" ? size + QUIET * 2 : size;
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const cssSize = canvas.clientWidth || 420;
    canvas.width = Math.round(cssSize * dpr);
    canvas.height = Math.round(cssSize * dpr);
    const px = canvas.width / total;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const pxCss = (canvas.width / total) / dpr;

    // background
    ctx.fillStyle = state.mode === "scan" ? "#ffffff" : "#0d0f15";
    ctx.fillRect(0, 0, canvas.width / dpr, canvas.height / dpr);

    const off = state.mode === "scan" ? QUIET : 0;
    // reveal cutoff based on placement diagonal (r+c) for a nice wipe
    const maxRC = (size - 1) * 2;
    const cutoff = state.reveal * (maxRC + 1);

    for (let r = 0; r < size; r++) {
      for (let c = 0; c < size; c++) {
        const cell = qr.grid[r][c];
        const role = state.moduleRole[r][c];
        const shown = (r + c) <= cutoff || state.reveal >= 1;
        const x = (c + off) * pxCss, y = (r + off) * pxCss;

        if (state.mode === "scan") {
          ctx.fillStyle = cell.dark ? "#0a0b0e" : "#ffffff";
          if (shown && cell.dark) ctx.fillRect(x, y, pxCss + 0.5, pxCss + 0.5);
          continue;
        }
        // role mode: tint region; dark cells full colour, light cells faint
        const col = roleColor(role);
        if (!shown) {
          ctx.fillStyle = "rgba(255,255,255,0.02)";
          ctx.fillRect(x, y, pxCss + 0.5, pxCss + 0.5);
          continue;
        }
        if (cell.dark) {
          ctx.fillStyle = col;
          ctx.globalAlpha = 1;
        } else {
          ctx.fillStyle = col;
          ctx.globalAlpha = 0.13;
        }
        roundedCell(x, y, pxCss);
        ctx.globalAlpha = 1;

        // mark damaged cells with a subtle overlay in role mode
        if (state.damaged.has(r + "," + c)) {
          ctx.strokeStyle = "#ff5c5c";
          ctx.globalAlpha = 0.9;
          ctx.lineWidth = Math.max(1, pxCss * 0.12);
          ctx.beginPath();
          ctx.moveTo(x + pxCss * 0.2, y + pxCss * 0.2);
          ctx.lineTo(x + pxCss * 0.8, y + pxCss * 0.8);
          ctx.moveTo(x + pxCss * 0.8, y + pxCss * 0.2);
          ctx.lineTo(x + pxCss * 0.2, y + pxCss * 0.8);
          ctx.stroke();
          ctx.globalAlpha = 1;
        }
      }
    }
  }
  function roundedCell(x, y, s) {
    // tiny inset so cells read as a grid, not a blob
    const pad = s > 8 ? s * 0.06 : 0;
    ctx.fillRect(x + pad, y + pad, s - pad * 2 + 0.4, s - pad * 2 + 0.4);
  }

  /* ---- interaction: hover + damage ---- */
  function cellAt(evt) {
    const qr = state.qr; if (!qr) return null;
    const rect = canvas.getBoundingClientRect();
    const size = qr.size;
    const total = state.mode === "scan" ? size + QUIET * 2 : size;
    const px = rect.width / total;
    const off = state.mode === "scan" ? QUIET : 0;
    const clientX = (evt.touches ? evt.touches[0].clientX : evt.clientX);
    const clientY = (evt.touches ? evt.touches[0].clientY : evt.clientY);
    const c = Math.floor((clientX - rect.left) / px) - off;
    const r = Math.floor((clientY - rect.top) / px) - off;
    if (r < 0 || c < 0 || r >= size || c >= size) return null;
    return { r, c, px, off, rect };
  }

  function onHover(evt) {
    const hit = cellAt(evt);
    if (!hit || state.tool !== "inspect" || state.mode === "scan") { tip.style.opacity = 0; return; }
    const role = state.moduleRole[hit.r][hit.c];
    const meta = ROLE[role];
    const m = state.modByPos[hit.r + "," + hit.c];
    let extra = "";
    if (m && m.cw >= 0) {
      const cm = state.qr.cwMeta[m.cw];
      extra = ` · byte #${m.cw} · block ${cm.block + 1}`;
    }
    tip.innerHTML = `<b>${meta.label}</b> — ${meta.note}${extra}`;
    tip.style.left = ((hit.c + hit.off + 0.5) * hit.px) + "px";
    tip.style.top = ((hit.r + hit.off) * hit.px) + "px";
    tip.style.opacity = 1;
  }

  let painting = false;
  function onPaintStart(evt) {
    if (state.tool !== "damage") return;
    painting = true; paintAt(evt); evt.preventDefault();
  }
  function paintAt(evt) {
    if (!painting) return;
    const hit = cellAt(evt); if (!hit) return;
    const key = hit.r + "," + hit.c;
    if (state._painted && state._painted.has(key)) return;
    (state._painted = state._painted || new Set()).add(key);
    // flip the module
    state.qr.grid[hit.r][hit.c].dark = !state.qr.grid[hit.r][hit.c].dark;
    if (state.damaged.has(key)) state.damaged.delete(key); else state.damaged.add(key);
    draw();
    updateVerdict();
  }
  function onPaintEnd() { painting = false; state._painted = null; }

  /* ---- verdict + meter (the payoff) ---- */
  function correctableSymbols(qr) {
    // per block: floor(ecPerBlock/2); summed
    return qr.blocks.length * Math.floor(qr.ecPerBlock / 2);
  }
  function currentSymbolErrors(qr) {
    // compare live codewords to pristine, grouped: total symbol errors
    const live = QR.readCodewords(qr);
    let errs = 0;
    for (let i = 0; i < live.length; i++) if (live[i] !== qr.cleanCodewords[i]) errs++;
    return errs;
  }
  function updateVerdict() {
    const qr = state.qr; if (!qr) return;
    const cap = correctableSymbols(qr);
    const errs = currentSymbolErrors(qr);
    const dec = QR.decode(qr);
    // meter shows errors vs capacity
    const pct = cap ? Math.min(100, (errs / cap) * 100) : 0;
    $("#qr-meter").style.width = pct + "%";
    $("#qr-budget").innerHTML = `<b>${errs}</b> of <b>${cap}</b> recoverable symbols spent`;
    if (errs === 0) {
      setVerdict(true, "Pristine", "No damage. Decodes cleanly.");
    } else if (dec.ok && dec.text === qr.text) {
      setVerdict(true, "Recovered ✓", `${dec.totalCorrected} corrupted byte${dec.totalCorrected===1?"":"s"} detected and repaired. Still reads: “${trunc(dec.text)}”.`);
    } else if (dec.ok && dec.text) {
      setVerdict(false, "Misread", `Correction ran out of margin and locked onto the wrong message: “${trunc(dec.text)}”.`);
    } else {
      setVerdict(false, "Unrecoverable ✗", `Damage blew past the safety margin in at least one block. A scanner sees garbage.`);
    }
  }
  function setVerdict(ok, title, body) {
    const v = $("#qr-verdict");
    v.className = "verdict " + (ok ? "ok" : "dead");
    $("#qr-verdict-title").textContent = title;
    $("#qr-verdict-body").textContent = body;
  }
  const trunc = (s) => (s && s.length > 42 ? s.slice(0, 41) + "…" : s);

  /* ---- side panel ---- */
  function updateFacts() {
    const qr = state.qr;
    const dataBytes = qr.blocks.reduce((a, b) => a + b, 0);
    const eccBytes = qr.blocks.length * qr.ecPerBlock;
    $("#qr-facts").innerHTML =
      `Version <b>${qr.version}</b> · ${qr.size}×${qr.size} cells · mask <b>${qr.mask}</b><br>` +
      `<b>${dataBytes}</b> data bytes · <b>${eccBytes}</b> parity bytes` +
      (qr.blocks.length > 1 ? ` · <b>${qr.blocks.length}</b> blocks` : "");
  }
  function buildLegend() {
    const items = ["finder", "timing", "alignment", "format", "data", "ecc"];
    $("#qr-legend").innerHTML = items.map((k) =>
      `<span class="li"><span class="sw" style="background:${roleColor(k)}"></span>${ROLE[k].label}</span>`
    ).join("");
  }

  /* ---- wiring ---- */
  function bindChips(groupSel, key, cb) {
    const group = $(groupSel);
    group.addEventListener("click", (e) => {
      const chip = e.target.closest(".chip"); if (!chip) return;
      group.querySelectorAll(".chip").forEach((c) => c.setAttribute("aria-pressed", "false"));
      chip.setAttribute("aria-pressed", "true");
      cb(chip.dataset[key]);
    });
  }

  bindChips("#qr-ecl", "ecl", (v) => { state.ecl = v; rebuild(true); });
  bindChips("#qr-mode", "mode", (v) => { state.mode = v; tip.style.opacity = 0; draw(); });
  bindChips("#qr-tool", "tool", (v) => {
    state.tool = v;
    canvas.style.cursor = v === "damage" ? "cell" : "crosshair";
    $("#qr-hint").innerHTML = v === "damage"
      ? "Drag across the code to knock cells out. Aim for the <b style='color:var(--c-data)'>data</b> and <b style='color:var(--c-ecc)'>parity</b> zones."
      : "Hover a cell to see its job. Switch to <b>Damage</b> and drag across the code to knock cells out.";
    if (v === "damage") tip.style.opacity = 0;
  });

  $("#qr-text").addEventListener("input", debounce((e) => {
    state.text = e.target.value; rebuild(false);
  }, 180));
  $("#qr-build").addEventListener("click", () => rebuild(true));
  $("#qr-repair").addEventListener("click", () => { rebuild(false); });

  canvas.addEventListener("mousemove", onHover);
  canvas.addEventListener("mouseleave", () => (tip.style.opacity = 0));
  canvas.addEventListener("mousedown", onPaintStart);
  window.addEventListener("mousemove", (e) => { if (painting) paintAt(e); });
  window.addEventListener("mouseup", onPaintEnd);
  canvas.addEventListener("touchstart", onPaintStart, { passive: false });
  canvas.addEventListener("touchmove", (e) => { if (painting) { paintAt(e); e.preventDefault(); } }, { passive: false });
  canvas.addEventListener("touchend", onPaintEnd);

  window.addEventListener("resize", debounce(() => draw(), 120));

  function debounce(fn, ms) { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; }

  // go
  rebuild(true);
})();
