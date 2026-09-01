/* Figures and widget for "What ranking does to a correlation".
 *
 * Everything on this page is drawn from methylation-subset.json, the committed
 * extract of GEO series GSE81211 built by prepare-data.py. No numbers are typed
 * in here.
 *
 * The widget's central trick: the readout always reports Pearson's r computed on
 * whatever is currently plotted. Switch the plot to ranks and that number slides
 * until it lands exactly on Spearman's rho -- because that is all Spearman is.
 */
(function () {
  "use strict";

  var DATA = null;
  var SCALE_LO = 0.65;                       // bottom of the heatmap colour scale
  var STOPS = [100, 75, 50, 25, 10, 5, 2, 1]; // "keep the most variable N%"

  /* ---------- small statistics helpers ---------- */

  function pearson(x, y) {
    var n = x.length, sx = 0, sy = 0, i;
    if (n < 3) return NaN;
    for (i = 0; i < n; i++) { sx += x[i]; sy += y[i]; }
    var mx = sx / n, my = sy / n, num = 0, dx2 = 0, dy2 = 0, a, b;
    for (i = 0; i < n; i++) {
      a = x[i] - mx; b = y[i] - my;
      num += a * b; dx2 += a * a; dy2 += b * b;
    }
    var den = Math.sqrt(dx2 * dy2);
    return den === 0 ? NaN : num / den;
  }

  // Average ranks, so tied beta values share a rank. 450K betas tie constantly.
  function ranks(v) {
    var n = v.length, order = new Array(n), out = new Array(n), i, j;
    for (i = 0; i < n; i++) order[i] = i;
    order.sort(function (a, b) { return v[a] - v[b]; });
    i = 0;
    while (i < n) {
      j = i;
      while (j + 1 < n && v[order[j + 1]] === v[order[i]]) j++;
      var shared = (i + j) / 2 + 1;
      for (var k = i; k <= j; k++) out[order[k]] = shared;
      i = j + 1;
    }
    return out;
  }

  function extent(v) {
    var lo = Infinity, hi = -Infinity;
    for (var i = 0; i < v.length; i++) {
      if (v[i] < lo) lo = v[i];
      if (v[i] > hi) hi = v[i];
    }
    return [lo, hi];
  }

  function fmt(v) { return isFinite(v) ? v.toFixed(3) : "--"; }

  /* ---------- static figure: the two correlation heatmaps ---------- */

  function shortCode(sample) {
    if (sample.group === "cancer-line") return "HCT";
    return (sample.group === "normal" ? "N" : "U") + sample.label.split(" ").pop();
  }

  function buildHeatmaps(host) {
    var samples = DATA.samples, n = samples.length;
    var codes = samples.map(shortCode);
    var panels = [
      ["beta_pearson", "Pearson <em>r</em> — the values as measured"],
      ["beta_spearman", "Spearman <em>ρ</em> — the same values, replaced by ranks"]
    ];

    panels.forEach(function (spec) {
      var m = DATA.matrices[spec[0]];
      var wrap = document.createElement("div");
      wrap.className = "rk-heat";

      var title = document.createElement("p");
      title.className = "rk-heat-title";
      title.innerHTML = spec[1];
      wrap.appendChild(title);

      var grid = document.createElement("div");
      grid.className = "rk-grid";
      grid.style.setProperty("--n", n);
      grid.setAttribute("role", "img");
      grid.setAttribute("aria-label",
        n + " by " + n + " matrix of " + spec[0].split("_")[1] +
        " correlations between colonic mucosa samples; the eleven mucosa samples are " +
        "a single flat colour and only the HCT116 cancer cell line stands out");

      var corner = document.createElement("span");
      corner.className = "rk-corner";
      grid.appendChild(corner);
      codes.forEach(function (c) {
        var el = document.createElement("span");
        el.className = "rk-collab";
        el.textContent = c;
        grid.appendChild(el);
      });

      for (var i = 0; i < n; i++) {
        var rowLab = document.createElement("span");
        rowLab.className = "rk-rowlab";
        rowLab.textContent = codes[i];
        grid.appendChild(rowLab);
        for (var j = 0; j < n; j++) {
          var v = m[i][j];
          var t = Math.max(0, Math.min(1, (v - SCALE_LO) / (1 - SCALE_LO)));
          var cell = document.createElement("span");
          cell.className = "rk-cell";
          cell.style.setProperty("--t", t.toFixed(3));
          cell.title = samples[i].label + " vs " + samples[j].label + ": " + v.toFixed(3);
          grid.appendChild(cell);
        }
      }
      wrap.appendChild(grid);
      host.appendChild(wrap);
    });
  }

  /* ---------- static figure: the distribution of beta values ---------- */

  function buildHistogram(host) {
    var h = DATA.histogram, frac = h.fraction;
    var peak = Math.max.apply(null, frac);
    var bars = document.createElement("div");
    bars.className = "rk-bars";
    bars.setAttribute("role", "img");
    bars.setAttribute("aria-label",
      "Histogram of beta values with two peaks, one near 0.06 and a larger one near " +
      "0.84, separated by a shallow valley");
    frac.forEach(function (f, i) {
      var bar = document.createElement("span");
      bar.className = "rk-bar";
      bar.style.setProperty("--h", (f / peak * 100).toFixed(1) + "%");
      bar.title = "beta " + h.edges[i].toFixed(2) + "–" + h.edges[i + 1].toFixed(2) +
        ": " + (f * 100).toFixed(2) + "% of measurements";
      bars.appendChild(bar);
    });
    host.appendChild(bars);
  }

  /* ---------- the widget ---------- */

  function buildWidget(root) {
    var canvas = root.querySelector("canvas");
    var ctx = canvas.getContext("2d");
    var selX = root.querySelector("#rk-x");
    var selY = root.querySelector("#rk-y");
    var rankBtns = root.querySelectorAll("[data-rank]");
    var scaleBtns = root.querySelectorAll("[data-scale]");
    var slider = root.querySelector("#rk-filter");
    var outR = root.querySelector("#rk-r");
    var outRho = root.querySelector("#rk-rho");
    var outN = root.querySelector("#rk-n");
    var outKept = root.querySelector("#rk-kept");

    var beta = DATA.beta, sd = DATA.sd;
    var bScale = DATA.beta_scale, sScale = DATA.sd_scale;
    var clip = DATA.clip;

    var state = { x: 0, y: 3, ranked: false, scale: "beta", stop: 0 };
    var anim = { t: 0, raf: null, timer: null };
    var view = null;   // recomputed whenever the selection or filter changes

    DATA.samples.forEach(function (s, i) {
      [selX, selY].forEach(function (sel) {
        var opt = document.createElement("option");
        opt.value = i;
        opt.textContent = s.label;
        sel.appendChild(opt);
      });
    });
    selX.value = state.x;
    selY.value = state.y;
    slider.max = STOPS.length - 1;
    slider.value = state.stop;

    function toM(b) {
      var c = Math.min(1 - clip, Math.max(clip, b));
      return Math.log2(c / (1 - c));
    }

    /* Recompute the point cloud: which CpGs survive the filter, where they sit in
       value space and in rank space, and the two coefficients. Ranks come from the
       beta values in every case -- the M-value transform is monotone, so it cannot
       reorder anything, and that is precisely the point the figure is making. */
    function recompute() {
      var pct = STOPS[state.stop];
      var keep;
      if (pct >= 100) {
        keep = sd.map(function (_, i) { return i; });
      } else {
        var sorted = sd.slice().sort(function (a, b) { return b - a; });
        var cut = sorted[Math.max(0, Math.ceil(sorted.length * pct / 100) - 1)];
        keep = [];
        for (var i = 0; i < sd.length; i++) if (sd[i] >= cut) keep.push(i);
      }

      var bx = [], by = [];
      for (var k = 0; k < keep.length; k++) {
        bx.push(beta[keep[k]][state.x] / bScale);
        by.push(beta[keep[k]][state.y] / bScale);
      }

      var vx = state.scale === "beta" ? bx : bx.map(toM);
      var vy = state.scale === "beta" ? by : by.map(toM);
      var rx = ranks(bx), ry = ranks(by);

      var n = keep.length;
      view = {
        n: n,
        kept: pct,
        // normalised to the unit square so the two spaces are directly comparable
        valX: norm(vx), valY: norm(vy),
        rnkX: rx.map(function (v) { return (v - 0.5) / n; }),
        rnkY: ry.map(function (v) { return (v - 0.5) / n; }),
        vx: vx, vy: vy, rx: rx, ry: ry,
        rho: pearson(rx, ry)
      };
    }

    function norm(v) {
      var e = extent(v), span = e[1] - e[0] || 1;
      return v.map(function (t) { return (t - e[0]) / span; });
    }

    function themeColours() {
      var cs = getComputedStyle(root);
      return {
        line: cs.getPropertyValue("--rk-line").trim() || "#d9e2e8",
        muted: cs.getPropertyValue("--rk-muted").trim() || "#63737a",
        accent: cs.getPropertyValue("--rk-accent").trim() || "#167d78"
      };
    }

    /* Match the backing store to the box CSS gave the canvas, at device resolution,
       so the axis labels are drawn at the size they are displayed. The element's own
       height is never written from here -- CSS sets the aspect ratio, and measuring
       a box this code had just resized would feed back into the layout. */
    function sizeCanvas() {
      var box = canvas.getBoundingClientRect();
      if (!box.width || !box.height) return false;
      var dpr = window.devicePixelRatio || 1;
      var w = Math.round(box.width * dpr);
      var h = Math.round(box.height * dpr);
      if (canvas.width === w && canvas.height === h) return false;
      canvas.width = w;
      canvas.height = h;
      return true;
    }

    function draw(t) {
      var dpr = window.devicePixelRatio || 1;
      var w = canvas.width, h = canvas.height, pad = 46 * dpr;
      var col = themeColours();
      ctx.clearRect(0, 0, w, h);

      var iw = w - pad * 1.4, ih = h - pad * 1.4;
      var x0 = pad, y0 = h - pad;

      ctx.strokeStyle = col.line;
      ctx.lineWidth = dpr;
      ctx.strokeRect(x0, y0 - ih, iw, ih);

      ctx.save();
      ctx.setLineDash([4 * dpr, 4 * dpr]);
      ctx.strokeStyle = col.line;
      ctx.beginPath();
      ctx.moveTo(x0, y0);
      ctx.lineTo(x0 + iw, y0 - ih);
      ctx.stroke();
      ctx.restore();

      // Small translucent dots: with a few thousand CpGs the density is the message.
      ctx.fillStyle = col.accent;
      ctx.globalAlpha = view.n > 1500 ? 0.16 : 0.42;
      var px, py, r = (view.n > 1500 ? 1.6 : 2.3) * dpr;
      for (var i = 0; i < view.n; i++) {
        px = x0 + (view.valX[i] + (view.rnkX[i] - view.valX[i]) * t) * iw;
        py = y0 - (view.valY[i] + (view.rnkY[i] - view.valY[i]) * t) * ih;
        ctx.beginPath();
        ctx.arc(px, py, r, 0, 6.2832);
        ctx.fill();
      }
      ctx.globalAlpha = 1;

      var unit = t > 0.5
        ? "rank among the " + view.n.toLocaleString() + " CpGs shown"
        : (state.scale === "beta" ? "methylation fraction (beta)" : "M-value");
      ctx.fillStyle = col.muted;
      ctx.font = (12 * dpr) + "px system-ui, sans-serif";
      ctx.textAlign = "center";
      ctx.fillText(DATA.samples[state.x].label + " — " + unit, x0 + iw / 2, h - 10 * dpr);
      ctx.save();
      ctx.translate(16 * dpr, y0 - ih / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(DATA.samples[state.y].label, 0, 0);
      ctx.restore();
    }

    /* r is computed on the interpolated coordinates, so during the transition the
       number visibly travels from Pearson to Spearman. */
    function readout(t) {
      var xs = new Array(view.n), ys = new Array(view.n);
      for (var i = 0; i < view.n; i++) {
        xs[i] = view.valX[i] + (view.rnkX[i] - view.valX[i]) * t;
        ys[i] = view.valY[i] + (view.rnkY[i] - view.valY[i]) * t;
      }
      var r = t >= 1 ? view.rho : pearson(xs, ys);
      outR.textContent = fmt(r);
      outRho.textContent = fmt(view.rho);
      outN.textContent = view.n.toLocaleString();
      outKept.textContent = view.kept >= 100
        ? "every CpG in the sample"
        : "the most variable " + view.kept + "%";
    }

    function render(t) {
      sizeCanvas();
      draw(t);
      readout(t);
    }

    var reduceMotion = window.matchMedia &&
      window.matchMedia("(prefers-reduced-motion: reduce)").matches;

    function settle(target) {
      if (anim.raf) { cancelAnimationFrame(anim.raf); anim.raf = null; }
      if (anim.timer) { clearTimeout(anim.timer); anim.timer = null; }
      anim.t = target;
      render(target);
    }

    /* The transition is the argument, so it is worth animating -- but a hidden tab
       stops firing rAF, and someone who switches away mid-slide must not come back
       to a half-ranked plot and a coefficient that matches neither view. The timer
       below snaps to the final state if the frame loop never finishes. */
    function animateTo(target) {
      if (anim.raf) cancelAnimationFrame(anim.raf);
      if (anim.timer) clearTimeout(anim.timer);
      if (reduceMotion || document.hidden) { settle(target); return; }

      var from = anim.t, dur = 750, start = null;
      function step(ts) {
        if (start === null) start = ts;
        var p = Math.min(1, (ts - start) / dur);
        var eased = p < 0.5 ? 2 * p * p : 1 - Math.pow(-2 * p + 2, 2) / 2;
        anim.t = from + (target - from) * eased;
        render(anim.t);
        if (p < 1) anim.raf = requestAnimationFrame(step);
        else settle(target);
      }
      anim.raf = requestAnimationFrame(step);
      anim.timer = setTimeout(function () { settle(target); }, dur + 120);
    }

    function refresh(animate) {
      recompute();
      if (animate) animateTo(state.ranked ? 1 : 0);
      else { anim.t = state.ranked ? 1 : 0; render(anim.t); }
    }

    function setPressed(btns, key, value) {
      btns.forEach(function (b) {
        b.setAttribute("aria-pressed", String(b.dataset[key] === value));
      });
    }

    selX.addEventListener("change", function () {
      state.x = +selX.value; refresh(false);
    });
    selY.addEventListener("change", function () {
      state.y = +selY.value; refresh(false);
    });
    rankBtns.forEach(function (b) {
      b.addEventListener("click", function () {
        var want = b.dataset.rank === "rank";
        if (want === state.ranked) return;
        state.ranked = want;
        setPressed(rankBtns, "rank", b.dataset.rank);
        refresh(true);
      });
    });
    scaleBtns.forEach(function (b) {
      b.addEventListener("click", function () {
        if (b.dataset.scale === state.scale) return;
        state.scale = b.dataset.scale;
        setPressed(scaleBtns, "scale", b.dataset.scale);
        refresh(false);
      });
    });
    slider.addEventListener("input", function () {
      state.stop = +slider.value; refresh(false);
    });

    // Quarto's light/dark switch flips a class on <body>; redraw into the new palette.
    new MutationObserver(function () { render(anim.t); })
      .observe(document.body, { attributes: true, attributeFilter: ["class"] });

    var resizeTimer = null;
    window.addEventListener("resize", function () {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(function () { render(anim.t); }, 120);
    });

    refresh(false);
  }

  /* ---------- boot ---------- */

  function init(data) {
    DATA = data;
    DATA.samples.forEach(function (s) { s.code = shortCode(s); });

    var heat = document.getElementById("rk-heatpair");
    if (heat) buildHeatmaps(heat);
    var hist = document.getElementById("rk-histogram");
    if (hist) buildHistogram(hist);
    var widget = document.getElementById("rk-widget");
    if (widget) buildWidget(widget);
  }

  function boot() {
    fetch("methylation-subset.json")
      .then(function (r) {
        if (!r.ok) throw new Error("HTTP " + r.status);
        return r.json();
      })
      .then(init)
      .catch(function (err) {
        console.error("[ranking] could not load methylation-subset.json", err);
        var w = document.getElementById("rk-widget");
        if (w) w.innerHTML =
          '<p class="rk-hint">The methylation extract could not be loaded, so the ' +
          'figures on this page are unavailable.</p>';
      });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
