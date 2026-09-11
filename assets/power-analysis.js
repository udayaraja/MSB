/*
 * Multi-Omics Power Analysis tool.
 * Phase 1: general two-group (Cohen's d / two-sample t-test) engine, reused
 * across the Multi-Omics and Proteomics tabs (same underlying model).
 *
 * Method: normal approximation to the noncentral t-distribution (Cohen, 1988),
 * equal n per group. Closed-form approximation, adequate for typical n >~ 10-15
 * per group -- not an exact noncentral-t computation.
 *
 * Transcriptomics (RNA-seq DE) and Metagenomics tabs are Phase 2/3, NOT
 * implemented here -- they need precomputed lookup tables generated offline
 * in R (RNASeqPower / PROPER for RNA-seq), interpolated client-side. Their
 * panels are static placeholders; see markup for details. Planned dispersion
 * input for RNA-seq mode: biological CV (decided 2026-09, not yet built).
 */

(function () {
  "use strict";

  // ---- Standard normal helpers -------------------------------------------

  function normCDF(z) {
    const b1 = 0.319381530;
    const b2 = -0.356563782;
    const b3 = 1.781477937;
    const b4 = -1.821255978;
    const b5 = 1.330274429;
    const p = 0.2316419;
    const c = 0.39894228;

    if (z >= 0) {
      const t = 1.0 / (1.0 + p * z);
      return 1.0 - c * Math.exp((-z * z) / 2.0) * t *
        (t * (t * (t * (t * b5 + b4) + b3) + b2) + b1);
    } else {
      return 1.0 - normCDF(-z);
    }
  }

  function normInv(p) {
    if (p <= 0 || p >= 1) {
      throw new RangeError("normInv: p must be in (0, 1)");
    }
    const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
               6.680131188771972e+01, -1.328068155288572e+01];
    const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
               3.754408661907416e+00];

    const pLow = 0.02425;
    const pHigh = 1 - pLow;
    let q, r;

    if (p < pLow) {
      q = Math.sqrt(-2 * Math.log(p));
      return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
             ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else if (p <= pHigh) {
      q = p - 0.5;
      r = q * q;
      return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
             (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
    } else {
      q = Math.sqrt(-2 * Math.log(1 - p));
      return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
              ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
  }

  // ---- Power / sample-size core ------------------------------------------

  function powerTwoSampleT(n, d, alpha, tails) {
    const ncp = d * Math.sqrt(n / 2);
    if (tails === 2) {
      const zCrit = normInv(1 - alpha / 2);
      return (1 - normCDF(zCrit - ncp)) + normCDF(-zCrit - ncp);
    }
    const zCrit = normInv(1 - alpha);
    return 1 - normCDF(zCrit - ncp);
  }

  function sampleSizeTwoSampleT(targetPower, d, alpha, tails) {
    if (d === 0) return Infinity;
    const zAlpha = tails === 2 ? normInv(1 - alpha / 2) : normInv(1 - alpha);
    const zBeta = normInv(targetPower);
    const n = 2 * Math.pow((zAlpha + zBeta) / Math.abs(d), 2);
    return Math.ceil(n);
  }

  // ---- Reusable calculator instance ---------------------------------------
  // Wires up one calculator (direction toggle + form + result box) scoped to
  // a container element, using data-role attributes instead of global IDs so
  // multiple independent instances can coexist on one page (e.g. Multi-Omics
  // and Proteomics tabs, both backed by the same engine).

  function initCalculator(panel) {
    const form = panel.querySelector('[data-role="calcForm"]');
    if (!form) return;

    const q = (role) => panel.querySelector('[data-role="' + role + '"]');
    const directionButtons = panel.querySelectorAll('[data-role="directionToggle"] button');
    const fieldN = q("field-n");
    const fieldPower = q("field-power");
    const effectInput = q("effectInput");
    const fieldD = q("field-d");
    const fieldMean1 = q("field-mean1");
    const fieldMean2 = q("field-mean2");
    const fieldSd = q("field-sd");
    const resultBox = q("resultBox");
    const resultHeadline = q("resultHeadline");
    const resultDetail = q("resultDetail");

    let direction = "power";

    directionButtons.forEach(function (btn) {
      btn.addEventListener("click", function () {
        directionButtons.forEach(function (b) { b.classList.remove("active"); });
        btn.classList.add("active");
        direction = btn.dataset.direction;
        fieldN.style.display = direction === "power" ? "" : "none";
        fieldPower.style.display = direction === "n" ? "" : "none";
        resultBox.classList.add("hidden");
      });
    });

    function syncEffectFields() {
      const useMeans = effectInput.value === "means";
      fieldD.style.display = useMeans ? "none" : "";
      fieldMean1.style.display = useMeans ? "" : "none";
      fieldMean2.style.display = useMeans ? "" : "none";
      fieldSd.style.display = useMeans ? "" : "none";
    }
    effectInput.addEventListener("change", syncEffectFields);
    syncEffectFields();

    function getEffectSize() {
      if (effectInput.value === "d") {
        return parseFloat(q("d").value);
      }
      const m1 = parseFloat(q("mean1").value);
      const m2 = parseFloat(q("mean2").value);
      const sd = parseFloat(q("sd").value);
      if (sd <= 0) throw new Error("Pooled SD must be greater than 0.");
      return (m1 - m2) / sd;
    }

    form.addEventListener("submit", function (e) {
      e.preventDefault();
      resultBox.classList.add("hidden");

      try {
        const alpha = parseFloat(q("alpha").value);
        const tails = parseInt(q("tails").value, 10);
        const d = getEffectSize();

        if (!(alpha > 0 && alpha < 1)) throw new Error("Alpha must be between 0 and 1.");
        if (d === 0) throw new Error("Effect size is 0 -- power is undefined/uninformative.");

        if (direction === "power") {
          const n = parseInt(q("n").value, 10);
          if (!(n >= 2)) throw new Error("Sample size per group must be at least 2.");
          const power = powerTwoSampleT(n, d, alpha, tails);
          resultHeadline.textContent = "Power = " + (power * 100).toFixed(1) + "%";
          resultDetail.innerHTML =
            "n = " + n + " per group, d = " + d.toFixed(3) +
            ", &alpha; = " + alpha + " (" + (tails === 2 ? "two-sided" : "one-sided") + ")";
        } else {
          const targetPower = parseFloat(q("targetPower").value);
          if (!(targetPower > 0 && targetPower < 1)) throw new Error("Target power must be between 0 and 1.");
          const n = sampleSizeTwoSampleT(targetPower, d, alpha, tails);
          resultHeadline.textContent = "n = " + n + " per group";
          resultDetail.innerHTML =
            "Target power = " + (targetPower * 100).toFixed(0) + "%, d = " + d.toFixed(3) +
            ", &alpha; = " + alpha + " (" + (tails === 2 ? "two-sided" : "one-sided") + ")";
        }

        resultBox.classList.remove("hidden");
      } catch (err) {
        resultHeadline.textContent = "Input error";
        resultDetail.textContent = err.message;
        resultBox.classList.remove("hidden");
      }
    });
  }

  // ---- Sidebar tab switching ----------------------------------------------

  function initTabs() {
    const sidebar = document.getElementById("powerSidebar");
    if (!sidebar) return;

    const tabButtons = sidebar.querySelectorAll(".power-tab");
    tabButtons.forEach(function (btn) {
      btn.addEventListener("click", function () {
        const target = btn.dataset.tab;

        tabButtons.forEach(function (b) { b.classList.remove("active"); });
        btn.classList.add("active");

        document.querySelectorAll(".power-panel").forEach(function (panel) {
          panel.classList.toggle("active", panel.id === "panel-" + target);
        });
      });
    });
  }

  // ---- Init -----------------------------------------------------------------

  document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll('.power-panel[data-calc="mo"], .power-panel[data-calc="prot"]').forEach(initCalculator);
    initTabs();
  });

  // Exposed for future unit testing / Phase 2+ reuse.
  window.MSBPower = {
    normCDF: normCDF,
    normInv: normInv,
    powerTwoSampleT: powerTwoSampleT,
    sampleSizeTwoSampleT: sampleSizeTwoSampleT
  };
})();
