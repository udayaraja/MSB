/*
 * General two-group (Cohen's d / two-sample t-test) power calculator.
 * Phase 1 of the MSB power analysis tool.
 *
 * Method: normal approximation to the noncentral t-distribution (Cohen, 1988),
 * equal n per group. This is a closed-form approximation, not an exact
 * noncentral-t computation -- adequate for typical n >~ 10-15 per group.
 *
 * Phase 2/3 (RNA-seq negative-binomial, metagenomics overdispersed-count
 * models) are NOT implemented here and will use precomputed lookup tables
 * generated offline in R (RNASeqPower / PROPER), interpolated client-side.
 * Planned dispersion input for RNA-seq mode: biological CV (not raw
 * dispersion parameter) -- decided 2026-09, not yet built.
 */

(function () {
  "use strict";

  // ---- Standard normal helpers -------------------------------------------

  // Standard normal CDF via Abramowitz & Stegun 26.2.17 approximation.
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

  // Inverse standard normal CDF (quantile function), Acklam's algorithm.
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

  // Power for a two-sample comparison, equal n per group.
  function powerTwoSampleT(n, d, alpha, tails) {
    const ncp = d * Math.sqrt(n / 2);
    if (tails === 2) {
      const zCrit = normInv(1 - alpha / 2);
      return (1 - normCDF(zCrit - ncp)) + normCDF(-zCrit - ncp);
    }
    const zCrit = normInv(1 - alpha);
    return 1 - normCDF(zCrit - ncp);
  }

  // Required n per group for a target power (closed-form approximation).
  function sampleSizeTwoSampleT(targetPower, d, alpha, tails) {
    if (d === 0) return Infinity;
    const zAlpha = tails === 2 ? normInv(1 - alpha / 2) : normInv(1 - alpha);
    const zBeta = normInv(targetPower);
    const n = 2 * Math.pow((zAlpha + zBeta) / Math.abs(d), 2);
    return Math.ceil(n);
  }

  // ---- DOM wiring ----------------------------------------------------------

  document.addEventListener("DOMContentLoaded", function () {
    const form = document.getElementById("ttestForm");
    if (!form) return; // JS loaded on a page without this calculator

    const directionButtons = document.querySelectorAll("#directionToggle button");
    const fieldN = document.getElementById("field-n");
    const fieldPower = document.getElementById("field-power");
    const effectInput = document.getElementById("effectInput");
    const fieldD = document.getElementById("field-d");
    const fieldMean1 = document.getElementById("field-mean1");
    const fieldMean2 = document.getElementById("field-mean2");
    const fieldSd = document.getElementById("field-sd");
    const resultBox = document.getElementById("resultBox");
    const resultHeadline = document.getElementById("resultHeadline");
    const resultDetail = document.getElementById("resultDetail");

    let direction = "power"; // "power" = solve for power, "n" = solve for sample size

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
        return parseFloat(document.getElementById("d").value);
      }
      const m1 = parseFloat(document.getElementById("mean1").value);
      const m2 = parseFloat(document.getElementById("mean2").value);
      const sd = parseFloat(document.getElementById("sd").value);
      if (sd <= 0) throw new Error("Pooled SD must be greater than 0.");
      return (m1 - m2) / sd;
    }

    form.addEventListener("submit", function (e) {
      e.preventDefault();
      resultBox.classList.add("hidden");

      try {
        const alpha = parseFloat(document.getElementById("alpha").value);
        const tails = parseInt(document.getElementById("tails").value, 10);
        const d = getEffectSize();

        if (!(alpha > 0 && alpha < 1)) throw new Error("Alpha must be between 0 and 1.");
        if (d === 0) throw new Error("Effect size is 0 -- power is undefined/uninformative.");

        if (direction === "power") {
          const n = parseInt(document.getElementById("n").value, 10);
          if (!(n >= 2)) throw new Error("Sample size per group must be at least 2.");
          const power = powerTwoSampleT(n, d, alpha, tails);
          resultHeadline.textContent = "Power = " + (power * 100).toFixed(1) + "%";
          resultDetail.innerHTML =
            "n = " + n + " per group, d = " + d.toFixed(3) +
            ", &alpha; = " + alpha + " (" + (tails === 2 ? "two-sided" : "one-sided") + ")";
        } else {
          const targetPower = parseFloat(document.getElementById("targetPower").value);
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
  });

  // Exposed for future unit testing / Phase 2+ reuse.
  window.MSBPower = {
    normCDF: normCDF,
    normInv: normInv,
    powerTwoSampleT: powerTwoSampleT,
    sampleSizeTwoSampleT: sampleSizeTwoSampleT
  };
})();
