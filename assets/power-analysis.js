/*
 * Multi-Omics Power Analysis tool.
 *
 * Two statistical engines:
 *  1. Two-sample t-test / Cohen's d (normal approximation to noncentral t,
 *     Cohen 1988) -- used by Multi-Omics and Proteomics tabs.
 *  2. Negative-binomial Wald-test approximation (Hart et al. 2013,
 *     Bioinformatics; the closed-form method behind the RNASeqPower
 *     Bioconductor package) -- used by Transcriptomics and Metagenomics
 *     tabs. Metagenomics reuses this engine with no compositional/DM-specific
 *     correction; see the caveat text on that panel.
 *
 * Both engines optionally take a multiple-testing-corrected alpha, following
 * Tarazona et al. 2020 (Nat Commun 11:3092), Eq. 3:
 *   alpha* = (m1 * FDR) / ((m - m1) * (1 - FDR)),  m1 = m * (expected % DE)
 * where the entered "target FDR" plays the role of alpha in that equation.
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

  // ---- Engine 1: two-sample t-test / Cohen's d ----------------------------

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

  // ---- Engine 2: negative-binomial Wald test (Hart et al. 2013) ----------
  // Var(count) = mu + cv^2 * mu^2 (NB parameterized by biological CV).
  // Delta-method: Var(log count) ~= 1/mu + cv^2 per sample; comparing two
  // group log-means gives total variance 2*(1/mu + cv^2)/n.

  function powerNBWald(n, mu, cv, fc, alpha, tails) {
    const logFC = Math.log(fc);
    const se = Math.sqrt((2 * (1 / mu + cv * cv)) / n);
    const z = Math.abs(logFC) / se;
    if (tails === 2) {
      const zCrit = normInv(1 - alpha / 2);
      return (1 - normCDF(zCrit - z)) + normCDF(-zCrit - z);
    }
    const zCrit = normInv(1 - alpha);
    return 1 - normCDF(zCrit - z);
  }

  function sampleSizeNBWald(targetPower, mu, cv, fc, alpha, tails) {
    const logFC = Math.log(fc);
    if (logFC === 0) return Infinity;
    const zAlpha = tails === 2 ? normInv(1 - alpha / 2) : normInv(1 - alpha);
    const zBeta = normInv(targetPower);
    const n = 2 * (1 / mu + cv * cv) * Math.pow((zAlpha + zBeta) / Math.abs(logFC), 2);
    return Math.ceil(n);
  }

  // ---- Multiple-testing correction (Tarazona et al. 2020, Eq. 3) ---------

  function correctedAlpha(m, dePercentFraction, fdr) {
    const m1 = m * dePercentFraction;
    if (m1 >= m) throw new Error("Expected % DE must be less than 100%.");
    return (m1 * fdr) / ((m - m1) * (1 - fdr));
  }

  function getEffectiveAlpha(panel) {
    const checkbox = panel.querySelector('[data-role="mtcEnable"]');
    if (checkbox && checkbox.checked) {
      const m = parseFloat(panel.querySelector('[data-role="numFeatures"]').value);
      const dePercent = parseFloat(panel.querySelector('[data-role="dePercent"]').value);
      const fdr = parseFloat(panel.querySelector('[data-role="targetFDR"]').value);
      if (!(m > 1)) throw new Error("Number of features tested must be greater than 1.");
      if (!(dePercent > 0 && dePercent < 100)) throw new Error("Expected % DE must be between 0 and 100.");
      if (!(fdr > 0 && fdr < 1)) throw new Error("Target FDR must be between 0 and 1.");
      const alphaStar = correctedAlpha(m, dePercent / 100, fdr);
      return { alpha: alphaStar, corrected: true, m: m, dePercent: dePercent, fdr: fdr };
    }
    const alpha = parseFloat(panel.querySelector('[data-role="alpha"]').value);
    return { alpha: alpha, corrected: false };
  }

  function initMTC(panel) {
    const checkbox = panel.querySelector('[data-role="mtcEnable"]');
    if (!checkbox) return;
    const fieldAlpha = panel.querySelector('[data-role="field-alpha"]');
    const fieldMtc = panel.querySelector('[data-role="field-mtc"]');
    checkbox.addEventListener("change", function () {
      fieldAlpha.style.display = checkbox.checked ? "none" : "";
      fieldMtc.style.display = checkbox.checked ? "" : "none";
    });
  }

  function alphaDetailText(alphaInfo, tails) {
    const tailLabel = tails === 2 ? "two-sided" : "one-sided";
    if (alphaInfo.corrected) {
      return "&alpha;* = " + alphaInfo.alpha.toExponential(3) +
        " (FDR-adjusted, m=" + alphaInfo.m + ", " + alphaInfo.dePercent + "% DE, target FDR=" +
        alphaInfo.fdr + ", " + tailLabel + ")";
    }
    return "&alpha; = " + alphaInfo.alpha + " (" + tailLabel + ")";
  }

  // ---- Calculator 1: t-test / Cohen's d instance --------------------------

  function initTTestCalculator(panel) {
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

    initMTC(panel);

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
        const tails = parseInt(q("tails").value, 10);
        const alphaInfo = getEffectiveAlpha(panel);
        const alpha = alphaInfo.alpha;
        const d = getEffectSize();

        if (!(alpha > 0 && alpha < 1)) throw new Error("Alpha must be between 0 and 1.");
        if (d === 0) throw new Error("Effect size is 0 -- power is undefined/uninformative.");

        if (direction === "power") {
          const n = parseInt(q("n").value, 10);
          if (!(n >= 2)) throw new Error("Sample size per group must be at least 2.");
          const power = powerTwoSampleT(n, d, alpha, tails);
          resultHeadline.textContent = "Power = " + (power * 100).toFixed(1) + "%";
          resultDetail.innerHTML =
            "n = " + n + " per group, d = " + d.toFixed(3) + ", " + alphaDetailText(alphaInfo, tails);
        } else {
          const targetPower = parseFloat(q("targetPower").value);
          if (!(targetPower > 0 && targetPower < 1)) throw new Error("Target power must be between 0 and 1.");
          const n = sampleSizeTwoSampleT(targetPower, d, alpha, tails);
          resultHeadline.textContent = "n = " + n + " per group";
          resultDetail.innerHTML =
            "Target power = " + (targetPower * 100).toFixed(0) + "%, d = " + d.toFixed(3) + ", " +
            alphaDetailText(alphaInfo, tails);
        }

        resultBox.classList.remove("hidden");
      } catch (err) {
        resultHeadline.textContent = "Input error";
        resultDetail.textContent = err.message;
        resultBox.classList.remove("hidden");
      }
    });
  }

  // ---- Calculator 2: NB Wald instance (Transcriptomics / Metagenomics) ---

  function initNBCalculator(panel) {
    const form = panel.querySelector('[data-role="calcForm"]');
    if (!form) return;

    const q = (role) => panel.querySelector('[data-role="' + role + '"]');
    const directionButtons = panel.querySelectorAll('[data-role="directionToggle"] button');
    const fieldN = q("field-n");
    const fieldPower = q("field-power");
    const resultBox = q("resultBox");
    const resultHeadline = q("resultHeadline");
    const resultDetail = q("resultDetail");

    initMTC(panel);

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

    form.addEventListener("submit", function (e) {
      e.preventDefault();
      resultBox.classList.add("hidden");

      try {
        const tails = parseInt(q("tails").value, 10);
        const alphaInfo = getEffectiveAlpha(panel);
        const alpha = alphaInfo.alpha;
        const mu = parseFloat(q("mu").value);
        const cv = parseFloat(q("cv").value);
        const fc = parseFloat(q("fc").value);

        if (!(alpha > 0 && alpha < 1)) throw new Error("Alpha must be between 0 and 1.");
        if (!(mu > 0)) throw new Error("Mean count must be greater than 0.");
        if (!(cv >= 0)) throw new Error("Biological CV must be 0 or greater.");
        if (!(fc > 0) || fc === 1) throw new Error("Fold change must be > 0 and not equal to 1.");

        if (direction === "power") {
          const n = parseInt(q("n").value, 10);
          if (!(n >= 2)) throw new Error("Sample size per group must be at least 2.");
          const power = powerNBWald(n, mu, cv, fc, alpha, tails);
          resultHeadline.textContent = "Power = " + (Math.min(power, 1) * 100).toFixed(1) + "%";
          resultDetail.innerHTML =
            "n = " + n + " per group, &mu; = " + mu + ", CV = " + cv + ", fold change = " + fc +
            ", " + alphaDetailText(alphaInfo, tails);
        } else {
          const targetPower = parseFloat(q("targetPower").value);
          if (!(targetPower > 0 && targetPower < 1)) throw new Error("Target power must be between 0 and 1.");
          const n = sampleSizeNBWald(targetPower, mu, cv, fc, alpha, tails);
          resultHeadline.textContent = "n = " + n + " per group";
          resultDetail.innerHTML =
            "Target power = " + (targetPower * 100).toFixed(0) + "%, &mu; = " + mu + ", CV = " + cv +
            ", fold change = " + fc + ", " + alphaDetailText(alphaInfo, tails);
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

    const emptyState = document.getElementById("powerEmptyState");
    const tabButtons = sidebar.querySelectorAll(".power-tab");

    tabButtons.forEach(function (btn) {
      btn.addEventListener("click", function () {
        const target = btn.dataset.tab;

        tabButtons.forEach(function (b) { b.classList.remove("active"); });
        btn.classList.add("active");

        if (emptyState) emptyState.classList.remove("active");

        document.querySelectorAll(".power-panel").forEach(function (panel) {
          panel.classList.toggle("active", panel.id === "panel-" + target);
        });
      });
    });
  }

  // ---- Init -----------------------------------------------------------------

  document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll('.power-panel[data-calc="mo"], .power-panel[data-calc="prot"]').forEach(initTTestCalculator);
    document.querySelectorAll('.power-panel[data-calc="nb"]').forEach(initNBCalculator);
    initTabs();
  });

  // Exposed for future unit testing / reuse.
  window.MSBPower = {
    normCDF: normCDF,
    normInv: normInv,
    powerTwoSampleT: powerTwoSampleT,
    sampleSizeTwoSampleT: sampleSizeTwoSampleT,
    powerNBWald: powerNBWald,
    sampleSizeNBWald: sampleSizeNBWald,
    correctedAlpha: correctedAlpha
  };
})();
