# Quantile series and quantile-crossing processes

Companion repository for **Quantile Series: Finite-Sample Structure and Asymptotic Links to Quantile-Crossing Processes**, by Tianbo Chen, School of Big Data and Statistics, Anhui University.

## Paper

Quantile series (QSER) reconstruct a time-domain sequence from the complex coefficients underlying coefficient-based quantile periodograms (QP). Their ordinary periodograms equal QP. Across quantiles, these representations describe serial dependence that ordinary covariances can miss. The paper studies exact finite-sample identities, marginal information, and asymptotic links to scaled quantile-crossing processes. Its results include reconstruction bounds, Fourier limits, marginal identification, and consistency of covariance and spectral summaries under stated conditions.

The information comparison concerns marginal location and amplitude. QSER and crossings share the same normalized population spectral target. Averaging periodograms across independent records reduces sampling variation while leaving finite-record bias. It can assess spectral approximation; the paired scale comparison establishes the amplitude information absent from crossings.

Gaussian simulations examine reconstruction and spectral accuracy. Increasing marginal transformations separate amplitude from threshold dependence. A volatility model distinguishes processes with identical ordinary spectra through their tail spectra. An EEG application reports voltage-scale QSER amplitude and quantile-dependent spectral concentration.

[Download code and numerical summaries](dist/QSER_code.zip)

## Contents

| Path | Contents |
| --- | --- |
| `scripts/study.R` | Gaussian simulations, paired marginal transformations, and EEG analysis |
| `scripts/process_study.py` | Numerical summaries and the reconstruction, transformation, and EEG figures |
| `scripts/download_eeg.py` | Official EEG download, checksum verification, and O1 extraction in microvolts |
| `scripts/ensemble.R` | Volatility simulations, averaged periodograms, and paired scale calculations |
| `scripts/process_ensemble.py` | Analytical spectral targets, averaging risks, and averaging figures |
| `results/study/` | Gaussian and transformation summaries, EEG window statistics, and mean spectra |
| `results/ensemble/` | Mean spectra, averaging risks, marginal quantities, and numerical identity checks |
| `dist/` | Downloadable code and numerical summaries |

## Reproduction

R requires `quantreg`; `parallel` is distributed with R. Python dependencies are listed in `requirements.txt`. Run from the repository root:

```sh
Rscript -e 'install.packages("quantreg", repos="https://cloud.r-project.org")'
python -m pip install -r requirements.txt
python scripts/download_eeg.py
Rscript scripts/study.R
python scripts/process_study.py
Rscript scripts/ensemble.R
python scripts/process_ensemble.py
```

The scripts resolve paths from their own locations. Figures are written to `figures/`, and numerical results to `results/study/` and `results/ensemble/`. Replicate-level simulations and raw EEG data are generated locally and are not tracked. R calculations run sequentially on Windows.

## Experimental settings

The Gaussian designs use independent standard normal observations and a stationary AR(1) process with coefficient 0.6 and innovation standard deviation 0.8. Record lengths are 128, 256, and 512, with 400 independent records per model and length. Quantile levels are 0.1, 0.5, and 0.9. The paired transformation design uses 400 records of length 256 with normal, lognormal, and Student t margins with three degrees of freedom.

The volatility experiment uses a symmetric two-state chain with persistence parameters 0 and 0.8, conditional standard deviations 1/sqrt(5) and 3/sqrt(5), and standard normal innovations. Both models have unit ordinary spectrum. There are 400 records per model and length, with lengths 128 and 512. Periodograms are averaged in disjoint groups of 1, 5, 20, or 100 records. The scaled crossing comparator uses the population density. All sample sizes, persistence settings, and quantile levels are included in the saved comparisons. At 100 records per average, four groups contribute to each reported mean and Monte Carlo standard error.

## EEG data

The [EEG Motor Movement/Imagery Dataset on PhysioNet](https://doi.org/10.13026/C28G6P) supplies participant S001, baseline records R01 and R02, and channel O1 at 160 Hz. After the first five seconds, the next 8192 observations are divided into windows of 1024 or 512 observations and detrended separately. Spectra use Bartlett weights with maximum lag n/8 and denominator n/8 + 1. The 8–13 Hz shares integrate the fitted cosine polynomial over the exact band endpoints. Raw recordings are downloaded from PhysioNet and are not redistributed here.

## Citation

Chen, T. (2026). *Quantile Series: Finite-Sample Structure and Asymptotic Links to Quantile-Crossing Processes*. Working paper, Anhui University.
