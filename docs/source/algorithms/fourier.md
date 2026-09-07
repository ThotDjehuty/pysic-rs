# Fourier Analysis

Pysic-rs provides the discrete Fourier transform and its applications in spectral analysis —
power spectral density, Welch's averaged periodogram, and spectral derivatives — built on the
`rustfft` crate.

> **Python binding status:** `fft`, `ifft`, `spectral_derivative`, and `welch_psd` are bound.
> See the [API page](../api/fourier.md) for exact signatures.

---

## Mathematical Foundations

### The Discrete Fourier Transform (DFT)

For $x_0, \dots, x_{N-1} \in \mathbb{C}$,

$$
X_k = \sum_{n=0}^{N-1} x_n\, e^{-2\pi i\, nk/N}, \qquad k = 0, \dots, N-1
$$

with the inverse

$$
x_n = \frac{1}{N}\sum_{k=0}^{N-1} X_k\, e^{+2\pi i\, nk/N}
$$

**Proven properties:** linearity, the convolution theorem
$\mathcal{F}(x*y) = \mathcal{F}(x)\,\mathcal{F}(y)$, Parseval's theorem
$\sum |x_n|^2 = \frac1N\sum |X_k|^2$, and the shift/duality relations.

The **Fast Fourier Transform** reduces the naive $O(N^2)$ cost to $O(N\log N)$ via the
Cooley–Tukey decomposition (proven exact — the FFT computes the same DFT up to rounding).

### Power Spectral Density

The periodogram estimate of the power spectral density (Wiener–Khinchin theorem — the PSD
is the Fourier transform of the autocorrelation):

$$
\hat{S}(f) = \frac{\Delta t}{N}\left|X(f)\right|^2
$$

**Consistency caveat (proven):** the raw periodogram is *inconsistent* — its variance does
not vanish as $N\to\infty$. Welch's method restores consistency by segment averaging.

### Welch's Method (segment-averaged periodogram)

Split $x$ into overlapping segments of length `nperseg` (overlap `noverlap`), window each,
FFT, average the squared magnitudes:

$$
\hat S_w(f) = \frac{1}{K}\sum_{k=1}^{K} \frac{\Delta t}{U}\left|X_k(f)\right|^2,
\qquad U = \frac1{nperseg}\sum_n w_n^2
$$

**Proven:** the window normalization $U$ removes spectral-bias from the window; averaging
over $K$ independent segments reduces estimator variance by $\sim 1/\sqrt{K}$.

### Spectral Derivative

Because $\widehat{f'} = i k\, \hat f$, derivatives are exact (in the band-limited sense) in
Fourier space:

$$
f'(x) = \mathcal F^{-1}\!\big[i k\, \mathcal F(f)\big]
$$

This is the *spectral* (collocation) derivative — exponentially accurate for analytic
periodic data, far exceeding finite differences.

---

## Routines

| Routine | Description |
|---------|-------------|
| `fft(real, imag)` / `ifft(real, imag)` | complex DFT / IDFT (returns real+imag) |
| `fft_real(x)` / `ifft_real(real, imag)` | real input helpers |
| `power_spectral_density(real, imag)` | periodogram $|X_k|^2/N$ |
| `fft_frequencies(n, dx)` | frequency axis $f_k = k/(N\Delta t)$ |
| `spectral_derivative(f, dx, order)` | $f^{(n)}$ via $k^n$ multipliers |
| `welch_psd(x, fs, nperseg, noverlap)` | Welch averaged PSD → (freqs, psd) |
| `periodogram(x, fs)` | raw periodogram → (freqs, psd) |

---

## Usage Examples

### FFT round-trip

```python
from pysicrs import fft, ifft
import math

n = 256
x = [math.sin(2 * math.pi * 3 * k / n) for k in range(n)]   # frequency 3
X = fft(x, [0.0] * n)

peak = max(range(n), key=lambda k: X[0][k]**2 + X[1][k]**2)
print(f"peak bin = {peak} (expect 3 or n-3)")

xr, xi = ifft(*X)
print(max(abs(a - b) for a, b in zip(x, xr)))  # ≈ 1e-12
```

### Spectral derivative accuracy

```python
from pysicrs import spectral_derivative
import math

n = 256
x = [math.sin(0.7 * k) for k in range(n)]
d = spectral_derivative(x, dx=1.0, order=1)
err = max(abs(d[k] - 0.7 * math.cos(0.7 * k)) for k in range(n))
print(f"max error vs analytic: {err:.2e}")  # ≈ 1e-12
```

### Welch PSD of a noisy sine

```python
from pysicrs import welch_psd
import math, random

fs, n = 1000.0, 10000
x = [math.sin(2 * math.pi * 60.0 * k / fs) + 0.1 * (random.random() - 0.5) for k in range(n)]
freqs, psd = welch_psd(x, fs, nperseg=512, noverlap=256)
peak_f = freqs[max(range(len(freqs)), key=lambda i: psd[i])]
print(f"dominant frequency = {peak_f:.1f} Hz (expect ~60)")
```

---

## Numerical Notes

- The FFT length should be a smooth number (power of 2 is fastest, but any fine-grain
  factorization works with rustfft's mixed-radix kernels).
- Un-windowed periodograms are inconsistent — prefer `welch_psd` for statistical estimates.
- Derivative order multiplies by $(ik)^m$; high orders amplify high-frequency noise.

---

## Advantages & Limitations

✅ $O(N\log N)$ exact DFT, exponentially accurate spectral derivatives

✅ Welch's method gives consistent PSD estimates

✅ No external FFTW dependency — self-contained `rustfft`

❌ Assumes uniform sampling and (for plain FFT) periodicity

❌ Raw periodogram variance is high — always segment-average for statistics

---

## References

1. Cooley, J.W. & Tukey, J.W. (1965). "An algorithm for the machine calculation of complex Fourier series." *Math. Comp.* 19:297–301.
2. Welch, P.D. (1967). "The use of fast Fourier transform for the estimation of power spectra." *IEEE Trans. Audio Electroacoust.* AU-15:70–73.
3. Press, W. et al. (2007). *Numerical Recipes*, 3rd ed. Cambs.
4. Trefethen, L.N. (2000). *Spectral Methods in MATLAB*. SIAM.

---

## Related Topics

- [PDE](pde.md) – split-step and spectral-Poisson solvers need these FFT kernels
- [Special Functions](special_functions.md) – windowing pairs with analytic kernels
- [Statistics (Optimiz-rs)](https://optimiz-rs.readthedocs.io/) – Welch PSD feeds time-series analyses