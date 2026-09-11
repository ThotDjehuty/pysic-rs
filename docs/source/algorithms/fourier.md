# Fourier Analysis

Pysic-rs provides the discrete Fourier transform and its applications in spectral analysis —
power spectral density, Welch's averaged periodogram, and spectral derivatives — built on
the `rustfft` crate.

> **Python binding status:** `fft`, `ifft`, `spectral_derivative`, and `welch_psd` are bound.
> See the [API page](../api/fourier.md) for exact signatures.

---

## 1. The Discrete Fourier Transform

### 1.1 Definition

For $x_0,\ldots,x_{N-1}\in\mathbb{C}$, the **DFT** is

$$
\boxed{X_k = \sum_{n=0}^{N-1}x_n\,e^{-2\pi i\,nk/N}, \qquad k=0,\ldots,N-1}
$$

with the **inverse DFT**:

$$
x_n = \frac{1}{N}\sum_{k=0}^{N-1}X_k\,e^{+2\pi i\,nk/N}
$$

**Convention note:** these satisfy $\text{IDFT}(\text{DFT}(x))=x$ exactly (no factors of
$2\pi$). The physics convention uses $\hat{f}(\xi)=\int f(x)e^{-2\pi i\xi x}dx$; the
`fft_frequencies` routine returns the corresponding frequencies $f_k = k/(N\Delta t)$.

### 1.2 Key Properties (Proven)

1. **Linearity:** $\text{DFT}(\alpha x+\beta y) = \alpha\,\text{DFT}(x)+\beta\,\text{DFT}(y)$.

2. **Convolution theorem:** $\text{DFT}(x*y) = \text{DFT}(x)\cdot\text{DFT}(y)$, where
   $(x*y)_n = \sum_m x_m y_{n-m}$ is circular convolution.

3. **Parseval's theorem:** $\sum_n|x_n|^2 = \frac{1}{N}\sum_k|X_k|^2$ — the DFT is
   unitary up to a factor.

4. **Shift property:** if $y_n = x_{n-m}$, then $Y_k = e^{-2\pi imk/N}X_k$.

5. **Hermitian symmetry:** for real $x_n$, $X_{N-k} = X_k^*$.

### 1.3 The Fast Fourier Transform

The **Cooley–Tukey FFT** reduces the naive $O(N^2)$ cost to $O(N\log N)$ by exploiting
the recursive factorization of the DFT matrix. This is *exact* — the FFT computes the
same DFT up to rounding.

**Radix-2 decomposition:** split $x$ into even/odd indexed elements:

$$
X_k = E_k + e^{-2\pi ik/N}O_k, \qquad X_{k+N/2} = E_k - e^{-2\pi ik/N}O_k
$$

where $E_k,O_k$ are the $N/2$-point DFTs of the even and odd subsequences. Recurse:
cost $= 2\cdot(N/2)\log(N/2) + O(N) = O(N\log N)$.

**Mixed-radix:** `rustfft` supports any factorization $N = N_1 N_2 \cdots$, not just
powers of 2 — but power-of-2 lengths are fastest.

---

## 2. Power Spectral Density

### 2.1 Wiener–Khinchin Theorem

The **power spectral density** (PSD) $S_{xx}(f)$ of a wide-sense stationary process $x(t)$
is the Fourier transform of the autocorrelation:

$$
S_{xx}(f) = \int_{-\infty}^{\infty}R_{xx}(\tau)\,e^{-2\pi if\tau}\,d\tau,
\qquad R_{xx}(\tau) = \mathbb{E}[x(t)x^*(t+\tau)]
$$

### 2.2 Raw Periodogram

The **periodogram** estimate of the PSD is

$$
\hat{S}(f) = \frac{\Delta t}{N}|X(f)|^2
$$

where $X(f) = \sum_n x_n e^{-2\pi ifn\Delta t}$.

**Consistency caveat (proven):** the raw periodogram is **inconsistent** — its variance
does *not* vanish as $N\to\infty$. The reason is that each frequency bin estimates the
spectral density from a single data realization; the variance equals the squared spectral
density regardless of $N$.

### 2.3 Welch's Method

Welch's method restores consistency by **segment averaging**:

1. Split the signal into $K$ overlapping segments of length `nperseg`.
2. Window each segment: $x_k^{(w)} = x_k \cdot w_n$ (Hanning, Hamming, etc.).
3. Compute the periodogram of each windowed segment.
4. Average: $\hat{S}_w(f) = \frac{1}{K}\sum_{k=1}^K \hat{S}_k(f)$.

$$
\hat{S}_w(f) = \frac{1}{K}\sum_{k=1}^K\frac{\Delta t}{U}\left|X_k(f)\right|^2,
\qquad U = \frac{1}{nperseg}\sum_n w_n^2
$$

**Proven:** the window normalisation $U$ removes spectral bias; averaging over $K$
independent segments reduces estimator variance by $\sim 1/\sqrt{K}$.

**Bias-variance trade-off:** larger `nperseg` → finer frequency resolution (less bias)
but fewer segments (more variance), and vice versa.

---

## 3. Spectral Derivatives

### 3.1 Exactness in Fourier Space

Because $\widehat{f'} = 2\pi i\xi\,\hat{f}$, derivatives are exact (in the band-limited
sense) in Fourier space:

$$
\boxed{f^{(m)}(x) = \mathcal{F}^{-1}\!\bigl[(2\pi i\xi)^m\,\mathcal{F}(f)\bigr]}
$$

This is the **spectral (collocation) derivative** — exponentially accurate for analytic
periodic data, far exceeding finite differences.

**Error analysis:** for a function $f$ with Fourier coefficients $|\hat{f}_k| \le
C\,e^{-\alpha|k|}$ (analytic), the spectral derivative has error $O(e^{-\alpha N/2})$
(exponential convergence). For $C^\infty$ (non-analytic) functions, convergence is
algebraic: $O(N^{-p})$ where $p$ depends on the smoothness.

### 3.2 Connection to PDE Solvers

The spectral derivative is the spatial discretisation used by **split-step** and
**pseudospectral** methods (see [PDE](pde.md)). The advantage: a single FFT computes
all spatial derivatives simultaneously.

---

## 4. Routines

| Routine | Description |
|---------|-------------|
| `fft(real, imag)` / `ifft(real, imag)` | complex DFT / IDFT |
| `fft_real(x)` / `ifft_real(real, imag)` | real input helpers |
| `power_spectral_density(real, imag)` | periodogram $|X_k|^2/N$ |
| `fft_frequencies(n, dx)` | frequency axis $f_k = k/(N\Delta t)$ |
| `spectral_derivative(f, dx, order)` | $f^{(m)}$ via $(2\pi ik)^m$ multipliers |
| `welch_psd(x, fs, nperseg, noverlap)` | Welch averaged PSD → (freqs, psd) |
| `periodogram(x, fs)` | raw periodogram → (freqs, psd) |

---

## 5. Usage Examples

### FFT round-trip

```python
from pysicrs import fft, ifft
import math

n = 256
x = [math.sin(2*math.pi*3*k/n) for k in range(n)]
X = fft(x, [0.0]*n)

peak = max(range(n), key=lambda k: X[0][k]**2 + X[1][k]**2)
print(f"peak bin = {peak} (expect 3 or n-3)")

xr, xi = ifft(*X)
print(max(abs(a-b) for a,b in zip(x, xr)))  # ≈ 1e-12
```

### Spectral derivative accuracy

```python
from pysicrs import spectral_derivative
import math

n = 256
x = [math.sin(0.7*k) for k in range(n)]
d = spectral_derivative(x, dx=1.0, order=1)
err = max(abs(d[k] - 0.7*math.cos(0.7*k)) for k in range(n))
print(f"max error vs analytic: {err:.2e}")  # ≈ 1e-12
```

### Welch PSD of a noisy sine

```python
from pysicrs import welch_psd
import math, random

fs, n = 1000.0, 10000
x = [math.sin(2*math.pi*60.0*k/fs) + 0.1*(random.random()-0.5) for k in range(n)]
freqs, psd = welch_psd(x, fs, nperseg=512, noverlap=256)
peak_f = freqs[max(range(len(freqs)), key=lambda i: psd[i])]
print(f"dominant frequency = {peak_f:.1f} Hz (expect ~60)")
```

---

## 6. Numerical Notes

- The FFT length should be a smooth number (power of 2 is fastest).
- Un-windowed periodograms are inconsistent — always use `welch_psd` for statistical estimates.
- Derivative order multiplies by $(2\pi ik)^m$; high orders amplify high-frequency noise.
- For non-periodic data, apply a window (Hanning, Blackman) to reduce spectral leakage.

---

## 7. Advantages & Limitations

✅ $O(N\log N)$ exact DFT, exponentially accurate spectral derivatives

✅ Welch's method gives consistent PSD estimates

✅ No external FFTW dependency — self-contained `rustfft`

❌ Assumes uniform sampling and (for plain FFT) periodicity

❌ Raw periodogram variance is high — always segment-average for statistics

---

## 8. References

1. Cooley, J.W. & Tukey, J.W. (1965). "An algorithm for the machine calculation of complex Fourier series." *Math. Comp.* 19:297–301.
2. Welch, P.D. (1967). "The use of fast Fourier transform for the estimation of power spectra." *IEEE Trans. Audio Electroacoust.* AU-15:70–73.
3. Press, W. et al. (2007). *Numerical Recipes*, 3rd ed. Cambs.
4. Trefethen, L.N. (2000). *Spectral Methods in MATLAB*. SIAM.
5. Oppenheim, A.V. & Schafer, R.W. (2010). *Discrete-Time Signal Processing*, 3rd ed. Pearson.

---

## 9. Related Topics

- [PDE](pde.md) — split-step and spectral-Poisson solvers use these FFT kernels
- [Special Functions](special_functions.md) — windowing pairs with analytic kernels
- [Optimiz-rs](https://optimiz-rs.readthedocs.io/) — Welch PSD feeds time-series analyses
