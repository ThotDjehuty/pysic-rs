# API Reference: Fourier Analysis

The Python surface exposes Fourier routines at top level (`from pysicrs import fft`). See
also [Algorithms: Fourier Analysis](../algorithms/fourier.md).

```python
import pysicrs
```

## Transforms

### `fft`

```python
fft(real: list[float], imag: list[float]) -> (real: list[float], imag: list[float])
```

Complex DFT (mixed-radix FFT). Returns positive-frequency-packed amplitudes:
$X_k = \sum_n x_n e^{-2\pi i nk/N}$.

### `ifft`

```python
ifft(real: list[float], imag: list[float]) -> (real: list[float], imag: list[float])
```

Inverse DFT with $1/N$ normalization.

```python
from pysicrs import fft, ifft

x = [1.0, 0.0, 0.0, 0.0]
r, i = fft(x, [0.0]*4)
print(r, i)            # [1,1,1,1], [0,0,0,0]
xr, xi = ifft(r, i)
print(xr)              # [1,0,0,0]
```

## Spectral

### `spectral_derivative`

```python
spectral_derivative(f: list[float], dx: float, order: int) -> list[float]
```

Spectral derivative via $(ik)^m$ multipliers. Exponentially accurate for smooth periodic
data.

### `welch_psd`

```python
welch_psd(signal: list[float], nperseg: int, noverlap: int)
    -> (freqs: list[float], psd: list[float])
```

Welch's averaged periodogram (segment + window + average). `fs` is inferred from the array
index scaling; pass your sampling rate by scaling inputs accordingly.

```python
from pysicrs import welch_psd
import math

fs, n = 1000.0, 10000
x = [math.sin(2*math.pi*60.0*k/fs) for k in range(n)]
freqs, psd = welch_psd(x, nperseg=512, noverlap=256)
peak = max(range(len(freqs)), key=lambda i: psd[i])
print(f"dominant ≈ {freqs[peak]:.1f} relative units")
```

## Rust-only (not yet bound)

| function | purpose |
|----------|---------|
| `power_spectral_density(real, imag)` | raw periodogram |
| `fft_frequencies(n, dx)` | frequency axis |
| `fft_real` / `ifft_real` | real-input helpers |