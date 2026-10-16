# Experimento 3: punto fijo frente a bisección

Referencias: r_cos = 0.7390851332151607 (Newton en BigFloat de 256 bits); √2 en BigFloat de 256 bits.

## 3A. x = cos x, x0 = 1, L = |sin r| ≈ 0.6736, L/(1−L) ≈ 2.064

| tipo | tol | PF estado | PF iter | PF |x−r| | cota L/(1−L)·paso | Bis estado | Bis iter | Bis |m−r| |
|---|---|---|---|---|---|---|---|---|
| Float16 | 0.001 | tol_x | 15 | 3.156e-04 | 2.015e-03 | tol_x | 9 | 1.149e-03 |
| Float16 | 1.0e-6 | ciclo | 17 | 3.156e-04 | 1.008e-03 | estancado | 12 | 1.727e-04 |
| Float16 | 1.0e-10 | ciclo | 17 | 3.156e-04 | 1.008e-03 | estancado | 12 | 1.727e-04 |
| Float32 | 0.001 | tol_x | 17 | 3.248e-04 | 1.665e-03 | tol_x | 9 | 1.149e-03 |
| Float32 | 1.0e-6 | tol_x | 35 | 2.934e-07 | 1.476e-06 | tol_x | 19 | 8.894e-07 |
| Float32 | 1.0e-10 | ciclo | 40 | 6.423e-08 | 2.460e-07 | residuo_cero | 24 | 4.629e-09 |
| Float64 | 0.001 | tol_x | 17 | 3.248e-04 | 1.665e-03 | tol_x | 9 | 1.149e-03 |
| Float64 | 1.0e-6 | tol_x | 34 | 3.931e-07 | 2.016e-06 | tol_x | 19 | 8.894e-07 |
| Float64 | 1.0e-10 | tol_x | 58 | 2.995e-11 | 1.536e-10 | tol_x | 33 | 8.873e-11 |

### Iterando hasta que la aritmética no permite avanzar (tol_x = 0)

| tipo | estado | iter | periodo | valores del ciclo | |x − r| | ulp(x) |
|---|---|---|---|---|---|---|
| Float16 | ciclo | 17 | 2 | 0.7393 ↔ 0.739 | 3.156e-04 | 4.883e-04 |
| Float32 | ciclo | 40 | 2 | 0.7390851 ↔ 0.7390852 | 6.423e-08 | 5.960e-08 |
| Float64 | punto_fijo_flotante | 91 | 0 | — | 3.064e-17 | 1.110e-16 |

## 3B. Contracción lenta F(x) = x − λ(x² − 2), λ = 0.001, x0 = 1: L = F'(√2) ≈ 0.99717, L/(1−L) ≈ 352.6

El criterio |x_n − x_{n−1}| ≤ tol se cumple con un error real ≈ 350 veces mayor. Además en Float16/Float32 aparece un punto fijo flotante (x_{n+1} == x_n) lejos de √2.

| tipo | tol | estado | iter | último paso | |x − √2| real | cota L/(1−L)·paso |
|---|---|---|---|---|---|---|
| Float16 | 0.0001 | punto_fijo_flotante | 237 | 0.000e+00 | 1.837e-01 | 0.000e+00 |
| Float16 | 1.0e-6 | punto_fijo_flotante | 237 | 0.000e+00 | 1.837e-01 | 0.000e+00 |
| Float16 | 1.0e-8 | punto_fijo_flotante | 237 | 0.000e+00 | 1.837e-01 | 0.000e+00 |
| Float16 | 1.0e-10 | punto_fijo_flotante | 237 | 0.000e+00 | 1.837e-01 | 0.000e+00 |
| Float32 | 0.0001 | tol_x | 918 | 9.966e-05 | 3.560e-02 | 3.514e-02 |
| Float32 | 1.0e-6 | tol_x | 2547 | 9.537e-07 | 3.571e-04 | 3.362e-04 |
| Float32 | 1.0e-8 | punto_fijo_flotante | 3507 | 0.000e+00 | 2.101e-05 | 0.000e+00 |
| Float32 | 1.0e-10 | punto_fijo_flotante | 3507 | 0.000e+00 | 2.101e-05 | 0.000e+00 |
| Float64 | 0.0001 | tol_x | 917 | 9.997e-05 | 3.570e-02 | 3.525e-02 |
| Float64 | 1.0e-6 | tol_x | 2552 | 9.991e-07 | 3.523e-04 | 3.522e-04 |
| Float64 | 1.0e-8 | tol_x | 4178 | 9.990e-09 | 3.522e-06 | 3.522e-06 |
| Float64 | 1.0e-10 | tol_x | 5804 | 9.986e-11 | 3.521e-08 | 3.521e-08 |

### Predicción del error de estancamiento: |x − √2| ≈ ulp(√2)/(4λ√2)

| tipo | estado | iter | x final | |x − √2| observado | predicción |
|---|---|---|---|---|---|
| Float16 | punto_fijo_flotante | 237 | 1.23 | 1.837e-01 | 1.726e-01 |
| Float32 | punto_fijo_flotante | 3507 | 1.4141926 | 2.101e-05 | 2.107e-05 |
| Float64 | punto_fijo_flotante | 10604 | 1.4142135623730558 | 3.921e-14 | 3.925e-14 |

