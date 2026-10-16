# Experimento 1: residuo vs. error

## 1A. Residuo |f(1+δ)| frente al error δ (Float64)

Para atan(x-1) el residuo ≈ error (f'(1)=1). Para (x-1)^7 el residuo es δ^7: un residuo de 1e-14 corresponde a un error de 1e-2.

| δ = x-1 | |atan(x-1)| | |(x-1)^7| factorizado | |(x-1)^7| desarrollado |
|---|---|---|---|
| 1.000e-01 | 9.967e-02 | 1.000e-07 | 1.000e-07 |
| 1.000e-02 | 1.000e-02 | 1.000e-14 | 7.994e-15 |
| 1.000e-03 | 1.000e-03 | 1.000e-21 | 1.998e-15 |
| 1.000e-04 | 1.000e-04 | 1.000e-28 | 1.776e-15 |
| 1.000e-05 | 1.000e-05 | 1.000e-35 | 3.331e-16 |
| 1.000e-06 | 1.000e-06 | 1.000e-42 | 1.110e-16 |
| 1.000e-07 | 1.000e-07 | 1.000e-49 | 4.441e-16 |
| 1.000e-08 | 1.000e-08 | 1.000e-56 | 1.887e-15 |

Cambios de signo de f calculada en 4001 puntos de [0.98, 1.02] (la función real cambia de signo UNA vez):

| evaluación | cambios de signo | máx |f| observado |
|---|---|---|
| factorizada Float64 | 1 | 1.280e-12 |
| desarrollada Float64 | 717 | 1.282e-12 |
| desarrollada Float32 | 1955 | 3.338e-06 |

## 1B. Bisección parando por residuo |f(m)| ≤ ε_f en [0.5, 1.7], raíz r = 1

Estado `tol_f`: el algoritmo declara éxito porque el residuo es pequeño. Estado `estancado`: el punto medio ya coincide con un extremo.

| tipo | f | ε_f usado | estado | iter | x final | |f(x)| | |x - 1| |
|---|---|---|---|---|---|---|---|
| Float16 | atan(x-1) | 5.960e-08 | residuo_cero | 10 | 1.0 | 0.000e+00 | 0.000e+00 |
| Float16 | (x-1)^7 fact. | 5.960e-08 | residuo_cero | 3 | 0.9497 | 0.000e+00 | 5.029e-02 |
| Float16 | (x-1)^7 desarr. | 5.960e-08 | estancado | 12 | 0.6543 | 3.418e-03 | 3.457e-01 |
| Float32 | atan(x-1) | 1.000e-12 | residuo_cero | 23 | 1.0 | 0.000e+00 | 0.000e+00 |
| Float32 | (x-1)^7 fact. | 1.000e-12 | tol_f | 5 | 0.9875001 | 4.768e-14 | 1.250e-02 |
| Float32 | (x-1)^7 desarr. | 1.000e-12 | estancado | 26 | 0.9161725 | 5.960e-08 | 8.383e-02 |
| Float64 | atan(x-1) | 1.000e-12 | tol_f | 39 | 0.9999999999992725 | 7.275e-13 | 7.275e-13 |
| Float64 | (x-1)^7 fact. | 1.000e-12 | tol_f | 5 | 0.9875 | 4.768e-14 | 1.250e-02 |
| Float64 | (x-1)^7 desarr. | 1.000e-12 | tol_f | 5 | 0.9875 | 4.985e-14 | 1.250e-02 |

## 1C. Newton en f(x) = x e^{-x} desde x0 = 2, parando SOLO por residuo

La única raíz es r = 0, pero f(x) → 0 cuando x → ∞: la iteración se aleja y el residuo cae por debajo de la tolerancia. Falsa convergencia.

| tipo | estado | iter | x final | |f(x)| | |x - r| |
|---|---|---|---|---|---|
| Float16 | tol_f | 12 | 16.53 | 1.013e-06 | 1.653e+01 |
| Float32 | tol_f | 13 | 17.60673 | 3.973e-07 | 1.761e+01 |
| Float64 | tol_f | 13 | 17.606730006234756 | 3.973e-07 | 1.761e+01 |

## 1D. Escalamiento: la misma raíz √2 con f(x) = s·(x² − 2) (Float64, ε_f = 1e-8)

Multiplicar f por s no cambia la raíz pero multiplica el residuo por s. Un criterio de residuo absoluto depende de la escala de f: con s pequeño para casi enseguida (error grande); con s grande la tolerancia es inalcanzable.

| s | κ = 1/|f'(r)| | estado | iter | |f(x)| | |x - √2| |
|---|---|---|---|---|---|
| 1.0e-10 | 3.536e+09 | tol_f | 1 | 2.500e-11 | 8.579e-02 |
| 0.0001 | 3.536e+03 | tol_f | 13 | 8.200e-09 | 2.899e-05 |
| 1.0 | 3.536e-01 | tol_f | 27 | 5.237e-09 | 1.851e-09 |
| 1.0e8 | 3.536e-09 | estancado | 53 | 4.441e-08 | 1.254e-16 |

