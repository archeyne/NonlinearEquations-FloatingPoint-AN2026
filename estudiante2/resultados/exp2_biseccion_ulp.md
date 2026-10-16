# Experimento 2: bisección y la separación entre flotantes

## 2.0 Parámetros de los formatos (IEEE 754)

| tipo | bits mantisa p | eps = 2^(1-p) | u = eps/2 | floatmin | menor subnormal | floatmax |
|---|---|---|---|---|---|---|
| Float16 | 11 | 9.766e-04 | 4.883e-04 | 6.104e-05 | 5.960e-08 | 6.550e+04 |
| Float32 | 24 | 1.192e-07 | 5.960e-08 | 1.175e-38 | 1.401e-45 | 3.403e+38 |
| Float64 | 53 | 2.220e-16 | 1.110e-16 | 2.225e-308 | 4.941e-324 | 1.798e+308 |

## 2A. f(x) = x² − 2 en [1, 2], tol_x = 0, deteniendo al detectar m == a o m == b

| tipo | p | estado | iter | a final | b final | ¿a,b consecutivos? | b - a | ulp(a) | |m - √2| |
|---|---|---|---|---|---|---|---|---|---|
| Float16 | 11 | estancado | 11 | 1.413 | 1.414 | true | 9.766e-04 | 9.766e-04 | 1.511e-04 |
| Float32 | 24 | estancado | 24 | 1.4142135 | 1.4142137 | true | 1.192e-07 | 1.192e-07 | 9.501e-08 |
| Float64 | 53 | estancado | 53 | 1.414213562373095 | 1.4142135623730951 | true | 2.220e-16 | 2.220e-16 | 1.254e-16 |

## 2B. Tolerancia absoluta fija tol_x = 1e-10 SIN detectar estancamiento (maxit = 200)

Si tol_x < ulp cerca de la raíz, el ciclo solo termina por el tope de iteraciones. Se compara con un criterio relativo b − a ≤ 2·eps(T)·|m|, que sí es alcanzable.

| problema | tipo | ulp(r) | estado (tol abs 1e-10) | iter | |m − r| | iter criterio relativo |
|---|---|---|---|---|---|---|
| x² − 2 en [1,2] | Float16 | 9.766e-04 | residuo_cero | 7 | 1.511e-04 | — |
| x² − 2 en [1,2] | Float32 | 1.192e-07 | max_iter | 200 | 9.501e-08 | 23 |
| x² − 2 en [1,2] | Float64 | 2.220e-16 | tol_x | 34 | 4.705e-11 | 52 |
| x² − 2e12 en [1e6,2e6] | Float32 | 1.250e-01 | max_iter | 200 | 6.237e-02 | 23 |
| x² − 2e12 en [1e6,2e6] | Float64 | 2.328e-10 | max_iter | 200 | 1.432e-10 | 52 |

## 2C. Fórmula del punto medio: (a+b)/2 frente a a+(b−a)/2 cerca de floatmax

| tipo | intervalo | punto medio | estado | iter | x final |
|---|---|---|---|---|---|
| Float16 | [6.0e4, 6.5e4] | (a+b)/2 | no_finito | 1 | Inf |
| Float16 | [6.0e4, 6.5e4] | a+(b-a)/2 | residuo_cero | 7 | 6.2e4 |
| Float64 | [1.0e308, 1.7e308] | (a+b)/2 | no_finito | 1 | Inf |
| Float64 | [1.0e308, 1.7e308] | a+(b-a)/2 | residuo_cero | 51 | 1.5e308 |

## 2D. Prueba de signo: f(a)·f(m) < 0 frente a sign(f(a)) ≠ sign(f(m))

Si |f(a)| y |f(m)| son pequeños, el producto se desborda por abajo (underflow) a 0, la prueba `< 0` falla y el algoritmo descarta el subintervalo que contiene la raíz.

| tipo | f | prueba de signo | estado | iter | x final | |x − √2| |
|---|---|---|---|---|---|---|
| Float16 | 0.001·(x²−2) | f(a)·f(m) < 0 | estancado | 11 | 1.5 | 8.579e-02 |
| Float16 | 0.001·(x²−2) | sign | residuo_cero | 7 | 1.414 | 1.511e-04 |
| Float64 | 1.0e-160·(x²−2) | f(a)·f(m) < 0 | estancado | 53 | 1.421875 | 7.661e-03 |
| Float64 | 1.0e-160·(x²−2) | sign | estancado | 53 | 1.414213562373095 | 1.254e-16 |

## 2E. ¿Cuántos flotantes x cerca de √2 dan fl(x² − 2) == 0 exactamente?

La función CALCULADA puede anularse en varios puntos distintos de la raíz real: un residuo exactamente cero no certifica error cero.

| tipo | # flotantes con fl(x²−2)=0 | menor | mayor | máx |x − √2| |
|---|---|---|---|---|
| Float16 | 1 | 1.414 | 1.414 | 1.511e-04 |
| Float32 | 0 | — | — | — |
| Float64 | 0 | — | — | — |

