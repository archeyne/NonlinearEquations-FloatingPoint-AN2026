# Estudiante 2 — Residuo, error y parada

Pregunta: ¿por qué un residuo pequeño |f(xₙ)| no garantiza que xₙ esté cerca de una raíz?,
¿qué pasa en bisección cuando el punto medio ya no se distingue de un extremo?, y ¿cómo se
comparan bisección y punto fijo x = F(x) en cuanto a criterios de parada y precisión finita?

## Requisitos

- Julia ≥ 1.10 (probado con 1.12.6).
- Los experimentos solo usan la biblioteca estándar. Los gráficos usan CairoMakie
  (declarado en `Project.toml`).

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. exp1_residuo_error.jl
julia --project=. exp2_biseccion_ulp.jl
julia --project=. exp3_punto_fijo.jl
julia --project=. graficos.jl
```

## Archivos

| archivo | contenido |
|---|---|
| `src/metodos.jl` | bisección, punto fijo y Newton genéricos en el tipo `T` (Float16/32/64/BigFloat), con registro de la historia y del criterio de parada activado |
| `src/comun.jl` | error en BigFloat, tablas y escritura de resultados |
| `exp1_residuo_error.jl` | residuo vs error: atan(x−1), (x−1)⁷ factorizado y desarrollado, x·e⁻ˣ, escalamiento de f |
| `exp2_biseccion_ulp.jl` | estancamiento de bisección en 1 ulp, tolerancias inalcanzables, desbordamiento de (a+b)/2, underflow en f(a)·f(m) |
| `exp3_punto_fijo.jl` | x = cos x y contracción lenta x − λ(x²−2): cota L/(1−L), ciclos y puntos fijos flotantes |
| `graficos.jl` | figuras en `figuras/` (PDF y PNG) |
| `resultados/` | tablas `.md` y datos `.csv` de cada experimento |
| `informe/estudiante2.tex` | sección del informe |

## Precisión usada

Todas las operaciones de cada método se hacen en el tipo de los datos de entrada
(las constantes se convierten con `T(...)`). Julia implementa la aritmética de
`Float16` en software, calculando en Float32 y redondeando a Float16 en cada operación.
Para +, −, ×, ÷ y √ este doble redondeo da el mismo resultado que un redondeo directo,
porque 24 ≥ 2·11 + 2 (Figueroa 1995, "When is double rounding innocuous?"). Las funciones
elementales (`cos`, `exp`, `atan`) en Float16 también se evalúan en Float32 y se redondean,
así que no son necesariamente correctamente redondeadas; esto se declara en el informe.
Las raíces de referencia se calculan con Newton en `BigFloat` de 256 bits, y los errores
|x − r| se evalúan convirtiendo x exactamente a BigFloat.
