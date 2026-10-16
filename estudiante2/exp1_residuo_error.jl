# Experimento 1 — Residuo pequeño ≠ error pequeño
#
# Idea teórica (Taylor alrededor de una raíz simple r):
#     f(x) = f'(r)(x - r) + O((x - r)^2)   ⇒   |x - r| ≈ |f(x)| / |f'(r)|
# El factor 1/|f'(r)| es el NÚMERO DE CONDICIÓN de la raíz. Si f'(r) es pequeña
# (o nula: raíz múltiple), un residuo diminuto es compatible con un error grande.
#
# Ejecutar:  julia --project=. exp1_residuo_error.jl

include(joinpath(@__DIR__, "src", "comun.jl"))

# --- Funciones de prueba (genéricas en el tipo T de x) -----------------------
f_atan(x)     = atan(x - one(x))               # raíz simple, f'(1) = 1  (bien condicionada)
f7_fact(x)    = (x - one(x))^7                 # raíz de multiplicidad 7, forma factorizada
# (x-1)^7 desarrollado y evaluado con Horner: x^7 -7x^6 +21x^5 -35x^4 +35x^3 -21x^2 +7x -1
function f7_exp(x::T) where {T}
    c = T.((1, -7, 21, -35, 35, -21, 7, -1))
    p = c[1]
    for k in 2:8
        p = p * x + c[k]
    end
    p
end
f_xexp(x)  = x * exp(-x)                       # única raíz r = 0; f(x) → 0 si x → ∞
df_xexp(x) = (one(x) - x) * exp(-x)

con_reporte("exp1_residuo_error") do io
    println(io, "# Experimento 1: residuo vs. error\n")

    # ---------------------------------------------------------------------
    println(io, "## 1A. Residuo |f(1+δ)| frente al error δ (Float64)\n")
    println(io, "Para atan(x-1) el residuo ≈ error (f'(1)=1). Para (x-1)^7 el residuo es δ^7: ",
                "un residuo de 1e-14 corresponde a un error de 1e-2.\n")
    filas = []
    for p in 1:8
        δ = 10.0^(-p)
        x = 1.0 + δ
        push!(filas, (fmt(δ), fmt(abs(f_atan(x))), fmt(abs(f7_fact(x))), fmt(abs(f7_exp(x)))))
    end
    tabla(io, ["δ = x-1", "|atan(x-1)|", "|(x-1)^7| factorizado", "|(x-1)^7| desarrollado"], filas)

    # Barrido fino para la figura del "ruido" del polinomio desarrollado
    xs = range(0.98, 1.02; length = 4001)
    escribir_csv(joinpath(DIR_RES, "exp1_barrido_f7.csv"),
                 ["x", "f7_fact_f64", "f7_exp_f64", "f7_exp_f32"],
                 [(x, f7_fact(x), f7_exp(x), f7_exp(Float32(x))) for x in xs])

    cambios(v) = (w = filter(!iszero, v); count(i -> sign(w[i]) != sign(w[i+1]), 1:length(w)-1))
    println(io, "Cambios de signo de f calculada en 4001 puntos de [0.98, 1.02] (la función real cambia de signo UNA vez):
")
    tabla(io, ["evaluación", "cambios de signo", "máx |f| observado"],
          [("factorizada Float64", cambios(f7_fact.(xs)), fmt(maximum(abs, f7_fact.(xs)))),
           ("desarrollada Float64", cambios(f7_exp.(xs)), fmt(maximum(abs, f7_exp.(xs)))),
           ("desarrollada Float32", cambios(f7_exp.(Float32.(xs))), fmt(maximum(abs, f7_exp.(Float32.(xs)))))])

    # ---------------------------------------------------------------------
    println(io, "## 1B. Bisección parando por residuo |f(m)| ≤ ε_f en [0.5, 1.7], raíz r = 1\n")
    println(io, "Estado `tol_f`: el algoritmo declara éxito porque el residuo es pequeño. ",
                "Estado `estancado`: el punto medio ya coincide con un extremo.\n")
    filas = []
    for T in TIPOS, (nombre, f) in (("atan(x-1)", f_atan), ("(x-1)^7 fact.", f7_fact), ("(x-1)^7 desarr.", f7_exp))
        tolf = T(1e-12) == 0 ? nextfloat(zero(T)) : T(1e-12)   # en Float16, 1e-12 se redondea a 0
        R = biseccion(f, T(0.5), T(1.7); tol_f = tolf, maxit = 200)
        push!(filas, (T, nombre, fmt(tolf), R.estado, R.iter, string(R.x), fmt(abs(f(R.x))), fmt(err(R.x, big(1)))))
    end
    tabla(io, ["tipo", "f", "ε_f usado", "estado", "iter", "x final", "|f(x)|", "|x - 1|"], filas)

    # Historia completa en Float64 para graficar error y residuo por iteración
    for (nombre, f) in (("atan", f_atan), ("f7fact", f7_fact), ("f7exp", f7_exp))
        R = biseccion(f, 0.5, 1.7; tol_f = 0.0, maxit = 200)
        escribir_csv(joinpath(DIR_RES, "exp1_hist_biseccion_$(nombre).csv"),
                     ["k", "m", "residuo", "error", "ancho"],
                     [(h.k, h.m, abs(h.fm), err(h.m, big(1)), h.ancho) for h in R.hist])
    end

    # ---------------------------------------------------------------------
    println(io, "## 1C. Newton en f(x) = x e^{-x} desde x0 = 2, parando SOLO por residuo\n")
    println(io, "La única raíz es r = 0, pero f(x) → 0 cuando x → ∞: la iteración se aleja ",
                "y el residuo cae por debajo de la tolerancia. Falsa convergencia.\n")
    filas = []
    for T in TIPOS
        R = newton(f_xexp, df_xexp, T(2); tol_f = T(1e-6), maxit = 100)
        push!(filas, (T, R.estado, R.iter, string(R.x), fmt(abs(f_xexp(R.x))), fmt(err(R.x, big(0)))))
    end
    tabla(io, ["tipo", "estado", "iter", "x final", "|f(x)|", "|x - r|"], filas)
    R = newton(f_xexp, df_xexp, 2.0; tol_f = 1e-6, maxit = 100)
    escribir_csv(joinpath(DIR_RES, "exp1_newton_xexp.csv"), ["k", "x", "residuo", "error"],
                 [(h.k, h.x, abs(h.fx), abs(h.x)) for h in R.hist])

    # ---------------------------------------------------------------------
    println(io, "## 1D. Escalamiento: la misma raíz √2 con f(x) = s·(x² − 2) (Float64, ε_f = 1e-8)\n")
    println(io, "Multiplicar f por s no cambia la raíz pero multiplica el residuo por s. ",
                "Un criterio de residuo absoluto depende de la escala de f: con s pequeño para ",
                "casi enseguida (error grande); con s grande la tolerancia es inalcanzable.\n")
    r2 = sqrt(big(2))
    filas = []
    for s in (1e-10, 1e-4, 1.0, 1e8)
        g(x) = s * (x^2 - 2)
        R = biseccion(g, 1.0, 2.0; tol_f = 1e-8, maxit = 200)
        κ = 1 / (s * 2 * sqrt(2))   # número de condición 1/|f'(r)|
        push!(filas, (s, fmt(κ), R.estado, R.iter, fmt(abs(g(R.x))), fmt(err(R.x, r2))))
    end
    tabla(io, ["s", "κ = 1/|f'(r)|", "estado", "iter", "|f(x)|", "|x - √2|"], filas)
end
