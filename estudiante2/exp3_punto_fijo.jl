# Experimento 3 — Punto fijo x = F(x) frente a bisección
#
# Teoría (teorema de punto fijo / contracción): si |F'(x)| ≤ L < 1 cerca de r,
#     |x_n - r| ≤ L/(1-L) · |x_n - x_{n-1}|       (cota a posteriori)
# Con L cerca de 1 el factor L/(1-L) es grande: un PASO pequeño no implica un
# ERROR pequeño. Bisección, en cambio, garantiza |m - r| ≤ (b - a)/2 siempre.
#
# En aritmética finita: x_{n+1} = fl(F(x_n)). Si |F(x) - x| < ulp(x)/2 entonces
# fl(F(x)) == x: aparece un "punto fijo flotante" que NO es el punto fijo real.
# Para F(x) = x - λ g(x), eso pasa cuando λ|g(x)| < ulp(x)/2, aunque x esté lejos de r.
#
# Ejecutar:  julia --project=. exp3_punto_fijo.jl

include(joinpath(@__DIR__, "src", "comun.jl"))

# --- Problema A: x = cos x  (L = |sin r| ≈ 0.674) -------------------------------
F_cos(x) = cos(x)
g_cos(x) = x - cos(x)                                     # misma raíz, para bisección
r_cos = raiz_referencia(g_cos, x -> 1 + sin(x), 0.7)
L_cos = Float64(sin(r_cos))

# --- Problema B: contracción lenta F(x) = x - λ(x² - 2), raíz √2 ----------------
λ = 1e-3
F_lento(x) = x - oftype(x, λ) * (x^2 - oftype(x, 2))
r2 = sqrt(setprecision(() -> big(2), BigFloat, 256))
L_lento = Float64(1 - 2 * big(λ) * r2)                    # F'(√2) = 1 - 2λ√2 ≈ 0.99717

con_reporte("exp3_punto_fijo") do io
    println(io, "# Experimento 3: punto fijo frente a bisección\n")
    println(io, "Referencias: r_cos = ", Float64(r_cos), " (Newton en BigFloat de 256 bits); ",
                "√2 en BigFloat de 256 bits.\n")

    # ---------------------------------------------------------------------
    println(io, "## 3A. x = cos x, x0 = 1, L = |sin r| ≈ ", round(L_cos, digits = 4),
                ", L/(1−L) ≈ ", round(L_cos / (1 - L_cos), digits = 3), "\n")
    filas = []
    for T in TIPOS, tol in (1e-3, 1e-6, 1e-10)
        tolT = T(tol)
        P = punto_fijo(F_cos, T(1); tol_x = tolT, maxit = 5000)
        B = biseccion(g_cos, T(0), T(1); tol_x = 2tolT, maxit = 500)   # garantiza |m-r| ≤ tol
        cota = L_cos / (1 - L_cos) * Float64(P.hist[end].paso)
        push!(filas, (T, tol, P.estado, P.iter, fmt(err(P.x, r_cos)), fmt(cota),
                      B.estado, B.iter, fmt(err(B.x, r_cos))))
    end
    tabla(io, ["tipo", "tol", "PF estado", "PF iter", "PF |x−r|", "cota L/(1−L)·paso",
               "Bis estado", "Bis iter", "Bis |m−r|"], filas)

    for T in TIPOS
        P = punto_fijo(F_cos, T(1); tol_x = zero(T), maxit = 5000)
        escribir_csv(joinpath(DIR_RES, "exp3_hist_cos_$(T).csv"), ["k", "x", "paso", "error"],
                     [(h.k, h.x, h.paso, err(h.x, r_cos)) for h in P.hist])
        B = biseccion(g_cos, T(0), T(1); tol_x = zero(T), maxit = 500)
        escribir_csv(joinpath(DIR_RES, "exp3_hist_bis_cos_$(T).csv"), ["k", "m", "ancho", "error"],
                     [(h.k, h.m, h.ancho, err(h.m, r_cos)) for h in B.hist])
    end

    println(io, "### Iterando hasta que la aritmética no permite avanzar (tol_x = 0)\n")
    filas = []
    for T in TIPOS
        P = punto_fijo(F_cos, T(1); tol_x = zero(T), maxit = 5000)
        ciclo = P.estado === :ciclo ? join(string.(T[h.x for h in P.hist[end-P.periodo+1:end]]), " ↔ ") : "—"
        push!(filas, (T, P.estado, P.iter, P.periodo, ciclo, fmt(err(P.x, r_cos)), fmt(eps(P.x))))
    end
    tabla(io, ["tipo", "estado", "iter", "periodo", "valores del ciclo", "|x − r|", "ulp(x)"], filas)

    # ---------------------------------------------------------------------
    println(io, "## 3B. Contracción lenta F(x) = x − λ(x² − 2), λ = $λ, x0 = 1: L = F'(√2) ≈ ",
                round(L_lento, digits = 5), ", L/(1−L) ≈ ", round(L_lento / (1 - L_lento), digits = 1), "\n")
    println(io, "El criterio |x_n − x_{n−1}| ≤ tol se cumple con un error real ≈ 350 veces mayor. ",
                "Además en Float16/Float32 aparece un punto fijo flotante (x_{n+1} == x_n) lejos de √2.\n")
    filas = []
    for T in TIPOS, tol in (1e-4, 1e-6, 1e-8, 1e-10)
        tolT = T(tol)
        P = punto_fijo(F_lento, T(1); tol_x = tolT, maxit = 200_000)
        paso = Float64(P.hist[end].paso)
        push!(filas, (T, tol, P.estado, P.iter, fmt(paso), fmt(err(P.x, r2)),
                      fmt(L_lento / (1 - L_lento) * paso)))
    end
    tabla(io, ["tipo", "tol", "estado", "iter", "último paso", "|x − √2| real", "cota L/(1−L)·paso"], filas)

    # Predicción del estancamiento: fl(F(x)) == x cuando λ|x²−2| < ulp(x)/2,
    # es decir |x − √2| ≲ ulp(√2) / (4 λ √2)
    println(io, "### Predicción del error de estancamiento: |x − √2| ≈ ulp(√2)/(4λ√2)\n")
    filas = []
    for T in TIPOS
        P = punto_fijo(F_lento, T(1); tol_x = zero(T), maxit = 200_000)
        pred = Float64(eps(T(sqrt(2.0)))) / (4λ * sqrt(2.0))
        push!(filas, (T, P.estado, P.iter, string(P.x), fmt(err(P.x, r2)), fmt(pred)))
        escribir_csv(joinpath(DIR_RES, "exp3_hist_lento_$(T).csv"), ["k", "x", "paso", "error"],
                     [(h.k, h.x, h.paso, err(h.x, r2)) for h in P.hist])
    end
    tabla(io, ["tipo", "estado", "iter", "x final", "|x − √2| observado", "predicción"], filas)
end
