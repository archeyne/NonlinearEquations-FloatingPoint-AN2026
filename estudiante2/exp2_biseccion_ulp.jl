# Experimento 2 — Bisección cuando el punto medio ya no se distingue de un extremo
#
# Entre dos flotantes consecutivos a < b no existe ningún otro flotante, así que
# el punto medio calculado fl((a+b)/2) es a o b. Desde ese momento el intervalo
# no puede reducirse más: si la tolerancia pedida es menor que b - a (= 1 ulp),
# el criterio |b - a| ≤ tol_x NUNCA se cumple.
#
# En [1, 2) la separación entre flotantes es 2^(1-p), con p = bits de mantisa
# (11, 24, 53). Empezando con ancho 1, tras k pasos el ancho es 2^(-k), así que se
# predice el estancamiento en la iteración k = p (ancho 2^(1-p) = 1 ulp).
#
# Ejecutar:  julia --project=. exp2_biseccion_ulp.jl

include(joinpath(@__DIR__, "src", "comun.jl"))

f_raiz2(x) = x^2 - oftype(x, 2)
r2 = sqrt(setprecision(() -> big(2), BigFloat, 256))

con_reporte("exp2_biseccion_ulp") do io
    println(io, "# Experimento 2: bisección y la separación entre flotantes\n")

    println(io, "## 2.0 Parámetros de los formatos (IEEE 754)\n")
    tabla(io, ["tipo", "bits mantisa p", "eps = 2^(1-p)", "u = eps/2", "floatmin", "menor subnormal", "floatmax"],
          [(e.tipo, e.bits_mantisa, fmt(e.eps), fmt(e.u), fmt(e.floatmin), fmt(e.subnormal_min), fmt(e.floatmax))
           for e in eps_info.(TIPOS)])

    # ---------------------------------------------------------------------
    println(io, "## 2A. f(x) = x² − 2 en [1, 2], tol_x = 0, deteniendo al detectar m == a o m == b\n")
    filas = []
    for T in TIPOS
        R = biseccion(f_raiz2, T(1), T(2); tol_x = zero(T), maxit = 500, parar_cero = false)
        h = R.hist[end]
        consecutivos = nextfloat(h.a) == h.b
        push!(filas, (T, precision(T), R.estado, R.iter, string(h.a), string(h.b), consecutivos,
                      fmt(h.ancho), fmt(eps(h.a)), fmt(err(R.x, r2))))
        escribir_csv(joinpath(DIR_RES, "exp2_hist_$(T).csv"), ["k", "m", "ancho", "error", "residuo"],
                     [(g.k, g.m, g.ancho, err(g.m, r2), abs(g.fm)) for g in R.hist])
    end
    tabla(io, ["tipo", "p", "estado", "iter", "a final", "b final", "¿a,b consecutivos?", "b - a", "ulp(a)", "|m - √2|"], filas)

    # ---------------------------------------------------------------------
    println(io, "## 2B. Tolerancia absoluta fija tol_x = 1e-10 SIN detectar estancamiento (maxit = 200)\n")
    println(io, "Si tol_x < ulp cerca de la raíz, el ciclo solo termina por el tope de iteraciones. ",
                "Se compara con un criterio relativo b − a ≤ 2·eps(T)·|m|, que sí es alcanzable.\n")
    filas = []
    casos = (("x² − 2 en [1,2]", f_raiz2, 1.0, 2.0, r2),
             ("x² − 2e12 en [1e6,2e6]", x -> x^2 - oftype(x, 2e12), 1e6, 2e6, sqrt(big(2e12))))
    for (nombre, f, a, b, r) in casos, T in TIPOS
        (T == Float16 && a > 6e4) && continue          # 1e6 no es representable en Float16 (floatmax ≈ 65504)
        Ra = biseccion(f, T(a), T(b); tol_x = T(1e-10), maxit = 200, parar_estancado = false)
        Rr = biseccion(f, T(a), T(b); tol_x = zero(T), maxit = 200, parar_estancado = false)
        # criterio relativo: recorrer la historia sin estancamiento y ver cuándo se cumpliría
        krel = findfirst(h -> h.ancho <= 2 * eps(T) * abs(h.m), Rr.hist)
        push!(filas, (nombre, T, fmt(eps(T(Float64(r)))), Ra.estado, Ra.iter, fmt(err(Ra.x, r)),
                      krel === nothing ? "—" : string(krel)))
    end
    tabla(io, ["problema", "tipo", "ulp(r)", "estado (tol abs 1e-10)", "iter", "|m − r|", "iter criterio relativo"], filas)

    # ---------------------------------------------------------------------
    println(io, "## 2C. Fórmula del punto medio: (a+b)/2 frente a a+(b−a)/2 cerca de floatmax\n")
    filas = []
    for (T, a, b, c) in ((Float16, 6.0e4, 6.5e4, 6.2e4), (Float64, 1.0e308, 1.7e308, 1.5e308))
        f = x -> x - T(c)
        for medio in (:suma, :resta)
            R = biseccion(f, T(a), T(b); medio = medio, maxit = 300)
            push!(filas, (T, "[$(T(a)), $(T(b))]", medio === :suma ? "(a+b)/2" : "a+(b-a)/2",
                          R.estado, R.iter, string(R.x)))
        end
    end
    tabla(io, ["tipo", "intervalo", "punto medio", "estado", "iter", "x final"], filas)

    # ---------------------------------------------------------------------
    println(io, "## 2D. Prueba de signo: f(a)·f(m) < 0 frente a sign(f(a)) ≠ sign(f(m))\n")
    println(io, "Si |f(a)| y |f(m)| son pequeños, el producto se desborda por abajo (underflow) a 0, ",
                "la prueba `< 0` falla y el algoritmo descarta el subintervalo que contiene la raíz.\n")
    filas = []
    for (T, s) in ((Float16, 1e-3), (Float64, 1e-160))
        g = x -> T(s) * (x^2 - T(2))
        for signo in (:producto, :sign)
            R = biseccion(g, T(1), T(2); signo = signo, maxit = 300)
            push!(filas, (T, "$(s)·(x²−2)", signo === :producto ? "f(a)·f(m) < 0" : "sign",
                          R.estado, R.iter, string(R.x), fmt(err(R.x, r2))))
        end
    end
    tabla(io, ["tipo", "f", "prueba de signo", "estado", "iter", "x final", "|x − √2|"], filas)

    # ---------------------------------------------------------------------
    println(io, "## 2E. ¿Cuántos flotantes x cerca de √2 dan fl(x² − 2) == 0 exactamente?
")
    println(io, "La función CALCULADA puede anularse en varios puntos distintos de la raíz real: ",
                "un residuo exactamente cero no certifica error cero.
")
    filas = []
    for T in TIPOS
        x = T(sqrt(2.0))
        for _ in 1:2000; x = prevfloat(x); end        # 2000 flotantes a la izquierda
        ceros = T[]
        for _ in 1:4001
            f_raiz2(x) == 0 && push!(ceros, x)
            x = nextfloat(x)
        end
        maxerr = isempty(ceros) ? "—" : fmt(maximum(err(c, r2) for c in ceros))
        push!(filas, (T, length(ceros), isempty(ceros) ? "—" : string(first(ceros)),
                      isempty(ceros) ? "—" : string(last(ceros)), maxerr))
    end
    tabla(io, ["tipo", "# flotantes con fl(x²−2)=0", "menor", "mayor", "máx |x − √2|"], filas)
end
