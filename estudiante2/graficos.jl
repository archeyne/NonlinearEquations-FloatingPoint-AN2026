# graficos.jl — genera las figuras del informe a partir de resultados/*.csv
# Requiere haber corrido antes exp1, exp2 y exp3.
# Ejecutar:  julia --project=. graficos.jl

using CairoMakie

const RES = joinpath(@__DIR__, "resultados")
const FIG = joinpath(@__DIR__, "figuras")
mkpath(FIG)

# Lector CSV mínimo: devuelve Dict(columna => Vector{Float64})
function leer(nombre)
    lineas = readlines(joinpath(RES, nombre))
    cols = split(lineas[1], ",")
    datos = [parse.(Float64, split(l, ",")) for l in lineas[2:end]]
    Dict(c => [d[i] for d in datos] for (i, c) in enumerate(cols))
end
# En escala log los ceros no se pueden dibujar: se reemplazan por NaN
pos(v) = [x > 0 ? x : NaN for x in v]

const COL = Dict("Float16" => :firebrick, "Float32" => :darkorange, "Float64" => :royalblue)
guardar(nombre, fig) = (save(joinpath(FIG, nombre * ".pdf"), fig); save(joinpath(FIG, nombre * ".png"), fig, px_per_unit = 2))

# --- Figura 1: ruido de (x-1)^7 desarrollado -----------------------------------
let d = leer("exp1_barrido_f7.csv")
    fig = Figure(size = (900, 360))
    ax1 = Axis(fig[1, 1], title = "(x−1)⁷ en Float64", xlabel = "x", ylabel = "f(x) calculada")
    lines!(ax1, d["x"], d["f7_exp_f64"], color = (:gray, 0.8), label = "desarrollado (Horner)")
    lines!(ax1, d["x"], d["f7_fact_f64"], color = :royalblue, linewidth = 2, label = "factorizado")
    ylims!(ax1, -1e-14, 1e-14)                     # acercamiento para ver el ruido ~1e-15
    axislegend(ax1, position = :lt)
    ax2 = Axis(fig[1, 2], title = "(x−1)⁷ desarrollado en Float32", xlabel = "x")
    lines!(ax2, d["x"], d["f7_exp_f32"], color = :darkorange)
    guardar("fig1_ruido_f7", fig)
end

# --- Figura 2: residuo vs error en bisección (Float64) -------------------------
let fig = Figure(size = (900, 360))
    for (j, (arch, tit)) in enumerate((("atan", "atan(x−1)"), ("f7fact", "(x−1)⁷ factorizado"), ("f7exp", "(x−1)⁷ desarrollado")))
        d = leer("exp1_hist_biseccion_$(arch).csv")
        ax = Axis(fig[1, j], title = tit, xlabel = "iteración k", yscale = log10,
                  ylabel = j == 1 ? "magnitud" : "")
        scatterlines!(ax, d["k"], pos(d["error"]), label = "error |m−1|", markersize = 5)
        scatterlines!(ax, d["k"], pos(d["residuo"]), label = "residuo |f(m)|", markersize = 5)
        j == 1 && axislegend(ax, position = :lb)
    end
    guardar("fig2_residuo_vs_error", fig)
end

# --- Figura 3: ancho del intervalo en bisección y estancamiento ---------------
let fig = Figure(size = (700, 400))
    ax = Axis(fig[1, 1], xlabel = "iteración k", ylabel = "ancho b − a", yscale = log2,
              title = "Bisección de x² − 2 en [1,2]: el ancho se detiene en 1 ulp")
    # Float64 primero y los de menor precisión encima (los anchos coinciden hasta estancarse)
    for (T, p, w) in (("Float64", 53, 1.5), ("Float32", 24, 4), ("Float16", 11, 7))
        d = leer("exp2_hist_$(T).csv")
        scatterlines!(ax, d["k"], d["ancho"], color = COL[T], label = "$T (p = $p)",
                      markersize = 2w, linewidth = w)
        hlines!(ax, [2.0^(1 - p)], color = COL[T], linestyle = :dash)
    end
    axislegend(ax, position = :lb)
    guardar("fig3_biseccion_ulp", fig)
end

# --- Figura 4: x = cos x, punto fijo vs bisección -----------------------------
let fig = Figure(size = (700, 400))
    ax = Axis(fig[1, 1], xlabel = "iteración k", ylabel = "error |x_k − r|", yscale = log10,
              title = "x = cos x: punto fijo (línea continua) y bisección (discontinua)")
    for T in ("Float16", "Float32", "Float64")
        d = leer("exp3_hist_cos_$(T).csv")
        lines!(ax, d["k"], pos(d["error"]), color = COL[T], label = T)
        b = leer("exp3_hist_bis_cos_$(T).csv")
        lines!(ax, b["k"], pos(b["error"]), color = COL[T], linestyle = :dash)
    end
    axislegend(ax, position = :rt)
    guardar("fig4_cos_pf_vs_biseccion", fig)
end

# --- Figura 5: contracción lenta, paso vs error -------------------------------
let fig = Figure(size = (900, 380))
    ax = Axis(fig[1, 1], xlabel = "iteración k", yscale = log10, xscale = log10,
              title = "F(x) = x − 10⁻³(x² − 2): paso (discontinua) y error (continua)")
    for (T, w) in (("Float64", 1.5), ("Float32", 4), ("Float16", 4))
        d = leer("exp3_hist_lento_$(T).csv")
        lines!(ax, d["k"], pos(d["error"]), color = COL[T], label = T, linewidth = w)
        lines!(ax, d["k"], pos(d["paso"]), color = COL[T], linestyle = :dash, linewidth = w)
    end
    axislegend(ax, position = :lb)
    guardar("fig5_contraccion_lenta", fig)
end

# --- Figura 6: Newton en x e^{-x} --------------------------------------------
let d = leer("exp1_newton_xexp.csv")
    fig = Figure(size = (700, 360))
    ax = Axis(fig[1, 1], xlabel = "iteración k", yscale = log10,
              title = "Newton en x·e⁻ˣ desde x₀ = 2 (Float64): la raíz es r = 0")
    scatterlines!(ax, d["k"], d["error"], label = "error |x_k − 0|")
    scatterlines!(ax, d["k"], d["residuo"], label = "residuo |f(x_k)|")
    hlines!(ax, [1e-6], color = :gray, linestyle = :dash, label = "tolerancia ε_f = 10⁻⁶")
    axislegend(ax, position = :lb)
    guardar("fig6_newton_xexp", fig)
end

println("Figuras guardadas en ", normpath(FIG))
