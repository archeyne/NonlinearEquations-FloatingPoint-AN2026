# metodos.jl — Estudiante 2: residuo, error y parada
#
# Implementaciones propias (sin funciones de caja negra) de bisección, punto fijo
# y Newton. Todas son genéricas en el tipo de punto flotante T: si los datos de
# entrada son Float16, TODAS las operaciones internas (punto medio, resta,
# comparación) se hacen en Float16. Las constantes se convierten con T(...) para
# no promover accidentalmente a Float64.
#
# Cada método devuelve un NamedTuple con:
#   hist   : vector de filas (NamedTuple) con la historia de la iteración
#   estado : Symbol con el criterio de parada que se activó
#   x      : última aproximación
#   iter   : número de iteraciones realizadas

module Metodos

using Printf
export biseccion, punto_fijo, newton, raiz_referencia, escribir_csv, eps_info

"""
    biseccion(f, a, b; tol_x, tol_f, maxit, medio, signo, parar_estancado)

Bisección en el tipo T = typeof(a).

- `medio = :suma`  usa m = (a+b)/2   (puede desbordarse si a+b > floatmax(T))
  `medio = :resta` usa m = a + (b-a)/2
- `signo = :producto` decide el subintervalo con f(a)*f(m) < 0 (puede dar 0 por underflow)
  `signo = :sign`     decide con sign(f(a)) != sign(f(m))
- `parar_estancado = true` detiene el ciclo cuando m == a o m == b, es decir, cuando
  el punto medio ya no se distingue de un extremo en la precisión usada.
- `parar_cero = false` NO se detiene cuando fl(f(m)) == 0 (sign(0) = 0 se trata
  como cambio de signo y se toma b = m).

Estados posibles: :tol_x (|b-a| <= tol_x), :tol_f (|f(m)| <= tol_f), :residuo_cero,
:estancado (m coincide con un extremo), :no_finito (Inf/NaN), :max_iter.
"""
function biseccion(f, a::T, b::T; tol_x = zero(T), tol_f = zero(T), maxit = 200,
                   medio = :resta, signo = :sign, parar_estancado = true,
                   parar_cero = true) where {T<:AbstractFloat}
    a < b || throw(ArgumentError("se requiere a < b"))
    fa, fb = f(a), f(b)
    sign(fa) * sign(fb) < 0 || throw(ArgumentError("f(a) y f(b) deben tener signos opuestos"))
    dos = T(2)
    hist = NamedTuple[]
    estado = :max_iter
    m = a
    for k in 1:maxit
        m = medio === :suma ? (a + b) / dos : a + (b - a) / dos
        fm = f(m)
        push!(hist, (k = k, a = a, b = b, m = m, fm = fm, ancho = b - a))
        if !isfinite(m) || !isfinite(fm)
            estado = :no_finito; break
        end
        if parar_cero && fm == 0
            estado = :residuo_cero; break
        end
        if tol_f > 0 && abs(fm) <= tol_f
            estado = :tol_f; break
        end
        if parar_estancado && (m == a || m == b)
            estado = :estancado; break
        end
        cambio = signo === :producto ? (fa * fm < 0) : (sign(fa) != sign(fm))
        if cambio
            b, fb = m, fm
        else
            a, fa = m, fm
        end
        if b - a <= tol_x
            estado = :tol_x; break
        end
    end
    return (hist = hist, estado = estado, x = m, iter = length(hist))
end

"""
    punto_fijo(F, x0; tol_x, maxit, detectar_ciclo)

Iteración x_{n+1} = F(x_n) en el tipo T = typeof(x0). Se detiene cuando
|x_{n+1} - x_n| <= tol_x. Como el conjunto de flotantes es finito, la sucesión
calculada termina siendo periódica: con `detectar_ciclo = true` se detecta el
regreso a un valor ya visitado.

Estados: :tol_x, :punto_fijo_flotante (x_{n+1} == x_n exactamente),
:ciclo (regresó a un valor visitado sin cumplir tol_x; `periodo` = longitud del ciclo),
:no_finito, :max_iter.
"""
function punto_fijo(F, x0::T; tol_x = zero(T), maxit = 1000,
                    detectar_ciclo = true) where {T<:AbstractFloat}
    hist = NamedTuple[]
    vistos = Dict{T,Int}(x0 => 0)
    x = x0
    estado = :max_iter
    periodo = 0
    for k in 1:maxit
        xn = F(x)
        paso = abs(xn - x)
        push!(hist, (k = k, x = xn, paso = paso))
        if !isfinite(xn)
            estado = :no_finito; x = xn; break
        end
        if xn == x
            estado = :punto_fijo_flotante; x = xn; break
        end
        if paso <= tol_x
            estado = :tol_x; x = xn; break
        end
        if detectar_ciclo && haskey(vistos, xn)
            estado = :ciclo; periodo = k - vistos[xn]; x = xn; break
        end
        vistos[xn] = k
        x = xn
    end
    return (hist = hist, estado = estado, x = x, iter = length(hist), periodo = periodo)
end

"""
    newton(f, df, x0; tol_x, tol_f, maxit)

Newton x_{n+1} = x_n - f(x_n)/f'(x_n) en el tipo T. Se usa solo para ilustrar
el peligro de parar únicamente con el residuo.
"""
function newton(f, df, x0::T; tol_x = zero(T), tol_f = zero(T), maxit = 100) where {T<:AbstractFloat}
    hist = NamedTuple[]
    x = x0
    estado = :max_iter
    for k in 1:maxit
        fx, dfx = f(x), df(x)
        if dfx == 0
            estado = :derivada_cero; break
        end
        xn = x - fx / dfx
        fxn = f(xn)
        push!(hist, (k = k, x = xn, fx = fxn, paso = abs(xn - x)))
        x = xn
        if !isfinite(xn) || !isfinite(fxn)
            estado = :no_finito; break
        end
        if abs(fxn) <= tol_f
            estado = :tol_f; break
        end
        if abs(hist[end].paso) <= tol_x
            estado = :tol_x; break
        end
    end
    return (hist = hist, estado = estado, x = x, iter = length(hist))
end

"""
    raiz_referencia(f, df, x0; bits = 256)

Raíz de referencia calculada con Newton en BigFloat de `bits` bits de mantisa
(≈ 77 dígitos decimales para 256 bits). Se itera hasta que el paso es menor que
2^(-bits+10); el resultado es "exacto" para comparar con Float16/32/64, pero NO
es la raíz exacta: se documenta como referencia numérica.
"""
function raiz_referencia(f, df, x0; bits = 256)
    setprecision(BigFloat, bits) do
        x = BigFloat(x0)
        for _ in 1:200
            dx = f(x) / df(x)
            x -= dx
            abs(dx) < big(2.0)^(-bits + 10) * max(1, abs(x)) && break
        end
        x
    end
end

"Información básica del formato T: épsilon de máquina, redondeo unitario, mínimos y máximo."
eps_info(T) = (tipo = string(T), eps = Float64(eps(T)), u = Float64(eps(T)) / 2,
               floatmin = Float64(floatmin(T)), subnormal_min = Float64(nextfloat(zero(T))),
               floatmax = Float64(floatmax(T)), bits_mantisa = precision(T))

"Escribe un CSV sencillo (sin dependencias externas)."
function escribir_csv(ruta, encabezado::Vector{String}, filas)
    open(ruta, "w") do io
        println(io, join(encabezado, ","))
        for fila in filas
            println(io, join((string(v) for v in fila), ","))
        end
    end
    return ruta
end

end # module
