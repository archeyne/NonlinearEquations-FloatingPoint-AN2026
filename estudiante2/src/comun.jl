# comun.jl — utilidades compartidas por los scripts de experimentos
include(joinpath(@__DIR__, "metodos.jl"))
using .Metodos
using Printf

const TIPOS = (Float16, Float32, Float64)
const DIR_RES = joinpath(@__DIR__, "..", "resultados")
mkpath(DIR_RES)

# Error hacia adelante |x - r| calculado en BigFloat (x se convierte EXACTAMENTE
# a BigFloat; r es la referencia de 256 bits) y reportado como Float64.
err(x, r) = Float64(abs(big(x) - r))

fmt(x) = @sprintf("%.3e", Float64(x))

# Tabla en formato Markdown en la consola (y opcionalmente en un archivo)
function tabla(io, encabezado, filas)
    println(io, "| ", join(encabezado, " | "), " |")
    println(io, "|", join(fill("---", length(encabezado)), "|"), "|")
    for f in filas
        println(io, "| ", join(string.(f), " | "), " |")
    end
    println(io)
end

# Escribe a la vez en consola y en resultados/<nombre>.md
function con_reporte(f, nombre)
    ruta = joinpath(DIR_RES, nombre * ".md")
    buf = IOBuffer()
    f(buf)
    s = String(take!(buf))
    print(s)
    write(ruta, s)
    println("→ guardado en ", normpath(ruta))
end
