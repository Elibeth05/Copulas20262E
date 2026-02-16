# Estimation of Distribution Algorithm (EDA) in Julia
# Author: Arturo Erdely
# Last update: 2026-01-01

"""
LA SIGUIENTE ES UNA FUNCION QU EMINIMIZA
* SI QUIERES maximizar, 
    Ssolo lo multiplicas por -1
* Si quieres que LA FUNCIÓN ALCANCE cierto valor- SE LO RESTAS ESE VALOR, LO PONES ELEVADO AL CUADRADO, LO MINIMIZAS Y ES LO MISMO
    EDA(fobj, valmin, valmax; iEnteros = zeros(Int, 0), tamgen = 1000, propselec = 0.3, difmax = 0.00001, maxiter = 1000)

`fobj` A real function of several variables to be minimized, where its argument is a vector or 1-D array.

`valmin, valmax` vectors or 1-D arrays of minimum and maximum values for the first generation.
DE AUI PARA ADELANTE SON PARÁMETROS OPCIONALES, PORQUE 
`iEnteros` index of variables that must take integer values., ESPECIFICAS, CUALES VAN A SER ENTEROS, SI NO le pasas el vector, asume que todos  van a ser números de punto flotante

`tamgen` size of each generation.

`propselec` proportion of population to be selected.

`difmax` error tolerance.

`maxiter` maximum number of iterations.

# Example 1
```
f(x) = (x[1] - 5)^4 - 16(x[1] - 5)^2 + 5(x[1] - 5) + 120
//los corchetes [] al lado de una variable indican que estás accediendo a un índice o elemento específico dentro de una estructura de datos,
//como un vector, una lista o un arreglo.
EDA(f, [0], [9])//que funcion y que intervalos son el mínimo y el máximo 
```

# Example 2
This non-negative function is clearly minimized at (5,-2).
```
//función de dos variables
f(z) = abs(z[1] - 5) + abs(z[2] + 2)
//el punto en el que se minimiza la fn no es deribable, ni continua
//pero ni continua necesito que sea
EDA(f, [-10, -10], [10, 10])
```

# Example 3
The same function but only allowing integer values:
```
EDA(f, [-10, -10], [10, 10], iEnteros = [1, 2])
```
"""
function EDA(fobj, valmin, valmax; iEnteros = zeros(Int, 0), tamgen = 1000,
             propselec = 0.3, difmax = 0.00001, maxiter = 1000)
    numiter = 1
    println("Iterating... ")
    numvar = length(valmin)
    nselec = Int(round(tamgen * propselec))
    G = zeros(tamgen, numvar)
    Gselec = zeros(nselec, numvar)
    for j ∈ 1:numvar
        G[:, j] = valmin[j] .+ (valmax[j] - valmin[j]) .* rand(tamgen)
    end
    if length(iEnteros) > 0
        for j ∈ iEnteros
            G[:, j] = round.(G[:, j])
        end
    end
    d(x, y) = sqrt(sum((x .- y) .^ 2))
    rnorm(n, μ, σ) = μ .+ (σ .* randn(n))
    promedio(x) = sum(x) / length(x)
    desvest(x) = sqrt(sum((x .- promedio(x)) .^ 2) / (length(x) - 1))
    fG = zeros(tamgen)
    maxGselec = zeros(tamgen)
    minGselec = zeros(tamgen)
    media = zeros(numvar)
    desv = zeros(numvar)
    while numiter < maxiter
        # evaluate objective function in current generation
        print(numiter, "\r")
        for i ∈ 1:tamgen
            fG[i] = fobj(G[i, :])
        end
        # selecting from current generation:
        umbral = sort(fG)[nselec]
        iSelec = findall(fG .≤ umbral)
        Gselec = G[iSelec, :]
        for j ∈ 1:numvar
            maxGselec[j] = maximum(Gselec[:, j])
            minGselec[j] = minimum(Gselec[:, j])
            media[j] = promedio(Gselec[:, j])
            desv[j] = desvest(Gselec[:, j])
        end
        # break cicle if stoping criteria is reached:
        if d(minGselec, maxGselec) < difmax 
            break
        end
        # otherwise draw a new generation:
        numiter += 1
        for j ∈ 1:numvar
            G[:, j] = rnorm(tamgen, media[j], desv[j])
        end
        if length(iEnteros) > 0
            for j ∈ iEnteros
                G[:, j] = round.(G[:, j])
            end
        end
    end
    println("...done")
    fGselec = zeros(nselec)
    for i ∈ eachindex(fGselec)
        fGselec[i] = fobj(Gselec[i, :])
    end
    xopt = Gselec[findmin(fGselec)[2], :]
    if length(iEnteros) > 0
        for j ∈ iEnteros
            xopt[j] = round(xopt[j])
        end
    end
    fxopt = fobj(xopt)
    r = (x = xopt, fx = fxopt, iter = numiter)
    if numiter == maxiter
        @warn "Maximum number of $maxiter iterations before reaching stopping criteria."
    end
    return r
end
