### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 9193e94c-46df-4534-89f0-6a5699c2725a
begin
	using PlutoUI, PlutoTeachingTools
	
	# PlutoTeachingTools looks up language based on ENV["LANG"]
	# Uncomment a line below to override default language
	#set_language!(PlutoTeachingTools.EnglishUS())      # default
	#set_language!(PlutoTeachingTools.GermanGermany())  
end

# ╔═╡ 08a9d538-1eba-11ed-1ac4-c1b15c5bb280
using DataFrames

# ╔═╡ a28594c7-c134-4548-9ac1-1049bbe14f4b
using CSV

# ╔═╡ 903341d9-ea5e-4b9a-ac7e-f0b1a7c3b98b
using PyPlot

# ╔═╡ 54ee7295-27c4-4fcc-8c4f-27407b0c9ddc
TableOfContents()   # from PlutoUI

# ╔═╡ 4d39b023-b2e7-4569-9dde-9d1aab571c21
md"""
## Introduction
"""

# ╔═╡ ceb4fa2c-3981-466d-b0f8-0468bb264c68
md"""

This notebook is built with Pluto Markdown[^2] using $\LaTeX$[^3]

It looks to build the math around this video on [Pi hiding in prime regularities](https://www.youtube.com/watch?v=NaL_Cb42WyY&t=209s)

*Primitive* Pythagorean Triples (or **PPT** for short)[^1]

A *primitive* Pythagorean triple is as follows

$a^2 + b^2 = c^2$

"""

# ╔═╡ 5b5ff634-c6bf-4fe9-b442-8e9821156e75
md"""

In additive number theory, Fermat's theorem on sums of two squares states that an odd prime p can be expressed as


${\displaystyle p=x^{2}+y^{2},}$https://rosettacode.org/wiki/Modular_exponentiation#Julia
with x and y integers, if and only if

$p \equiv 1 \pmod{4}.$https://rosettacode.org/wiki/Modular_exponentiation#Julia
The prime numbers for which this is true are called *Pythagorean primes*. For example, the primes $5, 13, 17, 29, 37$ and $41$ are all congruent to 1 modulo 4, and they can be expressed as sums of two squares in the following ways:

$5 = 1^2 + 2^2, \quad 13 = 2^2 + 3^2, \quad 17 = 1^2 + 4^2, \quad 29 = 2^2 + 5^2, \quad 37 = 1^2 + 6^2, \quad 41 = 4^2 + 5^2.$

On the other hand, the primes $3, 7, 11, 19, 23$ and $31$ are all congruent to 3 modulo 4, and none of them can be expressed as the sum of two squares. This is the easier part of the theorem and follows immediately from the observation that all squares are congruent to 0 or 1 modulo 4.

See here: 
- https://www.wikiwand.com/en/Proofs_of_Fermat%27s_theorem_on_sums_of_two_squares
- https://www.had2know.org/academics/gaussian-prime-factorization-calculator.html
- https://stackoverflow.com/questions/2269810/whats-a-nice-method-to-factor-gaussian-integers
- https://en.wikipedia.org/wiki/Table_of_Gaussian_integer_factorizations
- [Modular exponentiation](https://rosettacode.org/wiki/Modular_exponentiation#Julia)

"""

# ╔═╡ 5ec4c34e-3c37-4e33-ad32-f920f00341cf
protip(md"- [Fermat's theorem on sums of two squares](https://www.wikiwand.com/en/ProofsofFermat%27stheoremonsumsoftwosquares)","Invitation to learn more")

# ╔═╡ ea560410-5c83-4fe2-afd5-03d09af0685d
function primeFactors(number, list = Int[] )
    for n in 2:number
		if (number % n == 0)
			return primeFactors(number / n, push!(list, n))
		end
	end
	list
end

# ╔═╡ 7463bf07-e164-4580-9a0f-581625994762
md"""

Let's test our function

"""

# ╔═╡ 72460280-bbfa-492d-8270-e44293266e09
df_pi = DataFrame(sqrt_radius= 1:1000) 

# ╔═╡ 2792710d-b98b-462f-bf03-120d75b7eee6
aside(tip(md"Extra information to consider.") )

# ╔═╡ 3c888d25-fc21-4c4f-a2bf-341ed6aaf787
set_aside_width(400)

# ╔═╡ bb47b3f0-b4d7-43d8-945f-ee94e9def8f6
function gcd(z1::Complex{T}, z2::Complex{V}) where {T<:Integer,V<:Integer}
    R = promote_type(T, V)
    a = Complex{R}(z1)
    b = Complex{R}(z2)
    if abs(a) < abs(b)
        a, b = b, a
    end
    while b != 0
        # TODO: Create rem(::Complex{<:Integer}, ::Complex{<:Integer})
        # TODO: Create div(::Complex{<:Integer}, ::Complex{<:Integer})
        b̅ = conj(b)
        # TODO: Handle overflow when calculating a*b̅
        t = a * b̅
        # TODO: Create checked_abs2(::Complex{<:Integer})
        # TODO: Handle overflow when calculating the norm of complex numbers
        abs2_b = abs2(b)
        abs2_b < 0 && __throw_gcd_overflow(z1, z2)
        r = a - b * complex(div(real(t), abs2_b, RoundNearest),
            div(imag(t), abs2_b, RoundNearest))
        a = b
        b = r
    end
    ar, ai = reim(a)
    if ar == 0
        complex(abs(ai), zero(ar))
    elseif ai == 0
        complex(abs(ar), zero(ai))
    elseif ar > 0 && ai >= 0   # gcd is already in first quadrant
        a
    elseif ar < 0 && ai > 0     # In second quadrant
        complex(ai, -ar)
    elseif ar < 0 && ai < 0     # In third quadrant
        -a
    else                               # In fourth quadrant
        complex(-ai, ar)
    end
end

# ╔═╡ e4af41d5-754a-466f-a26c-cb2a6c4e160c
dict = Dict()

# ╔═╡ 72bcf653-3ec8-43ad-90a1-e08f7b063b33
int(x) = floor(Int, x)

# ╔═╡ af516c3f-6fcf-4db7-a89f-30b1e813863c
function gprime(arr)
    x = Complex[]
    for p in arr
        if p == 2
            push!(x, 1 + 1im)
            push!(x, 1 - 1im)
        elseif mod(p, 4) == 3 # p = 3 mod 4, q = p.
            push!(x, p)
        elseif mod(p, 4) == 1 # p = 1 mod 4
            if haskey(dict, p)
                push!(x, dict[p])
                push!(x, conj(dict[p]))
            else
                for k in 2:(p-1)
                   # if mod(mod(BigInt(k)^((p - 1) / 2), p), p) == mod(-1, p)
					y = powermod(k,int((p-1)/2),p)
					if mod(y, p) == mod(-1, p)
                        #real = BigInt(k)^((p - 1) / 4)
						real = powermod(k,int((p - 1) / 4),p)
                        factor = gcd(Complex(p), (real) + 1im)
                        if !in(factor, x)
                            factor_complex_cong = conj(factor)
                            push!(x, factor)
                            push!(x, factor_complex_cong)
                            dict[p] = factor
                            break
                        end
                    end
                end
            end
        end
    end
    x
end

# ╔═╡ 83fee348-2ee5-4612-b909-f73a800b9cb6
function mod_4(x)
	if mod1(x,4) == 1
		1
	elseif mod1(x, 4) == 3
		3
	end
end

# ╔═╡ 2de72be6-3ebd-4282-8ab0-d51bbe727aea
# https://oodlescoop.com/tutorials/julia/programs/julia-program-to-check-if-a-number-is-prime-number-or-not-;jsessionid=2F142D01C14896A498C338A5BAE587FF
function prime(num)    
    isPrime = true
        for i in 2:convert(Int64, round(num/2, digits = 0))
            if num % i == 0
                isPrime = false
                break
            end
        end
    return isPrime
end

# ╔═╡ 9d7518a1-95a5-42fc-b268-87b100c3e96d
df_pi.isPrime = prime.(df_pi.sqrt_radius)

# ╔═╡ c056d31e-6351-4d76-b484-fb8013c20b67
df_pi.mod4 = mod_4.(df_pi.sqrt_radius)

# ╔═╡ 2c8e6f69-f0b9-43a1-b18e-94a4e5de993d
df_pi.factors = primeFactors.(df_pi.sqrt_radius)

# ╔═╡ b899d763-5636-4142-b93b-30b578144d01
df_pi.cc = gprime.(df_pi.factors)

# ╔═╡ e2c0fbd0-1cdd-4d60-a12c-ba78d90748a1
# https://stackoverflow.com/questions/67698311/how-to-get-product-of-all-elements-in-a-row-of-matrix-in-julia
function cartesian(array)
	x = Complex[]
	array_size = size(array)
	len = length(array)
	#@show(len)
	if len == 1
		push!(x, array[1])
	elseif(isodd(len))
		x
	else
		# TODO this needs work
		len2 = Int(len/2)
		P = reshape(array, 2, len2)	
		c1 = prod.(eachrow(P))
		push!(x, c1[1])
		push!(x, c1[2])
		Q = reshape(array, len2, 2)
		c2 = prod.(eachcol(Q))
		push!(x, c2[1])
		push!(x, c2[2])
	end
	unique(x, dims=1)
end

# ╔═╡ 593279c7-1757-4c08-a19d-a3d269423040
df_pi.cartesian = cartesian.(df_pi.cc)

# ╔═╡ 7ec2f302-bd04-4e0e-97be-56e93ecc361c
df_pi

# ╔═╡ dd397413-3822-4efc-816d-2fe4ad5e1629
function complex_filter(mod4, prime_factors)::Bool # https://juliadatascience.io/filter_subset
	interesting_mod4 = 3 == mod4
	#interesting_prime_factors = 1 == size(prime_factors)
	interesting_prime_factors = length(prime_factors) == 1
	interesting_mod4 && interesting_prime_factors
end

# ╔═╡ cfe930f7-bcd4-44c3-8b69-c25e0cb6ef13
filter([:mod4, :factors] => !complex_filter, df_pi)

# ╔═╡ 9ca7013a-3115-4931-b9da-034934e50295
function img_ops(collection)
	unique(vcat(collection, collection * (0+1im), collection * (0-1im), collection * (-1+0im)), dims=1)
end

# ╔═╡ 09f4b07a-877c-4357-bbe1-4d46606e8dec
df_pi.img_ops = img_ops.(df_pi.cartesian)

# ╔═╡ e3cbf4ec-8e7b-44ce-887f-b978524c5cbf
df_pi.size = size.(df_pi.img_ops)

# ╔═╡ af40960b-1d79-42ed-a69c-2f24465cfc02
df_pi

# ╔═╡ a9d4a058-7ab0-42eb-b879-3787dcedf1c0
# https://discourse.julialang.org/t/how-to-convert-all-nothings-in-dataframe-to-missing/54004/10
df_pi.mod4 = replace(df_pi.mod4, nothing => missing)

# ╔═╡ c193d0f7-44d5-44ab-9a0e-1271e23b06f6

CSV.write("prime_pi.csv", df_pi)


# ╔═╡ 3858bddb-34c8-476a-9d27-b2997b54271b
twenty_five = [-3+4im,
0+5im,
4+3im,
3+4im,
5+0im,
3-4im]

# ╔═╡ 1e010dfd-e31a-426d-99bc-09cb5f43ed6c
twenty_five * 1im

# ╔═╡ 3bde1246-2ecc-4167-9629-75a1a3e5b221
nums = ComplexF64.([1,2,4],[2,2,-1])

# ╔═╡ 1a7ef4a2-50e0-4a72-804e-d07bfb6d2fdc
# ╠═╡ disabled = true
#=╠═╡
polar.(Base.vect.(0.0,angle.(nums)),Base.vect.(0.0,abs.(nums)),marker="o")
  ╠═╡ =#

# ╔═╡ c999e688-6387-4b3d-951e-3ab87b230d06
d = [-3+4im,
0+5im,
4+3im,
3+4im,
5+0im,
3-4im]

# ╔═╡ 0f279169-a831-406d-ae6c-975ef16848a7
# ╠═╡ disabled = true
#=╠═╡
plot(d)
  ╠═╡ =#

# ╔═╡ b0cc8956-47de-4f5f-a55f-4893fe84cb49
# ╠═╡ disabled = true
#=╠═╡
plot(real(d),imag(d)) # or directly with plot(d)
  ╠═╡ =#

# ╔═╡ 68913bd8-92f9-41d5-ad90-5e9b7ac75678
# ╠═╡ disabled = true
#=╠═╡
plot.show()
  ╠═╡ =#

# ╔═╡ bab47f30-e938-4403-8c4d-0677fcd25433
function cartesian_old(array)
	x = Complex[]
	array_size = size(array)
	if array_size[1] == 1
		push!(x, array[1])
	else
		for m in array, n in array # Corteston product
			push!(x, m .* n)
		end
	end
	unique(x, dims=1)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.3.4"
PlutoTeachingTools = "~0.2.3"
PlutoUI = "~0.7.40"
PyPlot = "~2.11.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0"
manifest_format = "2.0"
project_hash = "dedbf61e0507577c6a41f1761824a61c6547c88d"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "1833bda4a027f4b2a1c984baddcf755d77266818"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.1.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "78bee250c6826e1cf805a88b7f1e86025275d208"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.46.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "316daa94fad0b7a008ebd573e002efd6609d85ac"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.19"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "d19f9edd8c34760dca2de2b503f969d8700ed288"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "0f960b1404abb0b244c1ece579a0ec78d056a5d1"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.15"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "0e8bcc235ec8367a8e9648d48325ff00e4b0a545"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.5"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "d8be3432505c2febcea02f44e5f4396fae017503"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "a602d7b0babfca89005da04d89223b867b55319f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.40"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "53b8b07b721b77144a0fbbbc2675222ebf40a02d"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.94.1"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "f9d953684d4d21e947cb6d642db18853d43cb027"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "dad726963ecea2d8a81e26286f625aee09a91b7c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.4.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─9193e94c-46df-4534-89f0-6a5699c2725a
# ╟─54ee7295-27c4-4fcc-8c4f-27407b0c9ddc
# ╟─4d39b023-b2e7-4569-9dde-9d1aab571c21
# ╟─ceb4fa2c-3981-466d-b0f8-0468bb264c68
# ╠═5b5ff634-c6bf-4fe9-b442-8e9821156e75
# ╠═5ec4c34e-3c37-4e33-ad32-f920f00341cf
# ╠═08a9d538-1eba-11ed-1ac4-c1b15c5bb280
# ╠═ea560410-5c83-4fe2-afd5-03d09af0685d
# ╠═7463bf07-e164-4580-9a0f-581625994762
# ╠═72460280-bbfa-492d-8270-e44293266e09
# ╠═2792710d-b98b-462f-bf03-120d75b7eee6
# ╠═3c888d25-fc21-4c4f-a2bf-341ed6aaf787
# ╠═bb47b3f0-b4d7-43d8-945f-ee94e9def8f6
# ╠═e4af41d5-754a-466f-a26c-cb2a6c4e160c
# ╠═72bcf653-3ec8-43ad-90a1-e08f7b063b33
# ╠═af516c3f-6fcf-4db7-a89f-30b1e813863c
# ╠═83fee348-2ee5-4612-b909-f73a800b9cb6
# ╠═2de72be6-3ebd-4282-8ab0-d51bbe727aea
# ╠═9d7518a1-95a5-42fc-b268-87b100c3e96d
# ╠═c056d31e-6351-4d76-b484-fb8013c20b67
# ╠═2c8e6f69-f0b9-43a1-b18e-94a4e5de993d
# ╠═b899d763-5636-4142-b93b-30b578144d01
# ╠═e2c0fbd0-1cdd-4d60-a12c-ba78d90748a1
# ╠═593279c7-1757-4c08-a19d-a3d269423040
# ╠═7ec2f302-bd04-4e0e-97be-56e93ecc361c
# ╠═dd397413-3822-4efc-816d-2fe4ad5e1629
# ╠═cfe930f7-bcd4-44c3-8b69-c25e0cb6ef13
# ╠═9ca7013a-3115-4931-b9da-034934e50295
# ╠═09f4b07a-877c-4357-bbe1-4d46606e8dec
# ╠═e3cbf4ec-8e7b-44ce-887f-b978524c5cbf
# ╠═af40960b-1d79-42ed-a69c-2f24465cfc02
# ╠═a9d4a058-7ab0-42eb-b879-3787dcedf1c0
# ╠═a28594c7-c134-4548-9ac1-1049bbe14f4b
# ╠═c193d0f7-44d5-44ab-9a0e-1271e23b06f6
# ╠═3858bddb-34c8-476a-9d27-b2997b54271b
# ╠═1e010dfd-e31a-426d-99bc-09cb5f43ed6c
# ╠═903341d9-ea5e-4b9a-ac7e-f0b1a7c3b98b
# ╠═3bde1246-2ecc-4167-9629-75a1a3e5b221
# ╠═1a7ef4a2-50e0-4a72-804e-d07bfb6d2fdc
# ╠═c999e688-6387-4b3d-951e-3ab87b230d06
# ╠═0f279169-a831-406d-ae6c-975ef16848a7
# ╠═b0cc8956-47de-4f5f-a55f-4893fe84cb49
# ╠═68913bd8-92f9-41d5-ad90-5e9b7ac75678
# ╠═bab47f30-e938-4403-8c4d-0677fcd25433
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
