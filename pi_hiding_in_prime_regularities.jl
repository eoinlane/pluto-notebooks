### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ eae0d1c0-8b19-43e4-97f8-21f89708018e
begin
	using PlutoUI, PlutoTeachingTools, HypertextLiteral
	
	# PlutoTeachingTools looks up language based on ENV["LANG"]
	# Uncomment a line below to override default language
	#set_language!(PlutoTeachingTools.EnglishUS())      # default
	#set_language!(PlutoTeachingTools.GermanGermany())  
end

# ╔═╡ fe677804-cf1f-4413-abe4-684943d68d5f
using CairoMakie, Makie.GeometryBasics, LaTeXStrings

# ╔═╡ 08a9d538-1eba-11ed-1ac4-c1b15c5bb280
using DataFrames, StatsBase

# ╔═╡ a28594c7-c134-4548-9ac1-1049bbe14f4b
using CSV

# ╔═╡ 4bb7d51f-ef33-4228-a39f-0f39a709899a
TableOfContents()   # from PlutoUI

# ╔═╡ 4e5fccf3-d7d3-4e7f-be3f-4eb24f5a4239
@htl("""
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

<article class="learning">
	<h4>
		Learning HTML and CSS
	</h4>
	<p>
		The story here is to tie together &pi;, prime numbers, complex numbers and pi and come up with a formula for &pi; using an alternating infinite sum. When you see &pi; see up in math, there’s a circle hiding somewhere, sometimes very sneakily. So the goal here is to discover alternating infinite sum, and find the circle hiding behind it.
	</p>
	<p>
		An outline for the approach is as follows:
	</p>
<p>
<ol>
<li>Count lattics points </li>
<li>Things like 17 = 4<sup>2</sup> + 1<sup>2</sup> </li>
<li>Things like 17 = (4 + <i>i</i>)(4 - i<i>i</i>)  </li>
<li>Introduce  	&#967; </li>
</ol>
</p>
</article>


<style>

	article.learning {
		background: #f996a84f;
		padding: 1em;
		border-radius: 5px;
	}

	article.learning h4::before {
		content: "☝️";
	}

	article.learning p::first-letter {
		font-size: 1.5em;
		font-family: cursive;
	}

</style>
""")

# ╔═╡ ca796498-1218-4109-b90f-5038ab76c629
tip(md" A “lattice point” is a point (a, b) in the plane where a and b are both integers, a point where grid lines cross. .")

# ╔═╡ 4c72470d-729d-480d-ac59-f65f195ec4ce
md"""

$\frac{π}{4}= \frac{1}{1}+ \frac{0}{2}+\frac{-1}{3}+\frac{0}{4}+\frac{1}{5}+\ldots$

"""

# ╔═╡ 08df6171-cc09-4239-b458-2a8ba87d69ca
aside(tip(md"This notebook is built with [Pluto Markdown](https://www.juliafordatascience.com/first-steps-5-pluto/) using $\LaTeX$") )

# ╔═╡ 73b1b527-d820-4531-92a1-73df0a492b0b
md"""
## Introduction
"""

# ╔═╡ ceb4fa2c-3981-466d-b0f8-0468bb264c68
md"""
It looks to build the math around this video on [Pi hiding in prime regularities](https://www.youtube.com/watch?v=NaL_Cb42WyY&t=209s)

*Primitive* Pythagorean Triples (or **PPT** for short)[^1]

A *primitive* Pythagorean triple is as follows

$a^2 + b^2 = c^2$

"""

# ╔═╡ 5b5ff634-c6bf-4fe9-b442-8e9821156e75
md"""

In additive number theory, Fermat's theorem on sums of two squares states that an odd prime p can be expressed as


${\displaystyle p=x^{2}+y^{2},}$
with x and y integers, if and only if

$p \equiv 1 \pmod{4}.$
The prime numbers for which this is true are called *Pythagorean primes*. For example, the primes $5, 13, 17, 29, 37$ and $41$ are all congruent to 1 modulo 4, and they can be expressed as sums of two squares in the following ways:

$5 = 1^2 + 2^2, \quad 13 = 2^2 + 3^2, \quad 17 = 1^2 + 4^2, \quad 29 = 2^2 + 5^2, \quad 37 = 1^2 + 6^2, \quad 41 = 4^2 + 5^2.$

On the other hand, the primes $3, 7, 11, 19, 23$ and $31$ are all congruent to 3 modulo 4, and none of them can be expressed as the sum of two squares. This is the easier part of the theorem and follows immediately from the observation that all squares are congruent to 0 or 1 modulo 4.

See here: 
- https://www.wikiwand.com/en/Proofs_of_Fermat%27s_theorem_on_sums_of_two_squares
- https://www.had2know.org/academics/gaussian-prime-factorization-calculator.html
- https://stackoverflow.com/questions/2269810/whats-a-nice-method-to-factor-gaussian-integers
- https://en.wikipedia.org/wiki/Table_of_Gaussian_integer_factorizations

"""

# ╔═╡ f89f32a4-2d1d-4905-90d9-6746576bf432
# https://docs.makie.org/stable/tutorials/basic-tutorial/
begin
	f = Figure()
	ax = Axis(f[1, 1], 
    	title = "A Makie Axis",
    	xlabel = "The x label",
    	ylabel = "The y label"
	)
#x = range(0, 10, length=100)
#y = sin.(x)
	r = 5
	Θ = LinRange(0,2*π, 500)
	x = r*sin.(Θ)
	y = r*cos.(Θ)
	# https://discourse.julialang.org/t/how-to-add-grid-lines-on-top-of-a-heatmap-in-makie/77578
	f, ax, l1 = lines(x, y, linewidth = .5, color = :red, label = "cicle";
		figure = (; resolution = (500, 500)),	
		axis = (; title = L"\frac{\sin{x}}{x}", xlabel = L"\Re(z)", ylabel = L"\Im(z)", aspect = DataAspect(), xgridcolor = :black, ygridcolor = :black, xgridwidth = 0.5, ygridwidth = 0.5, xminorgridcolor = :grey,
    	yminorgridcolor = :grey,
		xminorgridvisible = true,
		yminorgridvisible = true,
		xminorticks = IntervalsBetween(3),
		yminorticks = IntervalsBetween(3),
		#https://github.com/MakieOrg/Makie.jl/issues/158
		backgroundcolor = :transparent,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
		)
	)
	scatter!(ax, 4, 3, color = :yellow, label = "point")
	scatter!(ax, 3, 4, color = :yellow)
	scatter!(ax, -3, 4, color = :yellow)
	scatter!(ax, -4, 3, color = :yellow)
	scatter!(ax, -3, -4, color = :yellow)
	scatter!(ax, -4, -3, color = :yellow)
	scatter!(ax, 3, -4, color = :yellow)
	scatter!(ax, 4, -3, color = :yellow)
	#\textcolor{blue}
	text!(L"(-3+4i)", position=(-3,4))
	text!(L"\sqrt{25}", position=(2,0))
	axislegend()
	xs = [Point2f(0,0)]
	ys = [Point2f(5,0)]
	arrows!(ax, xs, ys, linewidth = 2, arrowsize = 15, color = :red)
	#tightlimits!(ax)
f
end

# ╔═╡ a921b52a-9922-4346-92d6-c5cd03540381
# https://docs.makie.org/stable/tutorials/basic-tutorial/
begin
	f1 = Figure()
#x = range(0, 10, length=100)
#y = sin.(x)
	r1 = 5
	Θ1 = LinRange(0,2*π, 500)
	x1 = r1*sin.(Θ1)
	y1 = r1*cos.(Θ1)
	# https://discourse.julialang.org/t/how-to-add-grid-lines-on-top-of-a-heatmap-in-makie/77578
	lines(x1, y1, linewidth = .5, color = :red, label = "cicle";
		figure = (; resolution = (500, 500)),	
		axis = (; title = L"\frac{\sin{x}}{x}", xlabel = L"\Re(z)", ylabel = L"\Im(z)", aspect = DataAspect(), xgridcolor = :black, ygridcolor = :black, xgridwidth = 0.5, ygridwidth = 0.5, xminorgridcolor = :grey,
    	yminorgridcolor = :grey,
		xminorgridvisible = true,
		yminorgridvisible = true,
		xminorticks = IntervalsBetween(3),
		yminorticks = IntervalsBetween(3),
		#https://github.com/MakieOrg/Makie.jl/issues/158
		backgroundcolor = :transparent,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
		)
	)
	for i in 0:r1
		scatter!(4, r, color = :red)
	end
	
	#\textcolor{blue}
	text!(L"(-3+4i)", position=(-3,4))
	text!(L"\sqrt{25}", position=(2,0))
	axislegend()
	xs1 = [Point2f(0,0)]
	ys1 = [Point2f(5,0)]
	#arrows!(ax, xs, ys, linewidth = 2, arrowsize = 15, color = :red)
	#tightlimits!(ax)
f1
end

# ╔═╡ ea560410-5c83-4fe2-afd5-03d09af0685d
begin
	"""
    primeFactors(number, list = Int[] )

Compute the prime factor of a number

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function primeFactors(number, list = Int[] )
    for n in 2:number
		if (number % n == 0)
			return primeFactors(number / n, push!(list, n))
		end
	end
	list
end
end

# ╔═╡ 7463bf07-e164-4580-9a0f-581625994762
md"""

Let's test our function

"""

# ╔═╡ 72460280-bbfa-492d-8270-e44293266e09
df_pi = DataFrame(sqrt_radius= 1:100) #73 is the magic number

# ╔═╡ bb47b3f0-b4d7-43d8-945f-ee94e9def8f6
begin
	"""
    gcd(z1, z2)

Greatest common divisor operator. `x * y * z *...` calls this function with a pair of complex number arguments, i.e. `gcd(z1, z2)`.
"""
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
end

# ╔═╡ e4af41d5-754a-466f-a26c-cb2a6c4e160c
dict = Dict()

# ╔═╡ 72bcf653-3ec8-43ad-90a1-e08f7b063b33
int(x) = floor(Int, x)

# ╔═╡ af516c3f-6fcf-4db7-a89f-30b1e813863c
begin
	"""
    bar(x[, y])

Compute the Bar index between `x` and `y`.

If `y` is unspecified, compute the Bar index between all pairs of columns of `x`.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
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
						end

# ╔═╡ 83fee348-2ee5-4612-b909-f73a800b9cb6
begin
	"""
    mod_4(x)

For inputs 1 above a multiple of 4, the output is 1. For inputs 3 above a multiple of 4, it outputs 3. 
	
# Examples
```julia-repl
julia>  mod_4(7) = 3
```
""" 
function mod_4(x)
	if mod1(x,4) == 1
		1
	elseif mod1(x, 4) == 3
		3
	end
end
end

# ╔═╡ 763c51c1-a1a8-4262-8d40-14799603279a
begin
	"""
    Χ(n)

For inputs 1 above a multiple of 4, the output of Χ is 1. For inputs 3 above a multiple of 4, it outputs -1. And for even numbers, it gives 0. if you evaluate Χ on the natural numbers, it gives this nice cyclic pattern 1, 0, -1, 0, and repeat.
	
# Examples
```julia-repl
julia>  Χ(7) = -1
```
""" 
function Χ(n)
		if mod_4(n) == 1
			1
		elseif mod_4(n) == 3
			-1
		else iseven(n)
			0
		end
end

end

# ╔═╡ 46b2299a-3c88-4f5f-9f46-3c773c70b189
aside(tip(md"Extra information to consider. https://rosettacode.org/wiki/Modular_exponentiation#Julia2") )

# ╔═╡ 7bc5e477-37b0-435e-9224-152bc562ff04
function computeChi(array)
	D = Dict
	x = []
	collection = []
	prod = 1
		#@show(array)
		D = countmap(array)
		for (key, val) in D
			s = 0
			for i in 0:val
				s += Χ(key^i)
				#push!(collection,s)
			end
			prod = prod * s
		end
		#print(prod)
	4*prod
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

# ╔═╡ c3af66e7-7e7c-4af9-90cc-8ba7382b2791


# ╔═╡ bcd7bbc4-90c8-492c-903d-b9559d8a2eb5
df_pi

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

# ╔═╡ 10ed3abf-ae21-453c-ae7b-c27d0f106dae
df_pi.chi = computeChi.(df_pi.factors)

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

# ╔═╡ Cell order:
# ╟─eae0d1c0-8b19-43e4-97f8-21f89708018e
# ╟─4bb7d51f-ef33-4228-a39f-0f39a709899a
# ╠═4e5fccf3-d7d3-4e7f-be3f-4eb24f5a4239
# ╠═ca796498-1218-4109-b90f-5038ab76c629
# ╟─4c72470d-729d-480d-ac59-f65f195ec4ce
# ╠═a921b52a-9922-4346-92d6-c5cd03540381
# ╟─08df6171-cc09-4239-b458-2a8ba87d69ca
# ╟─73b1b527-d820-4531-92a1-73df0a492b0b
# ╠═ceb4fa2c-3981-466d-b0f8-0468bb264c68
# ╟─5b5ff634-c6bf-4fe9-b442-8e9821156e75
# ╠═fe677804-cf1f-4413-abe4-684943d68d5f
# ╠═08a9d538-1eba-11ed-1ac4-c1b15c5bb280
# ╠═f89f32a4-2d1d-4905-90d9-6746576bf432
# ╠═ea560410-5c83-4fe2-afd5-03d09af0685d
# ╠═7463bf07-e164-4580-9a0f-581625994762
# ╠═72460280-bbfa-492d-8270-e44293266e09
# ╠═bb47b3f0-b4d7-43d8-945f-ee94e9def8f6
# ╠═e4af41d5-754a-466f-a26c-cb2a6c4e160c
# ╠═72bcf653-3ec8-43ad-90a1-e08f7b063b33
# ╠═af516c3f-6fcf-4db7-a89f-30b1e813863c
# ╠═83fee348-2ee5-4612-b909-f73a800b9cb6
# ╠═763c51c1-a1a8-4262-8d40-14799603279a
# ╠═46b2299a-3c88-4f5f-9f46-3c773c70b189
# ╠═7bc5e477-37b0-435e-9224-152bc562ff04
# ╠═2de72be6-3ebd-4282-8ab0-d51bbe727aea
# ╠═9d7518a1-95a5-42fc-b268-87b100c3e96d
# ╠═c056d31e-6351-4d76-b484-fb8013c20b67
# ╠═2c8e6f69-f0b9-43a1-b18e-94a4e5de993d
# ╠═c3af66e7-7e7c-4af9-90cc-8ba7382b2791
# ╠═bcd7bbc4-90c8-492c-903d-b9559d8a2eb5
# ╠═b899d763-5636-4142-b93b-30b578144d01
# ╠═e2c0fbd0-1cdd-4d60-a12c-ba78d90748a1
# ╠═593279c7-1757-4c08-a19d-a3d269423040
# ╠═7ec2f302-bd04-4e0e-97be-56e93ecc361c
# ╠═dd397413-3822-4efc-816d-2fe4ad5e1629
# ╠═cfe930f7-bcd4-44c3-8b69-c25e0cb6ef13
# ╠═9ca7013a-3115-4931-b9da-034934e50295
# ╠═09f4b07a-877c-4357-bbe1-4d46606e8dec
# ╠═e3cbf4ec-8e7b-44ce-887f-b978524c5cbf
# ╠═10ed3abf-ae21-453c-ae7b-c27d0f106dae
# ╠═af40960b-1d79-42ed-a69c-2f24465cfc02
# ╠═a9d4a058-7ab0-42eb-b879-3787dcedf1c0
# ╠═a28594c7-c134-4548-9ac1-1049bbe14f4b
# ╠═c193d0f7-44d5-44ab-9a0e-1271e23b06f6
# ╠═3858bddb-34c8-476a-9d27-b2997b54271b
# ╠═1e010dfd-e31a-426d-99bc-09cb5f43ed6c
# ╠═3bde1246-2ecc-4167-9629-75a1a3e5b221
# ╠═1a7ef4a2-50e0-4a72-804e-d07bfb6d2fdc
# ╠═c999e688-6387-4b3d-951e-3ab87b230d06
# ╠═0f279169-a831-406d-ae6c-975ef16848a7
# ╠═b0cc8956-47de-4f5f-a55f-4893fe84cb49
# ╠═68913bd8-92f9-41d5-ad90-5e9b7ac75678
# ╠═bab47f30-e938-4403-8c4d-0677fcd25433
