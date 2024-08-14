include("filter_setup.jl")

multiwfs.reset!(cheb4_high)
for x in rand(10)
    @show output!(cheb4_high, x)
end