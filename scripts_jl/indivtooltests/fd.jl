using ForwardDiff: derivative

abstract type FunctionData end

function f(x, f::FunctionData)
    println(typeof(f.a))
    return f.a^2 * x^2
end

struct NonDiffable <: FunctionData
    a::Float64
end

# this is fine, we're doing input/output diffing
# print value is Float64
derivative(x -> f(x, NonDiffable(3)), 2)

# this isn't: ForwardDiff depends on type overwriting
derivative(a -> f(2, NonDiffable(a)), 3)

# we fix this by allowing our struct to take on whatever type we need it to
struct Diffable{D} <: FunctionData
    a::D
end

# and now we're good
# print value is ForwardDiff.Dual{ForwardDiff.Tag{var"#15#16", Int64}, Int64, 1}
derivative(a -> f(2, Diffable(a)), 3)