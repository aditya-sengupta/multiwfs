# Define the regex pattern
pattern = r"\(([\d.]+),\s*([\d.]+)\)"

# Test string
test_str = "(2.0, 50.0)"

# Extract numbers using the regex
m = match(pattern, test_str)
if !isnothing(m)
    num1 = parse(Float64, m.captures[1])
    num2 = parse(Float64, m.captures[2])
    println("Extracted numbers: ", num1, ", ", num2)
end