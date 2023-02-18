

x = [1,2,3,4,5,6]
println(maximum(x))
println(log10.(x))

function test()
@goto lol

println("non salto")


@label lol
println("salto")
end
test()

using Plots

gr()
h = plot(1:2) |>  display
delete!(h.parent, h)

