using Documenter
using IRBEM

makedocs(
    sitename="IRBEM.jl",
    format=Documenter.HTML(),
    modules=[IRBEM],
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports,
    doctest=true
)

deploydocs(
    repo="github.com/Beforerr/IRBEM.jl",
)
