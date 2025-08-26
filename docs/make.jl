using Documenter
using IRBEM

DocMeta.setdocmeta!(IRBEM, :DocTestSetup, :(using IRBEM; using PrettyPrinting); recursive=true)

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
    repo="github.com/JuliaSpacePhysics/IRBEM.jl",
)
