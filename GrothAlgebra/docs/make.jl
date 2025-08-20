using GrothAlgebra
using Documenter

DocMeta.setdocmeta!(GrothAlgebra, :DocTestSetup, :(using GrothAlgebra); recursive=true)

makedocs(;
    modules=[GrothAlgebra],
    authors="0xpantera <0xpantera@proton.me>",
    sitename="GrothAlgebra.jl",
    format=Documenter.HTML(;
        canonical="https://0xpantera.github.io/GrothAlgebra.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/0xpantera/GrothAlgebra.jl",
    devbranch="main",
)
