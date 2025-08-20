using GrothCurves
using Documenter

DocMeta.setdocmeta!(GrothCurves, :DocTestSetup, :(using GrothCurves); recursive=true)

makedocs(;
    modules=[GrothCurves],
    authors="0xpantera <0xpantera@proton.me>",
    sitename="GrothCurves.jl",
    format=Documenter.HTML(;
        canonical="https://0xpantera.github.io/GrothCurves.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/0xpantera/GrothCurves.jl",
    devbranch="main",
)
