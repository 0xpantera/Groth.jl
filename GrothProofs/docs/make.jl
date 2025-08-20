using GrothProofs
using Documenter

DocMeta.setdocmeta!(GrothProofs, :DocTestSetup, :(using GrothProofs); recursive=true)

makedocs(;
    modules=[GrothProofs],
    authors="0xpantera <0xpantera@proton.me>",
    sitename="GrothProofs.jl",
    format=Documenter.HTML(;
        canonical="https://0xpantera.github.io/GrothProofs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/0xpantera/GrothProofs.jl",
    devbranch="main",
)
