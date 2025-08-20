using GrothCrypto
using Documenter

DocMeta.setdocmeta!(GrothCrypto, :DocTestSetup, :(using GrothCrypto); recursive=true)

makedocs(;
    modules=[GrothCrypto],
    authors="0xpantera <0xpantera@proton.me>",
    sitename="GrothCrypto.jl",
    format=Documenter.HTML(;
        canonical="https://0xpantera.github.io/GrothCrypto.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/0xpantera/GrothCrypto.jl",
    devbranch="main",
)
