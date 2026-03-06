using Documenter
using Documenter.Remotes
using GrothAlgebra
using GrothCurves
using GrothProofs

DocMeta.setdocmeta!(GrothAlgebra, :DocTestSetup, :(using GrothAlgebra); recursive=true)
DocMeta.setdocmeta!(GrothCurves, :DocTestSetup, :(using GrothCurves); recursive=true)
DocMeta.setdocmeta!(GrothProofs, :DocTestSetup, :(using GrothProofs); recursive=true)

const PAGES = [
    "Home" => "index.md",
    "Getting Started" => "getting-started.md",
    "Groth16 End-to-End" => "groth16-e2e.md",
    "Packages" => [
        "Groth Algebra" => "algebra.md",
        "Groth Curves" => "curves.md",
        "Groth Proofs" => "proofs.md",
    ],
    "References" => [
        "Implementation Notes" => "implementation-notes.md",
        "Implementation vs Arkworks" => "implementation-vs-arkworks.md",
        "RareSkills Map" => "rareskills-map.md",
        "Benchmarks" => "benchmarks.md",
    ],
]

makedocs(
    sitename="Groth.jl",
    format=Documenter.HTML(
        prettyurls=true,
        collapselevel=1,
        inventory_version="dev",
    ),
    modules=[GrothAlgebra, GrothCurves, GrothProofs],
    pages=PAGES,
    clean=true,
    checkdocs=:none,
    repo=Remotes.GitHub("0xpantera", "Groth.jl"),
    meta=Dict(
        :description => "Groth.jl provides algebra, curve, and Groth16 proof tooling in Julia.",
    ),
)

deploydocs(
    repo="github.com/0xpantera/Groth.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    push_preview=true,
)
