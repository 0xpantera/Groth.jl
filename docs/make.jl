using Documenter
using DocumenterMermaid
using Documenter.Remotes
using GrothAlgebra
using GrothCurves
using GrothProofs

DocMeta.setdocmeta!(GrothAlgebra, :DocTestSetup, :(using GrothAlgebra); recursive=true)
DocMeta.setdocmeta!(GrothCurves, :DocTestSetup, :(using GrothCurves); recursive=true)
DocMeta.setdocmeta!(GrothProofs, :DocTestSetup, :(using GrothProofs); recursive=true)

const PAGES = [
    "Home" => "index.md",
    "Start Here" => [
        "Getting Started" => "getting-started.md",
        "Groth16 End-to-End" => "groth16-e2e.md",
    ],
    "Learning Path" => [
        "RareSkills ZK Book Map" => "rareskills-map.md",
        "Textbook To Optimized Code" => "textbook-to-optimized.md",
        "Architecture Map" => "architecture.md",
    ],
    "Packages" => [
        "Groth Algebra" => "algebra.md",
        "Groth Curves" => "curves.md",
        "Groth Proofs" => "proofs.md",
    ],
    "Engineering" => [
        "Benchmarks" => "benchmarks.md",
        "Benchmark Snapshots" => "benchmark-snapshots.md",
        "Implementation vs Arkworks" => "implementation-vs-arkworks.md",
        "Implementation Notes" => "implementation-notes.md",
    ],
    "API Reference" => "api.md",
]

makedocs(
    sitename="Groth.jl",
    format=Documenter.HTML(
        prettyurls=true,
        collapselevel=2,
        inventory_version="dev",
    ),
    modules=[GrothAlgebra, GrothCurves, GrothProofs],
    pages=PAGES,
    clean=true,
    checkdocs=:none,
    repo=Remotes.GitHub("0xpantera", "Groth.jl"),
    meta=Dict(
        :description => "Groth.jl is a Julia research platform for BN254 algebra, pairings, and Groth16.",
    ),
)

deploydocs(
    repo="github.com/0xpantera/Groth.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    push_preview=true,
)
