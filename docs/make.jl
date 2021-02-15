using SymPy, Documenter, FEMTools

makedocs(
    sitename="FEMTools.jl",
    format = Documenter.HTML(mathengine=Documenter.Writers.HTMLWriter.MathJax2()),
    
    pages=[
        "Home" => "index.md",
        "One dimensional examples" => "examples1.md",
        "Two dimensional examples" => "examples2.md"
        ]
    )

