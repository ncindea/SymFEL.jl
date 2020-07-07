using SymPy, Documenter, FEMTools

makedocs(
    sitename="FEMTools.jl",
    format = Documenter.HTML(mathengine=Documenter.Writers.HTMLWriter.MathJax()),
    
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md"
        ]
    )

