# An elliptic equation in one dimmension

Consider the following problem. Given \(f \in C([0, 1])\), find a function \(u\) satisfying
\(\displaystyle \left\{ \begin{array}{l} -u''(x) + u(x) = f(x) \textrm{ in } (0, 1) \\ u(0) = u(1) = 0.\end{array}\right. \)

```julia
using Literate
Literate.markdown("ex1-lagrange.jl", "."; documenter=false)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

