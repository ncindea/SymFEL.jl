# One dimensional examples

## A second order equation

Consider the following equation:
```math
	\left\{
	\begin{array}{l}
	-\dfrac{\partial^2}{\partial x^2}u(x) + u(x) = f(x) \textrm{ for } x \in (0, 1) \\
	u(0) = u(1) = 0.
	\end{array}
	\right.
```
```@example 0
using SymFEL
using SymPy
using LinearAlgebra
using SparseArrays

using PyPlot
close("all")

nothing #hide
```

### Discretization parameters
```@example 0
N = 101; # number of nodes
dx = 1 / (N - 1); # discretization step
nodes = range(0, stop=1, length=N);
bound_nodes = [1, N];
nothing #hide
```

### Exact solution
We take \\( f = (1 + \pi^2)x \sin(\pi x) - 2\pi \cos(\pi x)\\). For this choice of second member the exact solution is given by \\( u(x) = x \sin(\pi x)\\).

```@example 0
# exact solution
u_exact = nodes .* sin.(pi * nodes);
# right hand 
f = (1 + pi^2) * nodes .* sin.(pi * nodes) - 2 * pi * cos.(pi * nodes);
nothing #hide
```

### Finite element matrices
```@example 0
elem_K = SymFEL.get_lagrange_em(1, 1, 1);
elem_M = SymFEL.get_lagrange_em(1, 0, 0);
elem_K_dx = convert(Matrix{Float64}, elem_K.subs(h, dx));
elem_M_dx = convert(Matrix{Float64}, elem_M.subs(h, dx));

K = SymFEL.assemble_1d_FE_matrix(elem_K_dx, N, intNodes=0, dof1=1, dof2=1);
M = SymFEL.assemble_1d_FE_matrix(elem_M_dx, N, intNodes=0, dof1=1, dof2=1);

F = M * f;
nothing #hide
```
### Boundary conditions
```@example 0
tgv = 1e100
A = K + M
A[bound_nodes, bound_nodes] += tgv*sparse(Matrix{Float64}(I, 2, 2));

F[bound_nodes] = zeros(2);
nothing #hide
```

### Solve the linear system
```@example 0
u = A \ F;
nothing #hide
```

### Graphic representation
```@example 0
plot(nodes, u_exact);
plot(nodes, u);
legend(["Exact solution", "Approximate solution"]);
savefig("plots/ex1-lagrange-sol.png") #hide

err = u - u_exact;
println("L2 error = ", sqrt(err' * (M * err)))
println("H1 error = ", sqrt(err' * (K * err)))
```
![](plots/ex1-lagrange-sol.png)

### Study of the error 

```@example 0
NV = 10 * 2 .^ (0:10)
EL2 = zeros(11)
EH1 = zeros(11)
i = 1
for N = NV
    global i
    dx = 1 / (N - 1); # discretization step
    nodes = range(0, stop=1, length=N);
    bound_nodes = [1, N];

    # exact solution
    u_exact = nodes .* sin.(pi * nodes);
    # right hand 
    f = (1 + pi^2) * nodes .* sin.(pi * nodes) - 2 * pi * cos.(pi * nodes);

    elem_K_dx = convert(Matrix{Float64}, elem_K.subs(h, dx));
    elem_M_dx = convert(Matrix{Float64}, elem_M.subs(h, dx));

    K = SymFEL.assemble_1d_FE_matrix(elem_K_dx, N, intNodes=0, dof1=1, dof2=1);
    M = SymFEL.assemble_1d_FE_matrix(elem_M_dx, N, intNodes=0, dof1=1, dof2=1);

    F = M * f;

    # boundary conditions
    tgv = 1e100
    A = K + M
    A[bound_nodes, bound_nodes] += tgv*sparse(Matrix{Float64}(I, 2, 2));

    F[bound_nodes] = zeros(2);

    u = A \ F;

    err = u - u_exact;

    EL2[i] = sqrt(err' * M * err)
    EH1[i] = sqrt(err' * K * err)

    i += 1
end

figure()
loglog(1 ./ NV, EL2, "d")
loglog(1 ./ NV, EH1)
loglog(1 ./ NV, (1 ./ NV).^2)
xlabel(L"$\Delta x$")
legend(["L2 Error", "H1 Error", L"$\Delta x^2$"])
savefig("plots/ex1-lagrange-error.png") #hide
nothing #hide
```
![](plots/ex1-lagrange-error.png)

## A fourth order equation

Consider the following problem. Given \\(f \in C([0, 1])\\), find a function \\(u\\) satisfying
```math
	\left\{
	\begin{array}{l}
	\dfrac{\partial^4}{\partial x^4}u(x) = f(x) \textrm{ in } (0, 1) \\
	u(0) = u(1) = 0 \\
	u^\prime(0) = \pi, \qquad u^\prime(1) = -\pi.
	\end{array}
	\right.
```

```@example 1
using SymFEL
using SymPy
using LinearAlgebra
using SparseArrays

using PyPlot
close("all")

nothing # hide
```

### Discretization
```@example 1
N = 101; # number of nodes
dx = 1 / (N - 1); # discretization step
nodes = range(0, stop=1, length=N);
bound_nodes = [1, N];
nothing # hide
```

### Exact solution
```@example 1
# exact solution
u_exact = sin.(pi*nodes);
ud_exact = pi * cos.(pi*nodes)
# right hand 
f = pi^4 * sin.(pi * nodes);
nothing # hide
```

### Finite element matrices
```@example 1
elem_K = SymFEL.get_hermite_em(3, 2, 2);
elem_M = SymFEL.get_em(3, 1, 0, 0; fe1="Hermite", fe2="Lagrange")
elem_K_dx = convert(Matrix{Float64}, elem_K.subs(h, dx));
elem_M_dx = convert(Matrix{Float64}, elem_M.subs(h, dx));

K = SymFEL.assemble_1d_FE_matrix(elem_K_dx, N, intNodes=0, dof1=2, dof2=2);
M = SymFEL.assemble_1d_FE_matrix(elem_M_dx, N, intNodes=0, dof1=2, dof2=1);

# f is approached by an P1 function
F = M * f;
nothing # hide
```

### Boundary conditions
```@example 1
tgv = 1e100
KB = copy(K)
KB[2*bound_nodes .- 1, 2*bound_nodes .- 1] += tgv*sparse(Matrix{Float64}(I, 2, 2));
KB[2*bound_nodes, 2*bound_nodes] += tgv*sparse(Matrix{Float64}(I, 2, 2));

F[2*bound_nodes .- 1] = zeros(2);
F[2] = pi * tgv
F[2*N] = -pi * tgv
nothing # hide
```

### Solve the linear system
```@example 1
u = KB \ F;
nothing #hide
```

### Graphic representation
```@example 1
plot(nodes, u_exact);
plot(nodes, u[1:2:2*N-1]);
legend(["Exact solution", "Approximate solution"]);
savefig("plots/ex2-hermite.png"); nothing # hide
```
![](plots/ex2-hermite.png)

### Error evaluation
We choose the following norm for the error \\( \|e\| = \left(\int_0^1 |f^{\prime\prime}(x)|^2 dx\right)^\frac{1}{2}\\)
```@example 1
U_exact = zeros(2*N)
U_exact[1:2:2*N-1] = u_exact
U_exact[2:2:2*N] = ud_exact
err = u - U_exact;
println("H2 error = ", sqrt(err' * K * err))
```
