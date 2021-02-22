# # An equation of second order in one dimmension
#
# Consider the following problem. Given \\(f \in C([0, 1])\\), find a function \\(u\\) satisfying
# -u'' + u = f in (0, 1)
# u(0) = u(1) = 0


using SymFEL
using SymPy
using LinearAlgebra
using SparseArrays

using PyPlot
close("all")

# discretization parameters
N = 101; # number of nodes
dx = 1 / (N - 1); # discretization step
nodes = range(0, stop=1, length=N);
bound_nodes = [1, N];

# exact solution
u_exact = nodes .* sin.(pi * nodes);
# right hand 
f = (1 + pi^2) * nodes .* sin.(pi * nodes) - 2 * pi * cos.(pi * nodes);

# elementary matrices
elem_K = SymFEL.get_lagrange_em(1, 1, 1);
elem_M = SymFEL.get_lagrange_em(1, 0, 0);
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

plot(nodes, u_exact);
plot(nodes, u);
legend(["Exact solution", "Approximate solution"]);


err = u - u_exact;
println("L2 error = ", sqrt(err' * (M * err)))
println("H1 error = ", sqrt(err' * (K * err)))


# error convergence
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
