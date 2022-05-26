clear all;
close all;
%% Preallocs
%Data             = load("mycielskian8.mat");
Data             = load("GD96_c.mat");
%Data             = load("delaunay_n12.mat");
AdjacencyMatrix0 = Data.Problem.A;
Graph            = graph(AdjacencyMatrix0);

figure(1)
plot(Graph, 'NodeColor', 'k')

[ny0, nx0]    = size(AdjacencyMatrix0);
[U0 L0 V0]    = eigs(AdjacencyMatrix0);
%% Experiment: Katz centrality for each node as function of alpha
grid_size_K = 9;    
alphas      = linspace(0, 1/max(abs(diag(L0))), grid_size_K);
alphas      = alphas(2:end-1);

kvcent    = zeros(ny0, length(alphas));

for index = 1:length(alphas)
    kvcent(:, index) = KatzCentrality(AdjacencyMatrix0, alphas(index));
end 

figure(2)
for i = 1:length(AdjacencyMatrix0)
    plot(alphas, kvcent(i,:))
    hold on
end
%ax.XTickLabel = string(alphas);
grid on
hold off
%% Experiment: NBT centrality for each node as function of t
grid_size_NBT = 100;

[ny nx]  = size(AdjacencyMatrix0);
assert(nx == ny, 'Matrix is non-square');
    
I        = speye(ny, nx);
e        = ones(ny, 1);
D        = diag(sum(AdjacencyMatrix0));
Lam      = polyeig(I, AdjacencyMatrix0, (D-I));

min(abs(Lam))

ts            = linspace(0, min(abs(Lam)), grid_size_NBT);
ts            = ts(2:end-1);

nbtwcent      = zeros(ny0, length(ts));

for index = 1:length(ts)
    nbtwcent(:, index) = NBTCentrality(I, D, AdjacencyMatrix0, e, ts(index));
end

figure(3)
for i = 1:length(AdjacencyMatrix0)
    %index = i * 10 + 1;
    plot(ts, nbtwcent(i,:))
    hold on
end
%ax.XTickLabel = string(ts);
grid on
hold off
%% Katz-Centrality with given Adjacency matrix and alpha parameter
function katz_cen = KatzCentrality(A, alpha)
    [ny nx]  = size(A);
    assert(nx == ny, 'Matrix is non-square'); 

    I        = eye(ny, nx);
    e        = ones(ny, 1);
    katz_cen = (I - alpha*A)\e;
end
%% NBT centrality with given Adjacency matrix and t parameter
function nbt_cen = NBTCentrality(I, D, A, e, t)
    Mt       = I - t*A + (D - I) * t^2;
    nbt_cen  = (1 - t^2) * (Mt\e);
end
