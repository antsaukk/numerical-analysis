clear all;
close all;
%% Preallocs
Data             = load("mycielskian8.mat");
AdjacencyMatrix0 = Data.Problem.A;
Graph            = graph(AdjacencyMatrix0);

figure(1)
plot(Graph, 'NodeColor', 'k')

[ny0, nx0]    = size(AdjacencyMatrix0);
[U0 L0 V0]    = eigs(AdjacencyMatrix0);
%% Experiment: Katz centrality for each node as function of alpha
grid_size_K = 9;    
alphas      = linspace(0, 1/max(abs(diag(L0))), grid_size_K);
alphas      = alphas(2:end-1); %tai end

kvcent    = zeros(ny0, length(alphas));

for index = 1:length(alphas)
    kvcent(:, index) = KatzCentrality(AdjacencyMatrix0, alphas(index));
end 

figure(2)
for i = 1:191
    %index = i * 10 + 1;
    plot(alphas, kvcent(i,:), '--o')
    hold on
end
%ax.XTickLabel = string(alphas);
grid on
hold off
%% Experiment: NBT centrality for each node as function of t
grid_size_NBT = 9;    
ts            = linspace(0, 1, grid_size_NBT);
ts            = ts(2:end-1); %tai end

nbtwcent      = zeros(ny0, length(ts));

for index = 1:length(ts)
    nbtwcent(:, index) = NBTCentrality(AdjacencyMatrix0, ts(index));
end

figure(3)
for i = 1:191
    %index = i * 10 + 1;
    plot(ts, nbtwcent(i,:), '--o')
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
    %katz_cen = inv((I - alpha*A) - I)*e;
end
%% NBT centrality with given Adjacency matrix and alpha parameter
function nbt_cen = NBTCentrality(A, t)
    [ny nx]  = size(A);
    assert(nx == ny, 'Matrix is non-square');
    
    I        = eye(ny, nx);
    e        = ones(ny, 1);
    D        = diag(sum(A));
    Mt       = I - A*t + (D - I)*t^2;
    nbt_cen  = (1 - t^2) * (Mt\e);
end
