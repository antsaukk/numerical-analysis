clear all;
close all;
%% Preallocs

%% Compute Katz-Centrality with given Adjacency matrix and alpha parameter
function katz_cen = KatzCentrality(A, alpha)
    [ny nx]  = size(A);
    assert(nx == ny, 'Matrix is non-square'); 

    I        = eye(ny, nx);
    e        = ones(ny, 1);
    katz_cen = (I - alpha*A)\e;
end
%% Compute NBT centrality with given Adjacency matrix and alpha parameter
function nbt_cen = NBTCentrality(A, t)
    [ny nx]  = size(A);
    assert(nx == ny, 'Matrix is non-square');
    
    I        = eye(ny, nx);
    e        = ones(ny, 1);

    D        = diag(sum(A));
    Mt       = I - A*t + (D - I)*t^2;

    nbt_cen  = (1 - t^2) * (Mt\e);
end
