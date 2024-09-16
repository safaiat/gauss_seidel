function [A,b] = generate_Ab(res)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dm = @(n) spdiags([-ones(n,1) ones(n,1)], [0 1], n, n);
Nx = 8*res; Ny = 8*res; Nz = 8*res;
dx1 = dm(Nx); dx2 = kron(speye(Ny), dx1);
dy1 = dm(Ny); dy2 = kron(dy1, speye(Nx));
dz1 = dm(Nz);
dx3 = kron(speye(Nz), dx2);
dy3 = kron(speye(Nz), dy2);
dz3 = kron(dz1, speye(Nx*Ny));
N = Nx*Ny*Nz;
ZM = sparse(N,N);
a = [ZM, dz3, dy3; dz3, ZM, dx3; dy3, dx3, ZM];
A = a*a';
b = ones(size(A,1),1);
end