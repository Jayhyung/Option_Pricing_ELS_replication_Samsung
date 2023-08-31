function [ sol ] = TriD_Sol(diag, upper_diag, lower_diag, force, nsize )

% This function solves Tri-Diagonal linear system.


% Forward Substitution
for i=1:nsize-1
    diag(i+1) = diag(i+1) - upper_diag(i)*lower_diag(i)/diag(i);
    force(i+1) = force(i+1) - force(i)*lower_diag(i)/diag(i);
end

% Backward Substitution
sol = zeros(nsize,1);
sol(nsize) = force(nsize)/diag(nsize);

for i = nsize-1 : -1 : 1
    sol(i) = (force(i) - upper_diag(i)*sol(i+1))/diag(i);
end

