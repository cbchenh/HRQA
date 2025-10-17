function points = compute_points(S, perm, alpha, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes fractal points given a sequence S, permutation perm, alpha, and k
% s: a sequence of integers
% perm: a permutation array used to rearrange or map the elements of S to 
%       IFS
% alpha: the scalar of IFS
% k: the number of categories 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(S);
points = zeros(n,2);
C = [0;0];
for t = 1:n
    St = perm(S(t));
    C = alpha * C + [cos(2*pi*St/k); sin(2*pi*St/k)];
    points(t,:) = C';
end
end
