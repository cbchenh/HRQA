function R = RP(X, epsilon)
    N = length(X);
    R = zeros(N,N);
    for i = 1:N
       for j = i:N
          dij = norm(X(i)-X(j));
          if dij<epsilon
              R(i,j) = 1;
              R(j,i) = 1;
          end
       end
    end
    [x, y] = find(R == 1);
    figure()
    set(gcf,'Color',[1 1 1]);
    plot(x, y, '.');
    xlim([0,N]);
    ylim([0,N]);
end