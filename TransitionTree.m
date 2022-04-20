function [TTree, Tidx, TTreeS, TTreeSn] = TransitionTree(S, n)
    isaninteger = @(x)isfinite(x) & x==floor(x);
    if ~isaninteger(n)
        error('Number of transistions should be a non-negtive integer.')
    end
    if n < 0
       error('Number of transitions should be larger than 0.') 
    end
    l = length(S);
    if l < n
       error('Length of the sequence should be larger than the number of transitions.') 
    end
    M = zeros(l-n, n+1);
    if n > 0
        for i = 1:(n+1)
            M(:,i) = S(i:(l-n+i-1));
        end
    else
        M(:,1) = S;
    end
    TTree = unique(M, 'row');
    TTreeS = [];
    [Tidx, ~] = size(TTree);
    for i = 1:(l-n)
        [~, ~, TTreeS(i)] = intersect(M(i,:), TTree,'rows');
    end
    TTidx = 1:length(TTreeS);
    TTreeSn = TTreeS(TTidx(mod(1:length(TTreeS),(n+1))==0));
    return
end