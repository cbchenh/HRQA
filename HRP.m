%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heterogeneous Recurrence Plot                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Input:                                                               %
%         S: a sequence of states                                        %
%         k: the number of the state cardinality                         %
%         plot_idx: option to plot the HRP (1:=plot, 0:=not plot)        %
% - Output:                                                              %
%         hrp: the matrix of heterogeneous recurrence plot               %
%                                                                        %
% Authors:                                                               %
%         1. Cheng-Bang Chen    email: cbchen@chengbangchen.me           %
%         Copyright 2021, Cheng-Bang Chen, All rights reserved.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrp = HRP(S, k, plot_idx)
    S = S(:);
    dim = size(S);
    if nargin < 3
        plot_idx = 0
    end
    if nargin < 2
       k = max(S); 
    end
    if k == []
       k = max(S); 
    end
    hrp = zeros(dim(1),dim(1));
    for i = 1:k
        rec_idx = find(S == i);
        if length(rec_idx) > 1
            [r_idx, c_idx] = meshgrid(rec_idx,rec_idx);
            hrp(r_idx, c_idx) = i;
        end
    end
    if plot_idx > 0
        [x_idx, y_idx] = meshgrid(1:dim(1),1:dim(1));
        figure();
        plt = scatter(x_idx(hrp>0), y_idx(hrp>0), [], hrp(hrp>0), '.'); 
        colormap jet;
        %set(gcf,'position',[200,400,400,400]);
    end
end