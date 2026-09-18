%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularized Heterogeneous Recurrence Analysis Quantification 2         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Description:                                                         %
%      Input a sequence of integers, then provide a set of HRQAs &       %
%      corresponding IFS address. Includes HJSD, HMedian, and HIQR,      %
%      with mathematically rigorous epsilon stability.                   %
% - Input:                                                               %
%         s: a sequence of integers                                      %
%         K: number of states (default: max(s))                          %
%         r: zoom in order (starts from 0, default: 1)                   %
%         a: iterated parameter (default:0.9999*sin(pi/K)/(1+sin(pi/K))) %
%         B: # of sections for calculating entropy & Gini (default: 20)  %
% - Output:                                                              %
%         IdxM: labels of states                                         %
%         HRR: Heterogeneous Recurrence Rate                             %
%         HMean: Heterogeneous Recurrence Mean                           %
%         HVar: Heterogeneous Recurrence Variance                        %
%         HSkew: Heterogeneous Recurrence Skewness                       %
%         HKurt: Heterogeneous Recurrence Kurtosis                       %
%         HEnt: Heterogeneous Recurrence Entropy                         %
%         HGini: Heterogeneous Recurrence Gini Index                     %
%         HJSD: Heterogeneous Recurrence Jensen-Shannon Divergence       %
%         HMedian: Heterogeneous Recurrence Median                       %
%         HIQR: Heterogeneous Recurrence Interquartile Range             %
%                                                                        %
% Authors:                                                               %
%         Cheng-Bang Chen    email: cbchen@chengbangchen.me              %
%         Copyright 2021-2026, Cheng-Bang Chen, All rights reserved.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini, HJSD, HMedian, HIQR] = RHRQA2(s, K, r, a, B)
    if isempty(s)
        error('The input sequence, s, is empty.')
    end
    s = reshape(s,[length(s),1]);
    
    % Default values of K, r, a, B
    if nargin < 5
        if nargin < 4
           if nargin < 3
              if nargin < 2
                 K = max(s); 
              end
              r = 1;
           end
           a = 0.9999*sin(pi/K)/(1+sin(pi/K));
        end
        B=20;
    end
    
    if K < max(s)
       error('Error: K < max(s)');
    end
    
    % Numerical stability constant for division and logarithms
    eps_val = 1e-12;
    
    % IFS address: Cv
    Cv = IFS(s, K, a);
    Cv = Cv(2:end,:);
    
    % Generate combinations of all local clusters: M
    IdxM = [];
    if(r == 0)
        IdxM = 0;
        Idx_set = 1:length(s);
        N = size(s,1);
        H_bar = [];
        for i = 1:K
            H_bar = [H_bar , sum(s==i)];
        end
        
        % HRR (Squared, consistent with 2D Recurrence Plot density)
        HRR = sum((H_bar./N).^2);
        
        LC_D = ((1/a)^(r))*pdist(Cv(Idx_set,:));
        M_k = nchoosek(length(s),2);
        
        HMean = mean(LC_D);
        HVar = sum((LC_D - HMean).^2) / M_k;
        
        % Skew & Kurt with epsilon stability
        std_val = sqrt(HVar) + eps_val;
        HSkew = sum(((LC_D - HMean) ./ std_val).^3) / M_k; 
        HKurt = sum(((LC_D - HMean) ./ std_val).^4) / M_k;
        
        % Entropy & Gini with epsilon stability
        tmp_counts = hist(LC_D, B);
        p = tmp_counts ./ sum(tmp_counts);
        HEnt = -sum(p .* log2(p + eps_val));
        HGini = 1 - sum(p.^2);
        
        % New Metrics for r = 0
        HMedian = median(LC_D);
        HIQR = prctile(LC_D, 75) - prctile(LC_D, 25);
        HJSD = 0; % Global against itself is 0
        
    else
        for i = 1:r
            tmp_c = [];
            tmp_e = [];
            for j = 1:K
                tmp_e = [tmp_e; j*ones(K^(r-i),1)];
            end
            for k = 1:(K^(i-1))
                tmp_c = [tmp_c; tmp_e]; 
            end
            IdxM(:,r-i+1) = tmp_c;
        end
        
        % Generate real local clusters
        LC = [];
        for i = 1:r
            LC(:,i) = s(i:end-r+i);
        end
        N = size(s,1); % Total sequence length
        
        % Unique local clusters
        UC = unique(LC, 'row');
        
        % Result matrix: [HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini, HJSD, HMedian, HIQR]
        UC_Res = NaN(length(UC), 10); 
        
        % Data storage for Pass 2 (HJSD)
        all_LC_D = cell(length(UC), 1); 
        g_min = inf;
        g_max = -inf;
        
        % === PASS 1: Calculate metrics 1-7, 9-10 & collect global min/max ===
        for i = 1:size(UC,1)
            Idx_set = find(ismember(LC, UC(i,:), 'row'));
            Idx_set = Idx_set + (r-1)*ones(size(Idx_set));
            
            % H_bar: Cardinality of local cluster
            H_bar = length(Idx_set);
            
            % HRR (Squared)
            UC_Res(i,1) = (H_bar/N)^2;
            
            LC_ifs = Cv(Idx_set,:);
            LC_D = ((1/a)^(r))*pdist(LC_ifs);
            
            if H_bar > 1        
                % Store distances for Pass 2 pooling
                all_LC_D{i} = LC_D;
                M_k = nchoosek(H_bar, 2);
                
                % HMean and HVar
                HMean_val = mean(LC_D);
                HVar_val = sum((LC_D - HMean_val).^2) / M_k;
                UC_Res(i,2) = HMean_val;
                UC_Res(i,3) = HVar_val;
                
                % HSkew and HKurtosis (with epsilon stability)
                std_val = sqrt(HVar_val) + eps_val;
                UC_Res(i,4) = sum(((LC_D - HMean_val) ./ std_val).^3) / M_k;
                UC_Res(i,5) = sum(((LC_D - HMean_val) ./ std_val).^4) / M_k;
                
                % HEnt and HGini (Original local bins, with epsilon stability)
                tmp_counts = hist(LC_D, B);
                p = tmp_counts ./ sum(tmp_counts);
                UC_Res(i,6) = -sum(p .* log2(p + eps_val));
                UC_Res(i,7) = 1 - sum(p.^2);
                
                % HMedian and HIQR
                UC_Res(i,9) = median(LC_D);
                UC_Res(i,10) = prctile(LC_D, 75) - prctile(LC_D, 25);
                
                % Track global bounds for HJSD shared support
                if min(LC_D) < g_min
                    g_min = min(LC_D);
                end
                if max(LC_D) > g_max
                    g_max = max(LC_D);
                end
            end
        end
        
        % === PASS 2: HJSD using global pooled reference histogram ===
        if g_max > g_min
            edges = linspace(g_min, g_max, B+1);
            total_counts = zeros(1, B);
            cluster_counts = zeros(length(UC), B);
            
            % Generate counts for each valid cluster on the global grid
            for i = 1:size(UC,1)
                if ~isempty(all_LC_D{i})
                    [counts, ~] = histcounts(all_LC_D{i}, edges);
                    cluster_counts(i, :) = counts;
                    total_counts = total_counts + counts;
                end
            end
            
            % Estimate global reference q_b
            q_ref = total_counts / sum(total_counts);
            
            % Calculate divergence for each valid cluster
            for i = 1:size(UC,1)
                if ~isempty(all_LC_D{i})
                    p_counts = cluster_counts(i, :);
                    p_sum = sum(p_counts);
                    if p_sum > 0
                        p = p_counts / p_sum;
                        m = 0.5 * (p + q_ref);
                        
                        % HJSD formula with epsilon
                        term1 = 0.5 * sum(p .* log2((p + eps_val) ./ (m + eps_val)));
                        term2 = 0.5 * sum(q_ref .* log2((q_ref + eps_val) ./ (m + eps_val)));
                        UC_Res(i,8) = term1 + term2;
                    else
                        UC_Res(i,8) = 0;
                    end
                end
            end
        else
            % Edge case: All pairwise distances across sequence are identical
            for i = 1:size(UC,1)
                if ~isempty(all_LC_D{i})
                    UC_Res(i,8) = 0;
                end
            end
        end

        % Map Unique Cluster results back to the full IdxM matrix
        M_Res = NaN(size(IdxM,1), 10);
        for i = 1:length(UC)
            ins_idx = find(ismember(IdxM, UC(i,:), 'row'));
            if ~isempty(ins_idx)
                M_Res(ins_idx,:) = UC_Res(i,:);
            end
        end
        
        M_Res(isnan(M_Res(:,1)), 1) = 0;
        HRR = M_Res(:,1);
        HMean = M_Res(:,2); 
        HVar = M_Res(:,3); 
        HSkew = M_Res(:,4); 
        HKurt = M_Res(:,5); 
        HEnt = M_Res(:,6);
        HGini = M_Res(:,7);  
        HJSD = M_Res(:,8);
        HMedian = M_Res(:,9);
        HIQR = M_Res(:,10);
    end
end