%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo_C: Windowed RHRQA for a Categorical / Integer Sequence             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline:                                                               %
%   1. Load (or generate) a categorical integer sequence  S              %
%   2. Validate S and set K                                               %
%   3. Slide a window over S                                              %
%   4. Compute RHRQA statistics for each window  ->  one row per window   %
%   5. Save results table to CSV + visualise key features                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% =========================================================================
%% 0.  User Parameters
%% =========================================================================
win_size  = 200;    % window length (# samples)
step_size = 100;    % step between windows (50 % overlap)
r_order   = 1;      % local-cluster order  (1 = 1-step, 2 = 2-step …)

% --- Input source: set INPUT_MODE to 'file' or 'generate' ----------------
INPUT_MODE = 'generate';   % <-- change to 'file' to load your own data

% If INPUT_MODE == 'file', set the path to your CSV (one column, integers):
INPUT_FILE = 'my_sequence.csv';

%% =========================================================================
%% 1.  Load or Generate the Categorical Sequence
%% =========================================================================
if strcmp(INPUT_MODE, 'file')
    % -----------------------------------------------------------------------
    % Load from CSV – expects a single column of positive integers
    % -----------------------------------------------------------------------
    fprintf('[1] Loading sequence from  %s ...\n', INPUT_FILE);
    raw = readmatrix(INPUT_FILE);
    S   = round(raw(:,1));          % take first column, force integers
    S   = S(S > 0);                 % drop any zero / negative entries
    S   = reshape(S, [], 1);

else
    % -----------------------------------------------------------------------
    % Generate a synthetic 4-state sequence that mimics two regimes:
    %   Regime A (first half)  : states {1,2} dominate
    %   Regime B (second half) : states {3,4} dominate
    % -----------------------------------------------------------------------
    fprintf('[1] Generating synthetic 4-state sequence ...\n');
    rng(42);
    N     = 2000;
    K_gen = 4;

    half = N / 2;

    % Regime A: Markov chain biased toward states 1-2
    PA = [0.45 0.35 0.10 0.10;   % from state 1
          0.35 0.40 0.15 0.10;   % from state 2
          0.15 0.15 0.40 0.30;   % from state 3
          0.10 0.10 0.30 0.50];  % from state 4

    % Regime B: Markov chain biased toward states 3-4
    PB = [0.50 0.30 0.10 0.10;
          0.30 0.45 0.15 0.10;
          0.10 0.10 0.35 0.45;
          0.10 0.10 0.40 0.40];
    PB = [0.10 0.10 0.40 0.40;
          0.10 0.15 0.35 0.40;
          0.40 0.35 0.15 0.10;
          0.40 0.30 0.20 0.10];

    S = zeros(N, 1);
    S(1) = randi(K_gen);
    for i = 2 : half
        P_row  = PA(S(i-1), :);
        S(i)   = find(rand < cumsum(P_row), 1, 'first');
    end
    for i = half+1 : N
        P_row  = PB(S(i-1), :);
        S(i)   = find(rand < cumsum(P_row), 1, 'first');
    end

    fprintf('   Generated %d samples, K = %d states\n', N, K_gen);
end

%% =========================================================================
%% 2.  Validate Sequence and Set K
%% =========================================================================
S = reshape(S, [], 1);
assert(all(S == floor(S)) && all(S >= 1), ...
    'S must contain positive integers (1-based state labels).');

K = max(S);
N = length(S);
fprintf('[2] Sequence length N = %d,  states K = %d\n', N, K);

% -- Figure 1 : raw symbolic sequence
figure('Color','w','Name','Categorical Sequence');
plot(S, '.', 'MarkerSize', 3, 'Color', [0.2 0.4 0.8]);
xlabel('Sample index');  ylabel('State label');
title(sprintf('Input Categorical Sequence  (N=%d, K=%d)', N, K));
ylim([0, K+1]);  grid on;

% -- Figure 2 : state histogram
figure('Color','w','Name','State Histogram');
histogram(S, K, 'FaceColor',[0.2 0.6 0.5], 'EdgeColor','w');
xlabel('State');  ylabel('Count');
title('State Frequency Distribution');
xticks(1:K);  grid on;

%% =========================================================================
%% 3.  Sliding-Window RHRQA
%% =========================================================================
fprintf('[3] Running windowed RHRQA  (win=%d  step=%d  r=%d) ...\n', ...
        win_size, step_size, r_order);

starts = 1 : step_size : (N - win_size + 1);
nWin   = length(starts);
fprintf('   Total windows: %d\n', nWin);

% --- dry run to learn feature width ---------------------------------------
[~, HRR0, ~, ~, ~, ~, ~, ~] = RHRQA(S(1:win_size), K, r_order);
nFeats = length(HRR0);      % = K^r_order

% --- pre-allocate: [WinID, StartIdx, EndIdx, 7*nFeats] -------------------
nCols   = 3 + 7 * nFeats;
results = zeros(nWin, nCols);

print_every = max(1, floor(nWin / 10));

for w = 1 : nWin
    seg = S(starts(w) : starts(w) + win_size - 1);

    [~, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini] = ...
        RHRQA(seg, K, r_order);

    row = [HRR', HMean', HVar', HSkew', HKurt', HEnt', HGini'];
    row(isnan(row)) = 0;

    results(w, :) = [w, starts(w), starts(w)+win_size-1,  row];

    if mod(w, print_every) == 0 || w == nWin
        fprintf('   Window %3d / %d  (start=%d)\n', w, nWin, starts(w));
    end
end

%% =========================================================================
%% 4.  Build Column Names & Save to CSV
%% =========================================================================
metrics    = ["HRR","HMean","HVar","HSkew","HKurt","HEnt","HGini"];
feat_names = cell(1, 7 * nFeats);
for m = 1:7
    for f = 1:nFeats
        feat_names{(m-1)*nFeats + f} = sprintf('%s_%d', metrics(m), f);
    end
end

col_names = [{'WinID','StartIdx','EndIdx'}, feat_names];
T         = array2table(results, 'VariableNames', col_names);
outfile   = 'demo_c_statistics.csv';
writetable(T, outfile, 'WriteRowNames', false);
fprintf('[4] Results saved  ->  %s   (%d rows x %d cols)\n', ...
        outfile, nWin, nCols);

% Preview in command window
fprintf('\n--- Preview: first 3 windows (cols 1–%d) ---\n', min(10,nCols));
disp(results(1:min(3,nWin), 1:min(10,nCols)));

%% =========================================================================
%% 5.  Visualise Key Features Across Windows
%% =========================================================================
win_centers = results(:,2) + win_size/2;

% ---- 5a : HRR heatmap (all states across windows) ------------------------
HRR_mat = results(:, 3+(1:nFeats));        % nWin x nFeats

figure('Color','w','Name','HRR Heatmap');
imagesc(HRR_mat');
colormap(hot);  colorbar;
xlabel('Window index');  ylabel('State');
title('HRR per State per Window');
yticks(1:K);  yticklabels(arrayfun(@(k) sprintf('s%d',k), 1:K, 'uni',false));
set(gca,'YDir','normal');

% ---- 5b : HEnt heatmap (all states across windows) -----------------------
HEnt_mat = results(:, 3+5*nFeats+(1:nFeats));   % nWin x nFeats

figure('Color','w','Name','HEnt Heatmap');
imagesc(HEnt_mat');
colormap(parula);  colorbar;
xlabel('Window index');  ylabel('State');
title('HEnt per State per Window');
yticks(1:K);  yticklabels(arrayfun(@(k) sprintf('s%d',k), 1:K, 'uni',false));
set(gca,'YDir','normal');

% ---- 5c : line plots of aggregate statistics (mean across states) --------
HRR_mean  = mean(HRR_mat,  2);
HMean_mat = results(:, 3+nFeats+(1:nFeats));
HMean_avg = mean(HMean_mat, 2);
HEnt_avg  = mean(HEnt_mat,  2);

figure('Color','w','Name','Aggregate RHRQA Over Windows');

subplot(3,1,1);
plot(win_centers, HRR_mean, 'b-o','MarkerSize',3,'LineWidth',1.5);
xlabel('Window centre (sample)');  ylabel('Mean HRR');
title('Mean HRR across all states');  grid on;

subplot(3,1,2);
plot(win_centers, HMean_avg, 'r-o','MarkerSize',3,'LineWidth',1.5);
xlabel('Window centre (sample)');  ylabel('Mean HMean');
title('Mean HMean across all states');  grid on;

subplot(3,1,3);
plot(win_centers, HEnt_avg, 'g-o','MarkerSize',3,'LineWidth',1.5);
xlabel('Window centre (sample)');  ylabel('Mean HEnt');
title('Mean HEnt across all states');  grid on;

sgtitle('Windowed RHRQA Statistics (averaged across states)');

fprintf('[5] Done.\n');
