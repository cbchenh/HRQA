clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Case 1:\n')
% Define your sequence, k, and alpha
S = [1,4,8,2,4,8,7,8,5,6,3,2,1,1,2,4,3,5];
k = 8;
alpha = 0.99 / (1 + sin(pi/k));

% Q-learning parameters
num_episodes = 5000;
epsilon = 0.01;
learning_rate = 0.05;
patience = 10000;
tol = 1e-4; %tolerance threshold for significant improvement
show_iter = 500;

% Run Q-learning to get the best permutation
[best_perm, best_reward] = q_learning_label_assignment(S, k, alpha, num_episodes, epsilon, learning_rate, patience, tol, show_iter);

% Display the results
fprintf('\nOptimal permutation found: \n');
disp(best_perm);
fprintf('Maximum reward (distance sum): %.4f\n', best_reward);
fprintf('-----------------------------------------------------\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 2:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Case 2:\n')
% Read sequence from file 'demo_seq.csv'
S1 = csvread('demo_seq.csv');
k1 = max(S1);
alpha = 0.99 / (1 + sin(pi/k));
show_iter = 1000;

% Q-learning parameters
num_episodes = 15000;
epsilon = 0.01;
learning_rate = 0.05;
patience = 5000;
tol = 1e-4; %tolerance threshold for significant improvement
% Run Q-learning to get the best permutation
[best_perm1, best_reward1] = q_learning_label_assignment(S1, k1, alpha, num_episodes, epsilon, learning_rate, patience, tol, show_iter);

% Display the results
fprintf('\nOptimal permutation found: \n');
disp(best_perm1);
fprintf('Maximum reward (distance sum): %.4f\n', best_reward1);
fprintf('-----------------------------------------------------\n')