function [best_perm, best_reward] = q_learning_label_assignment(S, k, alpha, max_episodes, epsilon, learning_rate, patience, tol, show_iter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
%   S: A sequence of data used for calculating points.
%   k: Number of labels (the size of the permutation to choose).
%   alpha: Scaling factor for computing points.
%   max_episodes: Maximum number of iterations allowed for learning.
%   epsilon: Probability of exploration (selecting random actions) during 
%            Q-learning.
%   learning_rate: Learning rate controlling how strongly new information
%                  overrides old information.
%   patience: Number of consecutive episodes without significant 
%             improvement after which the learning stops 
%             (convergence criterion).
%   tol: Tolerance for improvement, used to determine significant increases 
%        in reward.
%   show_iter: show progress every show_iter iteration
%
% Output:
%   best_perm: The best permutation found.
%   best_reward: The corresponding maximum reward.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    labels = 1:k;
    Q = containers.Map('KeyType','char','ValueType','double');
    best_reward = -inf;
    best_perm = [];
    no_improve_counter = 0;

    for episode = 1:max_episodes
        available = labels;
        perm = zeros(1,k);
        state = '';

        for step = 1:k
            state_key = state;
            if rand < epsilon || ~isKey(Q, state_key)
                action = datasample(available,1);
            else
                q_values = arrayfun(@(a) get_Q(Q,[state num2str(a)]), available);
                [~, idx] = max(q_values);
                action = available(idx);
            end
            perm(step) = action;
            available = setdiff(available, action);
            state = [state num2str(action)];
        end

        points = compute_points(S, perm, alpha, k);
        reward = compute_reward(points);

        state_key = state;
        if ~isKey(Q,state_key)
            Q(state_key) = 0;
        end
        Q(state_key) = Q(state_key) + learning_rate * (reward - Q(state_key));

        if reward > best_reward + tol
            best_reward = reward;
            best_perm = perm;
            no_improve_counter = 0; % reset counter when improvement is significant
        else
            no_improve_counter = no_improve_counter + 1;
        end

        if no_improve_counter >= patience
            fprintf('Convergence reached at episode %d.\n', episode);
            break;
        end

        if mod(episode,show_iter)==0
            fprintf('Episode %d, Current best reward: %.2f\n', episode, best_reward);
        end
    end
end
