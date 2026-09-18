%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. read data from csv file
data = readmatrix('Lorenz.csv');
disp('Finish reading the csv file!')

% Create a 3x1 panel for line plots
figure;
% Plot the first column
subplot(3, 1, 1);
plot(data(:, 1));
title('X');
xlabel('Index');
ylabel('Value');

% Plot the second column
subplot(3, 1, 2);
plot(data(:, 2));
title('Y');
xlabel('Index');
ylabel('Value');

% Plot the third column
subplot(3, 1, 3);
plot(data(:, 3));
title('Z');
xlabel('Index');
ylabel('Value');


[Idx, Ub, Lb] = HAS(data, size(data,1)/4);
S = SymbG(data, Ub, Lb);
figure();
scatter3(data(:,1),data(:,2),data(:,3),5,S,'filled');title('Trajectory of Data');hold on;
PlotCell(Ub, Lb);
view(55,25);

%HRP(S,max(Idx),1);
tic
output = [];
IFS(S, size(Ub,1), [] ,1);
[IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini] = RHRQA(S, size(Ub,1),1);
% Tmp_HRQA = horzcat(IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini);
Tmp_HRQA = [HRR', HMean', HVar', HSkew', HKurt', HEnt', HGini'];
output(1,:) = Tmp_HRQA;
toc
varname = ["HRR";"HMean";"HVar";"HSkew";"HKurt";"HEnt";"HGini"];
index = repelem(varname,length(HRR));
varname = index + "_" + repmat(1:length(HRR),1,7)';
output = array2table(output,'VariableNames',varname');
outname = strcat('demo_statistics.csv');
writetable(output,outname,'WriteRowNames',false);