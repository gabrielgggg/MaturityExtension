clear; clc;

noStates = 21;
uncondMean = 1;
coefLags = 0.9;
stds = 0.02; % / sqrt(1 - 0.9^2) * sqrt(1 - coefLags^2);

[Z, Zprob] = tauchen1d(noStates, uncondMean, coefLags, stds, 2.0);

dlmwrite('y.tab', Z, 'delimiter', '\t', 'precision', 12);
dlmwrite('PI.tab', Zprob, 'delimiter', '\t', 'precision', 12);

