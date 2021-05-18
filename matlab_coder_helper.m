%% Matlab coder helper
input = 1-2*rand(200,1);
[p,o] = setupFastNLS(numel(input), 20, [100, 1000], 44100, input);