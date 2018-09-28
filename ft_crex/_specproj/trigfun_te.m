function Strig = trigfun_te(~)
% Specific function for MEG_TE project to define the marker values and
% associated condition/name

% Photodiode signal
iod = 512;

Strig = [];

%-- Stim: "1", "2", "3"
Strig(1).name = 'S';
Strig(1).value = [10 20 30] + iod;
% Function to transform trigger code to be kept in trl definition
Strig(1).fun = @(x) (x-512)/10;

%-- Action: Thumb, index, middle, ring, little finger
Strig(2).name = 'A';
Strig(2).value = [128 256 512 1024 2048];
Strig(2).fun = @(x) log2(x) - 6;

%-- Reward: Incorrect, Correct, Late
Strig(3).name = 'R';
Strig(3).value = [200 210 220] + iod;
Strig(3).fun = @(x)(x-512)/10 - 20;

%-- Stim: "1", "2", "3"
Strig(4).name = 'SAR';
Strig(4).value = [10 20 30] + iod;
Strig(4).fun = @(x) (x-512)/10;

