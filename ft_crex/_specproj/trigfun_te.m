function Strig = trigfun_te(~)
% Specific function for MEG_TE project to define the marker values and
% associated condition/name

% Photodiode signal
iod = 512;

Strig = [];

%-- Stim: "1", "2", "3"
Strig(1).name = 'S';
Strig(1).value = [10 20 30] + iod;

%-- Action: Thumb, index, middle, ring, little finger
Strig(2).name = 'A';
Strig(2).value = [128 256 512 1024 2048];

%-- Reward: Correct, Incorrect, Late
Strig(3).name = 'R';
Strig(3).value = [200 210 220] + iod;

%-- Stim: "1", "2", "3"
Strig(4).name = 'SAR';
Strig(4).value = [10 20 30] + iod;


