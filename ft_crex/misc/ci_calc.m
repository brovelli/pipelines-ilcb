function ci = ci_calc(vals, mdim)
% Compute confidence interval for values vals in dimension mdim
%
%-CREx-180918
if nargin < 2
    mdim = 1;
end
% Confidence interval parameters
pci = 0.95;
alpha = 1 - pci;

% Number of not-NaN values use to compute the mean
Nnn = sum(~isnan(vals), mdim);

% Confidence intervals 
t_ci = tinv(1-alpha/2, Nnn-1);       

% Coef to apply to the std (ci = sigma.*coef)
cci = t_ci./sqrt(Nnn);

% Compute the std across dimension mdim 
ci = cci.*nanstd(vals, [], mdim);