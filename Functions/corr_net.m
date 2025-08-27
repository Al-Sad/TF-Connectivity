% Copyright (c) 2025 Mohammad Al-Sa'd
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% Email: ext-mohammad.al-sad@hus.fi, mohammad.al-sad@helsinki.fi
%        mohammad.al-sad@tuni.fi, alsad.mohamed@gmail.com
%
% The following reference should be cited whenever this script is used:
% Al-Sa'd, M., & Boashash, B. (2025). Time-Frequency Functional Brain 
% Networks: Application to Newborn Seizure Analysis. 
% IntechOpen. doi: 10.5772/intechopen.1011395
%
% Last Modification: 27-August-2025
%
% Description:
% This function computes the conventional Pearson's and Spearman's 
% correlations between two signals to estimate amplitude-based connectivity.

function [pearson, spearman] = corr_net(x, alpha)
Nc = size(x,2);
[pearson, p] = corr(x,x,"Type",'Pearson');
pearson(p > alpha) = nan;
[spearman, p] = corr(x,x,"Type",'Spearman');
spearman(p > alpha) = nan;
pearson(logical(eye(Nc)))  = nan;
spearman(logical(eye(Nc))) = nan;
pearson = abs(pearson);
spearman = abs(spearman);
end
