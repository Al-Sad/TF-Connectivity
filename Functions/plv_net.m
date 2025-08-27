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
% This function computes the following conventional single-trial 
% phase-synchronization measures between two signals: phase locking value
% (PLV), imaginary PLV (iPLV), and corrected imaginary PLV (ciPLV).

function [plv, iplv, ciplv] = plv_net(x)
[Nt, Nc] = size(x);
z = zeros(Nt, Nc);
for p = 1:Nc
    z(:,p) = angle(hilbert(x(:,p)));
end
plv   = nan(Nc);
iplv  = nan(Nc);
ciplv = nan(Nc);
for p = 1:(Nc-1)
    phi_1 = z(:,p);
    for q = (p+1):Nc
        phi_2 = z(:,q);
        e = exp(1i*(phi_1-phi_2));
        plv(p,q)   = abs(mean(e));
        iplv(p,q)  = abs(mean(imag(e)));
        ciplv(p,q) = iplv(p,q)./sqrt(1-abs(mean(real(e))).^2);
        plv(q,p)   = plv(p,q);
        iplv(q,p)  = iplv(p,q);
        ciplv(q,p) = ciplv(p,q);
    end
end
end
