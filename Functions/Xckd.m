%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 Boualem Boashash
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Prof. Boualem Boashash        (boualem.boashash@gmail.com)
%          Dr. Abdeldjalil Aissa-El-Bey  (abdeldjalil.aissaelbey@telecom-bretagne.eu)
%          RA: Md.F.A
%
% The following references should be cited whenever this script is used:
% [1] B. Boashash, A. Aissa-El-Bey, Multisensor Time-Frequency Signal Processing:
%     A tutorial review with illustrations in selected application areas, Digital
%     Signal Processing, In Press.
% [2] B. Boashash, A. Aissa-El-Bey, M. F. Al-Sa'd, Multisensor time-frequency
%     signal processing software Matlab package: An analysis tool for multichannel
%     non-stationary data , SoftwareX, In Press.
%
% Last Modification: 21-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Cross Compact Support Kernel Distribution (XCKD)
%
% Syntax : TFR = Xckd(s1, s2, varargin);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% s1, s2   : Input signals. These can be real or analytic.
% varargin : This is a dynamic input argument that has to be supplied in
%            the following sequence:
%            varargin{1} : The smoothing parameter C (default : 1).
%            varargin{2} : The smoothing parameter D (default : 0.1).
%            varargin{3} : The smoothing parameter E (default : 0.1).
%            varargin{4} : Number of Frequency bins.
% <OUTPUTs>
% TFR  : Reduced Interference Cross-TFD using Cross-Compact Support Kernel Distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% s1 = chirp(0:255,0.1,255,0.25);
% s2 = chirp(0:255,0.2,255,0.35);
% s = s1 + s2;
% TFD = Xckd(s, s, 0.1, 0.05, 0.05, 512);
% figure; imagesc(0:1/1023:1/2,0:255,abs(TFD')); axis xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TFR = Xckd(s1, s2, C, D, E)
%% Computing the Kernel in the Ambiguity Domain
N = length(s1);
L = N-1;
M = 2^nextpow2(N);
dopp = linspace(-0.5,0.5,N);
lagg = linspace(-M/2,M/2,M);
E = E*M;
G = zeros(M, N);
for i = 1:M
    for j = 1:N
        if(abs(dopp(j)) < D && abs(lagg(i)) < E)
            G(i,j) = exp(2*C)*exp(((C*D^2)/(dopp(j)^2-D^2)) + ((C*E^2)/(lagg(i)^2-E^2)));
        end
    end    
end
% G = G/norm(G);

%% Computing the Analytic Associate of the Input Signals S1 and S2
P = nextpow2(N); K = 2*(2^P);
z1 = fft(s1,K);           % Transforming S1 into the frequency domain
z1(1) = z1(1);            % Holding the zero frequency amplitude/energy
z1(2:K/2) = 2*z1(2:K/2);  % Accepting and doubling half of the frequencies
z1(K/2+1:K) = 0;          % Rejecting half of the frequencies
z1 = ifft(z1);            % Transforming Z into the time domain
z2 = fft(s2,K);           % Transforming S2 into the frequency domain
z2(1) = z2(1);            % Holding the zero frequency amplitude/energy
z2(2:K/2) = 2*z2(2:K/2);  % Accepting and doubling half of the frequencies
z2(K/2+1:K) = 0;          % Rejecting half of the frequencies
z2 = ifft(z2);            % Transforming Z into the time domain

%% Computing the Signal Kernel in the Doppler-Lag Domain
N = length(s1);
K_TL = zeros(K, N);
L_half = fix(L/2);
for n = 1:K
    for tau = -L_half:L_half
        Z1 = z1(1 + rem(2*K + n - 1 + tau, K));
        Z2 = z2(1 + rem(2*K + n - 1 - tau, K));
        mm = 1 + rem(2*K + tau, M);
        K_TL(mm,n) = Z1*conj(Z2);
    end
end
dump = zeros(M, N);
temp = fft(K_TL,M);
dump(:,2:N) = temp(:,1:N-1);
AD = fftshift(ifft(fft(fftshift(dump),[],1),[],2));

%% Computing the Cross Compact Support Kernel Distribution
GAD = AD.*G;
TFR = fftshift(ifft(fft(fftshift(GAD),[],2),[],1));
