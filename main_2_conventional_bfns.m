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
% This main script generates the conventional connectivity measures for 
% each segment and each frequency band to build a tempo-spectral FBN for 
% each subject in the dataset. The generated brain networks are saved in 
% "Brain Networks\Conventional".

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Main script
input_path  = fullfile(pwd,"Processed EEG");
output_path = fullfile(pwd,"Brain Networks","Conventional");
if ~isfolder(output_path), mkdir(output_path); end
logfilename = 'log_brain_network_generation.txt';
if exist(logfilename,'file'), delete(logfilename); end
file_logging(logfilename, 'Logging Starts');
M = length(dir(input_path + "/*.mat"));
for m = 1:M % iterate through the patients
    % Signal loading
    file_logging(logfilename, "Processing subject " + string(m) + '/' + string(M));
    load(fullfile(input_path,"eeg" + m + ".mat"));
    sig = eeg_sig;
    % Brain network computation [segment x frequency x channel x channel]
    [Nt, Ns, Nc, Nf] = size(sig);
    pear_network  = zeros(Ns,Nf,Nc,Nc);
    spear_network = zeros(Ns,Nf,Nc,Nc);
    plv_network   = zeros(Ns,Nf,Nc,Nc);
    iplv_network  = zeros(Ns,Nf,Nc,Nc);
    ciplv_network = zeros(Ns,Nf,Nc,Nc);
    pli_network   = zeros(Ns,Nf,Nc,Nc);
    wpli_network  = zeros(Ns,Nf,Nc,Nc);
    dwpli_network = zeros(Ns,Nf,Nc,Nc);
    for f = 1:Nf % iterate through all frequency bands
        file_logging(logfilename, "Frequency band " + string(f) + "/" + string(Nf));
        for n = 1:Ns % iterate through all segments
            [pear_network(n,f,:,:), spear_network(n,f,:,:)] = corr_net(squeeze(sig(:,n,:,f)), 0.05);
            [plv_network(n,f,:,:),  iplv_network(n,f,:,:), ciplv_network(n,f,:,:)] = plv_net(squeeze(sig(:,n,:,f)));
            [pli_network(n,f,:,:),  wpli_network(n,f,:,:), dwpli_network(n,f,:,:)] = pli_net(squeeze(sig(:,n,:,f)));
        end
    end
    % Saving
    file_logging(logfilename, "Saving");
    save(fullfile(output_path,"network" + m + ".mat"), ...
        'eeg_mask','pear_network','spear_network', ...
        'plv_network','iplv_network','ciplv_network', ...
        'pli_network','wpli_network','dwpli_network');
end
