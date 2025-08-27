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
% This main script generates the proposed time-frequency connectivity 
% measures to build a full tempo-spectral FBN for each subject in the dataset.
% The time-frequency brain networks are generated for each frequency band
% and averaged within each segment to allow associating connectivity with 
% the segment-level seizure annotations. The time-frequency brain networks
% are saved in "Brain Networks\Time Frequency".

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Main script
input_path  = fullfile(pwd,"Processed EEG");
output_path = fullfile(pwd,"Brain Networks","Time Frequency");
if ~isfolder(output_path), mkdir(output_path); end
logfilename = 'log_tf_brain_network_generation.txt';
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
        temp = squeeze(sig(:,:,:,f));
        for n = 1:Ns % iterate through all segments
            tf_pear  = nan(Nc,Nc);
            tf_spear = nan(Nc,Nc);
            tf_plv   = nan(Nc,Nc);
            tf_iplv  = nan(Nc,Nc);
            tf_ciplv = nan(Nc,Nc);
            tf_pli   = nan(Nc,Nc);
            tf_wpli  = nan(Nc,Nc);
            tf_dwpli = nan(Nc,Nc);
            for p = 1:Nc
                x1 = squeeze(temp(:,n,p));
                aTFD1 = real(Xckd(x1, x1, 1, 0.45, 0.05))';
                for q = (p+1):Nc
                    x2 = squeeze(temp(:,n,q));
                    aTFD2 = real(Xckd(x2, x2, 1, 0.45, 0.05))';
                    xTFD = Xckd(x1, x2, 1, 0.45, 0.05)';
                    [tf_pear(p,q), tf_spear(p,q)] = tf_corr_net(aTFD1, aTFD2, 0.05);
                    [tf_plv(p,q), tf_iplv(p,q), tf_ciplv(p,q)] = tf_plv_net(xTFD);
                    [tf_pli(p,q), tf_wpli(p,q), tf_dwpli(p,q)] = tf_pli_net(xTFD);
                    tf_pear(q,p)  = tf_pear(p,q);
                    tf_spear(q,p) = tf_spear(p,q);
                    tf_plv(q,p)   = tf_plv(p,q);
                    tf_iplv(q,p)  = tf_iplv(p,q);
                    tf_ciplv(q,p) = tf_ciplv(p,q);
                    tf_pli(q,p)   = tf_pli(p,q);
                    tf_wpli(q,p)  = tf_wpli(p,q);
                    tf_dwpli(q,p) = tf_dwpli(p,q);
                end
            end
            pear_network(n,f,:,:)  = tf_pear;
            spear_network(n,f,:,:) = tf_spear;
            plv_network(n,f,:,:)   = tf_plv;
            iplv_network(n,f,:,:)  = tf_iplv;
            ciplv_network(n,f,:,:) = tf_ciplv;
            pli_network(n,f,:,:)   = tf_pli;
            wpli_network(n,f,:,:)  = tf_wpli;
            dwpli_network(n,f,:,:) = tf_dwpli;
        end
    end
    % Saving
    file_logging(logfilename, "Saving network");
    save(fullfile(output_path,"network" + m + ".mat"), ...
        'eeg_mask','pear_network','spear_network', ...
        'plv_network','iplv_network','ciplv_network', ...
        'pli_network','wpli_network','dwpli_network');
end
