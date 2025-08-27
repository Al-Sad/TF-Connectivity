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
% This main script preprocesses the newborn EEG signals via filtering,
% downsampling, frequency-band decomposition, and time segmentation. The
% preprocessed signals are saved in the "Processed EEG" folder.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));
dataset_path = fullfile(pwd,"The Helsinki EEG dataset"); % This is the raw dataset path

%% Parameters
fs_old = 256;  % The old sampling rate
fs_new = 64;   % The new sampling rate
seg_T  = 16;   % Segment length in seconds
seg_P  = 4;    % Segment overlap length in seconds
seiz_L = 8;    % Decision aggregation length

%% EEG IIR bandpass filters
full_h = designfilt('bandpassiir','PassbandFrequency1',0.5,'StopbandFrequency1',0.4,...
    'PassbandFrequency2',30,'StopbandFrequency2',31.5,'StopbandAttenuation1',30,...
    'PassbandRipple',1,'StopbandAttenuation2',30,'match','stopband','SampleRate',fs_old);
delta_h = designfilt('bandpassiir','PassbandFrequency1',0.5,'StopbandFrequency1',0.4,...
    'PassbandFrequency2',4,'StopbandFrequency2',4.3,'StopbandAttenuation1',30,...
    'PassbandRipple',1,'StopbandAttenuation2',30,'match','stopband','SampleRate',fs_new);
theta_h = designfilt('bandpassiir','PassbandFrequency1',4,'StopbandFrequency1',3.8,...
    'PassbandFrequency2',8,'StopbandFrequency2',8.2,'StopbandAttenuation1',30,...
    'PassbandRipple',1,'StopbandAttenuation2',30,'match','stopband','SampleRate',fs_new);
alpha_h = designfilt('bandpassiir','PassbandFrequency1',8,'StopbandFrequency1',7.7,...
    'PassbandFrequency2',13,'StopbandFrequency2',13.3,'StopbandAttenuation1',30,...
    'PassbandRipple',1,'StopbandAttenuation2',30,'match','stopband','SampleRate',fs_new);
beta_h = designfilt('bandpassiir','PassbandFrequency1',13,'StopbandFrequency1',12.7,...
    'PassbandFrequency2',30,'StopbandFrequency2',30.3,'StopbandAttenuation1',30,...
    'PassbandRipple',1,'StopbandAttenuation2',30,'match','stopband','SampleRate',fs_new);
eeg_h = {delta_h, theta_h, alpha_h, beta_h};

%% Pre-processing the EEG signals
logfilename = 'log_preprocessing.txt';
if exist(logfilename,'file'), delete(logfilename); end
file_logging(logfilename, 'Logging Starts');
load(fullfile(dataset_path,'annotations_2017.mat'));
M = length(dir(dataset_path + "/*.edf")); % number of patients
if ~isfolder("Processed EEG"), mkdir("Processed EEG"); end
for m = 2:M % iterate through the patients
    file_logging(logfilename, "Processing subject " + string(m) + '/' + string(M));
    % Get the seizure masks
    mask = sum(annotat_new{m}) > 1;
    % Load the EEG signal
    file_logging(logfilename, "Data Loading");
    data = edfread(fullfile(dataset_path,['eeg' num2str(m) '.edf']));
    % Take the EEG sensors only
    x_raw = table2array(data); x_raw = x_raw(:,1:19);
    x_raw = cell2mat(x_raw);
    % Signal bandpass filtering
    file_logging(logfilename, "Bandpass filtering");
    signal_ud = flipud(x_raw);
    signal_big = [signal_ud; x_raw; signal_ud];
    temp = filtfilt(full_h, signal_big);
    temp = temp(size(x_raw,1)+1:2*size(x_raw,1),:);
    % Resampling
    file_logging(logfilename, "Frequency down-sampling");
    temp = resample(temp,fs_new,fs_old);
    % EEG frequency band decomposition
    dump = zeros([size(temp), length(eeg_h)]);
    for i = 1:length(eeg_h)
        file_logging(logfilename, "Frequency band decomposition " + ...
            string(i) + "/" + string(length(eeg_h)));
        signal_ud = flipud(temp);
        dump1 = filtfilt(eeg_h{i}, [signal_ud; temp; signal_ud]);
        dump(:,:,i) = dump1(size(temp,1)+1:2*size(temp,1),:);
    end
    temp = dump; clear dump; clear dump1;
    % Segmentation
    file_logging(logfilename, "Signal Segmentation");
    eeg_sig  = zeros([seg_T*fs_new ceil((size(temp,1)-seg_P*fs_new)/(seg_T*fs_new-seg_P*fs_new)) size(temp,[2 3])]);
    eeg_mask = sum(buffer(mask,seg_T,seg_P,'nodelay')) >= seiz_L;
    for i = 1:size(temp,2)
        for j = 1:size(temp,3)
            eeg_sig(:,:,i,j) = buffer(squeeze(temp(:,i,j)),seg_T*fs_new,seg_P*fs_new,'nodelay');
        end
    end
    % Saving
    file_logging(logfilename, "Saving");
    save(fullfile("Processed EEG","eeg" + m + ".mat"),"eeg_sig","eeg_mask");
end
