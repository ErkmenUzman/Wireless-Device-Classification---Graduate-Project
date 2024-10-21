clear;
clc;
close all;

%-------------------------------------------------------------------------
% SİMÜLASYON AYARLARI
% Initial Parameters

tic

add_awgn = true;                    % AWGN eklensin mi
add_iq_imbalance = true;            % IQ Imbalance Eklensin mi
OverralCell = false;
IQcell = false;
SampleCell = true;
PacketCell = false;
RawData = false;
PlotAll = false;

iq_amp_values = [0.0864 0.3406 0.5876 0.8279 1.0615];
iq_phase_values = [1 4 7 10 13];

monte_carlo_size = length(iq_amp_values);   % Her bir SNR seviyesesnde simülasyonun kaç kez yapılacağını belirler

iq_snr_level = 20;              % Simülasyonda kullanılacak Monte Conte SNR seviyeleri (dB)
repeat_times = 8000;
% iq_amp_level_interval = [10 3];   genlik dengesizliği alt ve üst sınır (dB)
% iq_phase_level_interval = [8 2];    faz dengesizliği alt ve üst sınır (derece)

rxMatrixCell = cell(monte_carlo_size,1);
rxMatrixCellIQ = cell(monte_carlo_size,2);

iq_amp_values_matrix =   NaN(length(iq_snr_level),monte_carlo_size);    % grafik için tüm genlik dengesizliği değerlerini tutar
iq_phase_values_matrix = NaN(length(iq_snr_level),monte_carlo_size);    % grafik için tüm faz dengesizliği değerlerini tutar
%-------------------------------------------------------------------------

%%

% TRANSMITTER
% -------------------------------------------------------------------------

% BURLARA DOKUNMA

%Configure all the scopes and figures for the example.

% Setup handle for image plot
if ~exist('imFig','var') || ~ishandle(imFig) %#ok<SUSENS>
    imFig = figure;
    imFig.NumberTitle = 'off';
    imFig.Name = 'Image Plot';
    imFig.Visible = 'off';
else
    clf(imFig); % Clear figure
    imFig.Visible = 'off';
end

% Setup Spectrum viewer
spectrumScope = dsp.SpectrumAnalyzer( ...
    'SpectrumType','Power density', ...
    'SpectralAverages',10, ...
    'YLimits',[-90 -30], ...
    'Title','Received Baseband WLAN Signal Spectrum', ...
    'YLabel','Power spectral density', ...
    'Position',[69 376 800 450]);

% Setup the constellation diagram viewer for equalized WLAN symbols
refQAM = wlanReferenceSymbols('64QAM');
constellation = comm.ConstellationDiagram(...
    Title='Equalized WLAN Symbols',...
    ShowReferenceConstellation=true,...
    ReferenceConstellation=refQAM,...
    Position=[878 376 460 460]);

%%
% Prepare Image File

% Input an image file and convert to binary stream
fileTx = 'peppers.png';                          % Image file name
fData = imread(fileTx);                          % Read image data from file
scale = 0.3;                                     % Image scaling factor
origSize = size(fData);                          % Original input image size
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calculate new image size
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:);               % Resize image
imsize = size(fData);                            % Store new image size
txImage = fData(:);


% % Plot transmit image
% imFig.Visible = 'on';
% subplot(211);
% imshow(fData);
% title('Transmitted Image');
% subplot(212);
% title('Received image appears here...');
% set(gca,'Visible','off');
%
% set(findall(gca, 'type', 'text'), 'visible', 'on');

%%
%Fragment Transmit Data

msduLength = 2304; % MSDU length in bytes
numMSDUs = ceil(length(txImage)/msduLength);
padZeros = msduLength-mod(length(txImage),msduLength);
txData = [txImage;zeros(padZeros,1)];
txDataBits = double(reshape(de2bi(txData, 8)',[],1));

% Divide input data stream into fragments
bitsPerOctet = 8;
data = zeros(0,1);

for i=0:numMSDUs-1

    % Extract image data (in octets) for each MPDU
    frameBody = txData(i*msduLength+1:msduLength*(i+1),:);

    % Create MAC frame configuration object and configure sequence number
    cfgMAC = wlanMACFrameConfig(FrameType='Data',SequenceNumber=i);

    % Generate MPDU
    [psdu, lengthMPDU]= wlanMACFrame(frameBody,cfgMAC,OutputFormat='bits');

    % Concatenate PSDUs for waveform generation
    data = [data; psdu]; %#ok<AGROW>

end

%%
% Generate 802.11a Baseband WLAN Signal

nonHTcfg = wlanNonHTConfig;       % Create packet configuration
nonHTcfg.MCS = 6;                 % Modulation: 64QAM Rate: 2/3
nonHTcfg.NumTransmitAntennas = 1; % Number of transmit antenna
chanBW = nonHTcfg.ChannelBandwidth;
nonHTcfg.PSDULength = lengthMPDU; % Set the PSDU length

scramblerInitialization = randi([1 127],numMSDUs,1);

osf = 1.5;

sampleRate = wlanSampleRate(nonHTcfg); % Nominal sample rate in Hz

% Generate baseband NonHT packets separated by idle time
txWaveform = wlanWaveformGenerator(data,nonHTcfg, ...
    'NumPackets',numMSDUs,'IdleTime',20e-6, ...
    'ScramblerInitialization',scramblerInitialization, ...
    'OversamplingFactor',osf);


%%
% RECEIVER


fs = sampleRate*osf ;
CarrierFrequency = 2.4*10^9;

% Choose tgn_channel and define the tgn channel which you will work on
tgn_channel = wlanTGnChannel('SampleRate', fs,'LargeScaleFadingEffect','Pathloss and shadowing','CarrierFrequency',CarrierFrequency,'DelayProfile','Model-A');

process = 0; % simülasyonun yüzde kaçının tamamlandığını göstermek için gösterge

for monte_carlo_iter = 1:monte_carlo_size %This nested loop creates your dataset, size of the dataset can be adjusted by changing monte_carlo_size and iq_snr_iter variables

    iq_amp_values   = sort(iq_amp_values);
    iq_phase_values = sort(iq_phase_values);

    for iq_snr_iter = 1:repeat_times

        process = process + 1;  % simülasyonun yüzde kaçının tamamlandığını göstermek için gösterge
        is_successfull = true;

        % Add AWGN and IQ Imbalance
        if add_awgn == true && add_iq_imbalance == false
            rxWaveform = awgn(txWaveform,iq_snr_level(iq_snr_iter),'measured');    % Sinyal awgn kanalından geçirilir
        end

        if add_awgn == false && add_iq_imbalance == true
            rxWaveform = iqimbal(txWaveform, iq_amp_values(monte_carlo_iter),iq_phase_values(monte_carlo_iter));   % girilen değerlerde IQ Imbalance Eklenir
        end
            
        if add_awgn == true && add_iq_imbalance == true %tx tarafında iq imbalance eklenmesi için iki fonksiyonun yerini değiştirmek yeterli gibi ama awgn için kod sıkıntı verdiğinden doğrulayamadım.
            rxWaveform_temp = iqimbal(txWaveform, iq_amp_values(monte_carlo_iter),iq_phase_values(monte_carlo_iter));   % girilen değerlerde IQ Imbalance Eklenir
            rxWaveform_awgn = awgn(rxWaveform_temp,iq_snr_level,'measured');    % Sinyal awgn kanalından geçirilir
            rxWaveform = tgn_channel(rxWaveform_awgn);
        end

        if add_awgn == false && add_iq_imbalance == false
            rxWaveform = txWaveform;
        end
        
        if OverralCell == true
            rxMatrixCell{monte_carlo_iter} = num2cell(rxWaveform);
        end

        if IQcell == true
            rxMatrixCellIQ{monte_carlo_iter , 1} = num2cell(real(rxWaveform));
            rxMatrixCellIQ{monte_carlo_iter , 2} = num2cell(imag(rxWaveform));
        end

        if SampleCell == true

            % Convert rxWaveform to arrays
            realPart = real(rxWaveform(1:osf*160));
            imagPart = imag(rxWaveform(1:osf*160));
    
            % Assign the arrays to rxMatrixIQ
            rxMatrixRawIQ(:, 1) = realPart;
            rxMatrixRawIQ(:, 2) = imagPart;

            % Generate a unique filename based on monte_carlo_iter
            filename = sprintf('rxMatrixArrayIQ_Batch2_Class_%d_Sample_%d.mat', monte_carlo_iter , iq_snr_iter);
    
            % Save rxMatrixIQ to the unique filename
            save(filename, 'rxMatrixRawIQ');

        end

        if RawData == true

            % Convert rxWaveform to arrays
            realPart = real(rxWaveform);
            imagPart = imag(rxWaveform);
    
            % Assign the arrays to rxMatrixIQ
            rxMatrixRawIQ(:, 1) = realPart;
            rxMatrixRawIQ(:, 2) = imagPart;

            % Generate a unique filename based on monte_carlo_iter
            filename = sprintf('rxMatrixArrayIQ_Class_%d_Sample_%d.mat', monte_carlo_iter , iq_snr_iter);
    
            % Save rxMatrixIQ to the unique filename
            save(filename, 'rxMatrixRawIQ');

        end
        
        if PlotAll == true
            figure;
            plot(abs(rxWaveform(1:160*osf)));
            hold on;
        end

    end
    
    if OverralCell == true
        save('rxMatrixCell.mat', 'rxMatrixCell');
    end

    if IQcell == true
        save('rxMatrix_IQ.mat', 'rxMatrixCellIQ');
    end

    if SampleCell == true
        save('rxMatrix_Sampled.mat', 'rxMatrixCellIQ');
    end
    
end

toc