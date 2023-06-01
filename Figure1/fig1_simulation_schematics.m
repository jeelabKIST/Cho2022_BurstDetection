%% Configure Library Path
util_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/utils');
data_path = genpath('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection/data/simulation_data');
addpath(util_path);
addpath(data_path);
%% Set Colors
red = '#d13b3b';
orange = '#d98100';
bluel = '#377A98';
blued = '#244A80';
blue_shade = '#d7ebfa';
black = '#242621';
%% Simulate Example Signal
% [1] Set Parameters
method = 'tukey';
T = 10;
f = 9;
f_spans = [1,2,4];
nBurst = 8;
noise_snr = 3;
rng_pks = [3 12];
rnd_seed = 11;
sim_verbose = false;
% [2] Simulate a Signal
[tvec,simulated_signal,Fs,data,btime,burst,bfreq,cyc] = generate_simulation(method,T,f,f_spans,nBurst,noise_snr,rng_pks,rnd_seed,sim_verbose);
%% Plot Simulated Signal Schematics
fnt_sz = 54;
draw_schematics(tvec,data(:,2),{'θ Bursts';'(8-10Hz)'},fnt_sz,bluel);
draw_schematics(tvec,data(:,3),{'β Bursts';'(23-27Hz)'},fnt_sz,blued);
draw_schematics(tvec,data(:,4),{'γ Bursts';'(37-45Hz)'},fnt_sz,blued);
draw_schematics(tvec,data(:,5),{'Aperiodic';'Noise'},fnt_sz,orange);
draw_schematics(tvec,data(:,6),{'Gaussian';'Noise'},fnt_sz,red);
draw_schematics(tvec,simulated_signal,{'Simulated';'Signal'},fnt_sz,black,true);
%% Plot Welch's Power Spectral Density Estimate (Entire Dataset)
method = 'tukey';
T = 300;           % unit: s
f = 9;             % unit: Hz
f_spans = [1,2,4];
nBurst = 30;
rng_pks = [3 12];
sim_verbose = false;
% [1] Get Random Parameters
nSimulationPerNoise = 1000; % number of simulations per each SNR
listNoise = -10:2:10;
if isfile('randseed.mat')
    listSeed = load('randseed.mat').randseed;
else
    listSeed = zeros(length(listNoise),nSimulationPerNoise);
    for n = 1:length(listNoise)
        listSeed(n,:) = randsample(1:2^15,nSimulationPerNoise);
    end
end
nSeed = numel(listSeed);
% [2] Simulate Signals
fprintf("Starting simulation ... \n");
[psd_sig_data, psd_noise_data] = deal(cell(size(listSeed)));
for n = 1:length(listNoise)
    SNR_wn = listNoise(n);
    for s = 1:size(listSeed,2)
        [tvec,simulated_signal,Fs,data,~,~,~,~] = generate_simulation(method,T,f,f_spans,nBurst,SNR_wn,rng_pks,listSeed(s),sim_verbose);
        psd_window = Fs*30;
        noverlap = floor(psd_window/2); % 50 percent overlap
        nfft = psd_window;
        [pxx_sig,f_sig] = pwelch(simulated_signal,hanning(psd_window),noverlap,nfft,Fs);
        [pxx_noise,f_noise] = pwelch(data(:,end-1)+data(:,end),hanning(psd_window),noverlap,nfft,Fs);
        psd_sig_data{n,s} = pxx_sig;
        psd_noise_data{n,s} = pxx_noise;
    end
    fprintf([num2str(n) '/' num2str(length(listNoise)) ' complete. \n']);
end
psd_sig = sum([psd_sig_data{:}],2)./nSeed;
psd_noise = sum([psd_noise_data{:}],2)./nSeed;
% [3] Define y-axis Boundary
f_idx = find(f_sig == 60);
vmin = min(psd_sig(1:f_idx),[],'all')-1e-4;
vmax = max(psd_sig(1:f_idx),[],'all')+1e-3;
% [4] Plot Welch's PSD Estimate
figure(); hold on;
plot(f_sig,psd_sig,'Linewidth',5,'Color','#242621');
plot(f_noise,psd_noise,'Linewidth',5,'Color','#d98100');
xlim([-0.5 60]);
ylim([vmin,vmax]);
xlabel('Frequency (Hz)');
ylabel('PSD (mV^2/Hz)');
xticks(0:10:60);
lgnd = legend('Simulated Signal','Simulated Noise');
lgnd.FontSize = 35;
set(gca,'YScale','log','LineWidth',6,'TickDir','out','FontSize',45,'FontName','Helvetica','FontWeight','bold');
set(gcf,'Units','Normalized','Color','w','Position',[0.205729166666667,0,0.6546875,0.9037]);
%% Plot Welch's Power Spectral Density Estimate (Different SNRs)
% [1] Compute Welch's PSD Estimate with Different SNRs
psd_sig_snr = cell(1,length(listNoise));
seed = 1997;
fprintf("Starting simulation ... \n");
for n = 1:length(listNoise)
    [tvec,simulated_signal,Fs,data,~,~,~,~] = generate_simulation(method,T,f,f_spans,nBurst,listNoise(n),rng_pks,seed,sim_verbose);
    psd_window = Fs*10;
    noverlap = floor(psd_window/2); % 50 percent overlap
    nfft = psd_window;
    [pxx_sig,f_sig] = pwelch(simulated_signal,hanning(psd_window),noverlap,nfft,Fs);
    psd_sig_snr{n} = pxx_sig;
    fprintf([num2str(n) '/' num2str(length(listNoise)) ' complete. \n']);
end
% [2] Define y-axis Boundary
combined_psd = [psd_sig_snr{:}];
f_idx = find(f_sig == 60);
vmin = min(combined_psd(1:f_idx,:),[],'all')-5e-5;
vmax = 0.1;
% [3] Plot Welch's PSD
figure(); hold on;
xf = [9, 25, 41];
for i = 1:length(f_spans)
    rectangle('Position',[xf(i)-f_spans(i) vmin+1e-5 2*f_spans(i) 1],'EdgeColor',blue_shade,'FaceColor',blue_shade);
    text(xf(i),vmax,[num2str(xf(i)-f_spans(i)) '-' num2str(xf(i)+f_spans(i)) ' Hz'],'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontSize',35,'Color',blued);
end
colororder(hot(length(listNoise)+5));
[p, leg_lbl] = deal(cell(1,length(listNoise)));
for n = 1:length(listNoise)
    p{n} = plot(f_sig,psd_sig_snr{n},'LineWidth',3);
    leg_lbl{n} = [num2str(listNoise(n)) ' dB'];
end
xlim([-0.5 60]);
ylim([vmin,vmax]);
xlabel('Frequency (Hz)');
ylabel('PSD (mV^2/Hz)');
xticks(0:10:60);
% (A) Set Legend
leg = legend([p{:}],leg_lbl);
title(leg,'SNR Levels');
% (B) Figure Settings
set(leg,'FontSize',26,'Position',[0.8339,0.2107,0.1515,0.4790]);
set(gca,'YScale','log','LineWidth',6,'TickDir','out','FontName','Helvetica','FontSize',45,'FontWeight','bold','Position',[0.2109,0.2106,0.6014,0.7500]);
set(gcf,'Units','Normalized','Color','w','Position',[0.205729166666667,0,0.6546875,0.9037]);
%% Appendix: In-Script Functions
% Function #1: Plot Simulation Schematics
function draw_schematics(tvec,signal,label,fnt_sz,clr,scale_opt)
    if nargin < 6
        scale_opt = false;
    end
    figure(); hold on;
    plot(tvec,signal,'Color',clr,'LineWidth',4);
    text(3.52,0,label,'FontSize',fnt_sz,'FontWeight','bold','HorizontalAlignment','center');
    ylim([-2.7,2.8]);
    gcf_position = [0.07473544973545,0.395112016293279,0.830687830687831,0.343177189409369];
    if scale_opt
        line([6.0,6.5],[-2.6,-2.6],'LineWidth',6.0,'Color',clr);
        line([6.5,6.5],[-2.6,-1.6],'LineWidth',6.0,'Color',clr);
        text(6.52,-2.1,'1mV','FontSize',fnt_sz,'FontName','Helvetica','FontWeight','bold');
        text(6.01,-3.3,'500ms','FontSize',fnt_sz,'FontName','Helvetica','FontWeight','bold');
    end
    set(gca,'TickLength',[0 0],'XColor','none','YColor','none','xlim',[4.0, 6.5],'Color','none', ... 
        'Position',[0.244411692458855,0.228486646884273,0.64756606766823,0.638377258441171]);
    set(gcf,'Units','Normalized','Position',gcf_position,'Color','white');
end
