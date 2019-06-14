clear 
close all

[y,fs] = audioread('bgnoise_snr5_p20_02.wav');
[~, Nch] = size(y);
% [B,A] = adsgn(fs);    %A-weighting filter 
% y = filter(B,A,y);    

%% %%%% Manual Setting %%%%%
% You need to set STARTPOINT and NOISESEG by yourself first.
% Turn on DRAW_PLOT to help you decide the values.

draw_plot = 1;
startpoint = 1;
noiseSeg = 1:1e4;


% other parameters
noiseThres = 1.4;
segLen  = 512;

%% Initialize valid signal segment

if draw_plot
    figure(1)
    subplot(3,1,1)
    u = zeros(length(y),1);
    u(startpoint:end) = 0.9*max(max(abs(y)));
    if startpoint==1
        u(1) = 0;
    end
    plot(y(:,1)); grid on; hold on 
    plot(u,'LineWidth',1.5); axis tight;
    title('Check STARTPOINT to pick up valid signal segment. ')
end

%% convert signal sampling frequency to 16kHz 
if fs ~= 16000
    warning(['Input data is ', num2str(fs), 'Hz. Convert to 16000Hz.' ]);
    y = resample(y(startpoint:end,:),16000,fs);
else
    y = y(startpoint:end,:);
end

%% initialize noise power
if draw_plot
    subplot(3,1,2)
    u = zeros(length(y),1);
    u(noiseSeg) = 0.9*max(max(abs(y)));
    plot(y(:,1)); grid on; hold on 
    plot(u,'LineWidth',2); axis tight %axis([1,1e5,-0.006,0.016])
end
pow_n = diag(y(noiseSeg,:)'*y(noiseSeg,:))'/length(noiseSeg);
title('Check NOISESEG to pick up valid noise segment. ')

%% calculate SNR
yLen    = length(y);
Nframe  = floor(yLen/segLen);
yPow = zeros(Nframe,Nch);
for n = 1:Nframe
    yPow(n,:) = diag(y((n-1)*segLen+1:n*segLen,:)'*y((n-1)*segLen+1:n*segLen,:))/segLen;
end


for i=1:Nch
    ind = find(yPow(:,i)>pow_n(i)*noiseThres);
    sig_pow(i) = sum(yPow(ind,i))/length(ind);
    sig_powdB(i) = db(sig_pow(i), 'power');
    ind = find(yPow(:,i)<=pow_n(i)*noiseThres);
    noi_pow(i) = sum(yPow(ind,i))/length(ind);
    noi_powdB(i) = db(noi_pow(i), 'power');
    
    snr(i) = 10*log10(sig_pow(i)/noi_pow(i));
end
%print result
sig_powdB
noi_powdB
snr   

A = noiseThres;
if draw_plot
    subplot(3,1,3)
    pow_sort = sort(yPow(:,1));
    plot(pow_sort,'LineWidth',1);hold on
    u = zeros(size(pow_sort));
    ind = find(pow_sort<A*pow_n(1));
    u(ind) = max(yPow(:,1))/10;
    plot(u,'LineWidth',1);axis tight
    title('Sorted power denotes noise power and signal power')
end