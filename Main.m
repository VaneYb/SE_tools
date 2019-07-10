clc
clear
close all

filename = 'meetingroom_snr0.wav';

[data,fs] = audioread(filename);
if(fs == 16000)
    y = data;
else
    y = resample(data,16000,fs);
    fs = 16000;
end
clear data

x = do_proc_in_freq(single(y(:,1)));
% plot((y(:,1)));hold on;
% plot(x);% ;
% audiowrite(['out_' filename],x,16000)
audiowrite(['out.wav'],x,fs)
% plot(y(:,1)); hold on
% plot(x)

