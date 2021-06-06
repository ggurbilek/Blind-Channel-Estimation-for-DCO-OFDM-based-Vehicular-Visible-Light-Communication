clear all
%close all
clc

load('sumodata.mat')

figure
plot(vtime,dist)
xlabel('Time (s)')
ylabel('Distance (m)')
title('SUMO inter-vehicle distance in time')

t=0:0.0000005:vtime(end);
d=interp1(vtime,dist,t);

const=[4 8 16 32 64];
    
%% Local Variables: number of bits, samples, guard symbols, temp variables for interfunction calls
nfft=64; %fft size
subcar=48; %number of data subcarriers
guard = (nfft-subcar)/2;
nsym=10000; %number of ofdm symbols

ncp=4; %cyclic prefix size
nbitsym = nfft + ncp;


k=1:floor(length(t)/(nbitsym*nsym));

% Channel Effect 1
tt=t((nbitsym*nsym)*(k-1)+1:(k*nbitsym*nsym));
tt0(k)=t((nbitsym*nsym)*(k-1)+1);


        
load('sumober10.mat')

%% Plot SNR vs. BER for block type CE
figure
semilogy(qsnr(1,:),qber(1,:),'d')
hold on
semilogy(qsnr(2,:),qber(2,:),'s')
hold on
semilogy(qsnr(3,:),qber(3,:),'p')
hold on
semilogy(qsnr(4,:),qber(4,:),'h')
hold on
semilogy(qsnr(5,:),qber(5,:),'^')
xlabel('Eb/N0 (dB)')
ylabel('BER')
legend('4-QAM','8-QAM','16-QAM','32-QAM','64-QAM')
title('Proposed CE SUMO')
grid on

figure
semilogy(tt0,(1-qber(1,:))*(2*subcar/2)/(nfft+ncp),'-m','LineWidth',2)
hold on
semilogy(tt0,(1-qber(2,:))*(3*subcar/2)/(nfft+ncp),'-g','LineWidth',2)
hold on
semilogy(tt0,(1-qber(3,:))*(4*subcar/2)/(nfft+ncp),'-b','LineWidth',2)
hold on
semilogy(tt0,(1-qber(4,:))*(5*subcar/2)/(nfft+ncp),'-c','LineWidth',2)
hold on
semilogy(tt0,(1-qber(5,:))*(6*subcar/2)/(nfft+ncp),'-r','LineWidth',2)
xlabel('Time (s)')
ylabel('Spectral Efficiency (C [bits/s/Hz])')
legend('4-QAM','8-QAM','16-QAM','32-QAM','64-QAM')
%title('Proposed CE SUMO')
%ylim([0 6.5])
grid on




