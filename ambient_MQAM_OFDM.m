%% DCO-OFDM, Only Half Data Subcarriers are used for QAM Data for Hermitian Symmetry
clear all
clc
%close all

%% Local Variables: number of bits, samples, guard symbols, temp variables for interfunction calls
qam=16;
nfft=64; %fft size
subcar=48; %number of data subcarriers
nbitqam = log2(qam);
guard = (nfft-subcar)/2;
nsym=600*nbitqam; %number of ofdm symbols
nbitdata = subcar*nbitqam*nsym/2;

ncp=4; %cyclic prefix size
nbitsym = nfft + ncp;

data=(rand(1,nbitdata)>0.5)+0;

%channnel effect: AWGN at different SNRs
% alpha is from F_2M_20dB
alpha=[0.609981229934408,0.715459955559726,0.636965602683617,0.694127681700978,0.695098676789100,0.681699620219806,0.705628746895461,0.834320339106743,0.873378196662836,0.827668204581301,0.768804769140058,0.756148741911805,0.926978115212095,0.913129736817394,0.871793084001999,0.848536621661314,0.909045348037926,0.834381831138879,0.896913346571862,0.928182910840626,0.932323910122903,0.907741812468878,0.943650044868540,1];

% set PS(kb) for block type, PS(kc) for comb type
kb=150;
kc=4;

symbols=1:nsym*subcar/2;
bits=1:nsym*nbitqam*subcar/2;
Cbits=1:nsym*nbitqam*subcar/2*(1/2);
CCbits=1:nsym*nbitqam*subcar/2*(2/3);
CCCbits=1:nsym*nbitqam*subcar/2*(3/4);
CCCCbits=1:nsym*nbitqam*subcar/2*(kc/(kc+1));
Bbits=1:nsym*nbitqam*subcar/2*(1/2);
BBbits=1:nsym*nbitqam*subcar/2*(2/3);
BBBbits=1:nsym*nbitqam*subcar/2*(3/4);
BBBBbits=1:nsym*nbitqam*subcar/2*(kb/(kb+1));

%theta = logspace(-1,4,100); 

% theta for dynamic case
d0=1;
%theta_0=0.3513;
theta_0_dB = [-7.06432267663874;-8.70770993276061;-11.5993217774338;-13.9108884331651;-26.2641929929968;-25.9573144830053;-22.2466917438169;-12.7328959588067;-11.0212434828756;-15.3158572581427].'+40;% 40dB PD gain
theta_0=db2mag(theta_0_dB);
ambient_light={'a','b','c','d','e','f','g','h','i','j'};
v=30; %m/s2
a=0; %m/s2
dt=1/(2*1e6);

% set seed for random number generation
%rng default

for n=1:10%n=1:10
    rng(n);
    noise1=normrnd(0,1,[1,subcar*nsym/2])+i*normrnd(0,1,[1,subcar*nsym/2]);
    rng(n);
    noise2=normrnd(0,1,[1,subcar*nsym])+i*normrnd(0,1,[1,subcar*nsym]);
    rng(n);
    noise3=normrnd(0,1,[1,subcar*nsym*0.75])+i*normrnd(0,1,[1,subcar*nsym*0.75]);
    rng(n);
    noise4=normrnd(0,1,[1,subcar*nsym/2*(1+1/3)])+i*normrnd(0,1,[1,subcar*nsym/2*(1+1/3)]);
    rng(n);
    noisekb=normrnd(0,1,[1,round(subcar*nsym*(1+1/kb)/2)])+i*normrnd(0,1,[1,round(subcar*nsym*(1+1/kb)/2)]);
    rng(n);
    noisekc=normrnd(0,1,[1,subcar*nsym/2*(1+1/kc)])+i*normrnd(0,1,[1,subcar*nsym/2*(1+1/kc)]);
    
    for k=1:length(theta_0)
        %% THIS IS FOR THE PROPOSED ALGORITHM
        % TX
        % parallalize bin data and do qam modulation
        qam_data_tx = qammod(reshape(data, nbitqam,[]),qam,'InputType','bit','UnitAveragePower',true);
        
        % parallalize qam data for collective ifft
        par_data = reshape(qam_data_tx.', subcar/2, nsym).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        pilot_ins_data=[zeros(nsym,guard) par_data zeros(nsym,1) conj(flip(par_data.').') zeros(nsym,guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        ifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(pilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        cyclic_add_data = [ifft_data(:,((nfft - ncp +1):nfft)) ifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        ofdm_tx = reshape(cyclic_add_data.',nbitsym*nsym,1);

        %%
        % Channel Effect 1
        t=(0:(nbitsym*nsym-1))*dt;
        theta=(theta_0(k)./(d0+v*t+0.5*a*t.^2).^3.346).';
        ch_data=ofdm_tx.*theta;
                
        
        %% 
        % Receiver        
        %parallalize received data
        par_rec_data = reshape(ch_data.', nbitsym, nsym).';

        % remove cyclic prefix
        cyclic_pre_rem=par_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        fft_data = fftshift(fft(cyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        rem_pilot = fft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        ch_rem_pilot = reshape((rem_pilot.*alpha).',1,[])+noise1;
        
        % calculate SNR
        EsN0(n,k)=mean(abs(ch_rem_pilot(symbols)).^2)/mean(abs(noise1).^2);
        
        % parallalize ch data
        pardata=reshape(ch_rem_pilot, subcar/2,[]).';

        % estimate theta
        qamconst=qammod(0:(qam-1),qam,'UnitAveragePower',true);
        est_theta=sum(((abs(real(pardata))+abs(imag(pardata)))./alpha).')./(mean(abs(real(qamconst))+abs(imag(qamconst)))*subcar/2);

        % calculate theta*data
        td=(pardata./alpha).';

        % estimate qam_data_tx
        qam_data_rx=td./est_theta;

        % serialize qam data
        serdata=reshape(qam_data_rx, subcar*nsym/2,[]).';

        % demodulate the qam data in the ofdm symbols
        bin_data_rx = reshape(qamdemod(serdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        BER(n,k)=sum((data(bits)~=bin_data_rx(bits)))/length(bits);
        
        %% THIS IS FOR ST-BASED CE ALGORITHM (# of data symbols = # of pilot symbols)
        rng(n);
        % TX
        % optimal relative power of pilot symbols compared to total symbols
        % (qam symbol power is 1)
        Salpha = 0.748;
        
        % number of averaged OFDM symbols
        SM = 200;
        
        % parallalize bin data and do qam modulation
        Sqam_data_tx = reshape(qam_data_tx,subcar/2,[]); 
        
        % parallalize pilot symbols
        Spilot_sym = (1+i)*ones(subcar/2,nsym);
        
        % superimpose data and pilot symbols w.r.t. their relative power
        Spar_data = (Sqam_data_tx/mean(mean(abs(Sqam_data_tx).^2)) + sqrt(Salpha/(1-Salpha))*Spilot_sym/sqrt(mean(mean(abs(Spilot_sym).^2)))).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        Spilot_ins_data=[zeros(nsym,guard) Spar_data zeros(nsym,1) conj(flip(Spar_data.').') zeros(nsym,guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        Sifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(Spilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        Scyclic_add_data = [Sifft_data(:,((nfft - ncp +1):nfft)) Sifft_data];

        % arrange symbols row wise and serialize to transmit time domain signal
        Sofdm_tx = reshape(Scyclic_add_data.',nbitsym*nsym,1);
        
        %%
        % Channel Effect 1 
        St=(0:(nbitsym*nsym-1))*dt;
        Stheta=(theta_0(k)./(d0+v*St+0.5*a*St.^2).^3.346).';
        Sch_data=Sofdm_tx.*Stheta;                        
        
        %% 
        % Receiver
        %parallalize received data
        Spar_rec_data = reshape(Sch_data.', nbitsym, nsym).';

        % remove cyclic prefix
        Scyclic_pre_rem=Spar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        Sfft_data = fftshift(fft(Scyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        Srem_pilot = Sfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        Sch_rem_pilot = reshape((Srem_pilot.*alpha).',1,[])+noise1;
                
        % parallalize ch data
        Spardata=reshape(Sch_rem_pilot, subcar/2,[]).';

        % get received pilots by averaging SM OFDM symbols
        for s=1:nsym/SM
            Srecpilot(((s-1)*SM+1):(s*SM),:) = (real(repmat(mean(Spardata(((s-1)*SM+1):(s*SM),:)),SM,1)) + imag(repmat(mean(Spardata(((s-1)*SM+1):(s*SM),:)),SM,1)))/2;
        end
        
        % subtract pilots from the received signal
        Srecdata = Spardata - (1+i)*Srecpilot;
                
        % calculate SNR
        SEsN0(n,k)=mean(abs(Srecdata(symbols)).^2)/mean(abs(noise1).^2);

        % estimate qam_data_tx
        Sqam_data_rx=Srecdata./(Srecpilot*(sqrt(2)/sqrt(Salpha/(1-Salpha))));

        % serialize qam data
        Sserdata=reshape(Sqam_data_rx.', subcar*nsym/2,[]).';

        % demodulate the qam data in the ofdm symbols
        Sbin_data_rx = reshape(qamdemod(Sserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        SBER(n,k)=sum((data(bits)~=Sbin_data_rx(bits)))/length(bits);
        
       
        %% THIS IS FOR BLOCK TYPE CE ALGORITHM (1 data packet, 1 pilot packet)
        rng(n);
        % TX
        % parallalize bin data and do qam modulation
        Bqam_data_tx = reshape([(rand(subcar/2,nsym)>0.5)*2-1; reshape(qam_data_tx,subcar/2,[])],1,[]); 
        
        % parallalize qam data for collective ifft
        Bpar_data = reshape(Bqam_data_tx.', subcar/2, nsym*2).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        Bpilot_ins_data=[zeros(nsym*2,guard) Bpar_data zeros(nsym*2,1) conj(flip(Bpar_data.').') zeros(nsym*2,guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        Bifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(Bpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        Bcyclic_add_data = [Bifft_data(:,((nfft - ncp +1):nfft)) Bifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        Bofdm_tx = reshape(Bcyclic_add_data.',nbitsym*nsym*2,1);

        %%
        % Channel Effect 1 
        Bt=(0:(nbitsym*nsym*2-1))*dt;
        Btheta=(theta_0(k)./(d0+v*Bt+0.5*a*Bt.^2).^3.346).';
        Bch_data=Bofdm_tx.*Btheta;
        
        
        
        %% 
        % Receiver
        %parallalize received data
        Bpar_rec_data = reshape(Bch_data.', nbitsym, nsym*2).';

        % remove cyclic prefix
        Bcyclic_pre_rem=Bpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        Bfft_data = fftshift(fft(Bcyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        Brem_pilot = Bfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        Bch_rem_pilot = reshape((Brem_pilot.*alpha).',1,[])+noise2;
                
        % parallalize ch data
        Bpardata=reshape(Bch_rem_pilot, subcar,[]).';

        % separate data and training matrix
        Brecpilot=abs(Bpardata(:,1:subcar/2));
        Brecdata=Bpardata(:,subcar/2+1:end);
        
        % calculate SNR
        BEsN0(n,k)=mean(abs(Bch_rem_pilot(symbols)).^2)/mean(abs(noise2).^2);

        % estimate qam_data_tx
        Bqam_data_rx=Brecdata./Brecpilot;

        % serialize qam data
        Bserdata=reshape(Bqam_data_rx.', subcar*nsym/2,[]).';

        % demodulate the qam data in the ofdm symbols
        Bbin_data_rx = reshape(qamdemod(Bserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        BBER(n,k)=sum((data(Bbits)~=Bbin_data_rx(Bbits)))/length(Bbits);
        

        %% THIS IS FOR BLOCK TYPE CE ALGORITHM (2 data packet, 1 pilot packet)
        rng(n);
        % TX
        % parallalize bin data and do qam modulation
        BBqam_data_tx = reshape([(rand(subcar/2,nsym/2)>0.5)*2-1; reshape(qam_data_tx,subcar,[])],1,[]); % Pilot power = Mean Absolute Power of qam_data_tx
        
        % parallalize qam data for collective ifft
        BBpar_data = reshape(BBqam_data_tx.', subcar/2, nsym*1.5).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        BBpilot_ins_data=[zeros(nsym*1.5,guard) BBpar_data zeros(nsym*1.5,1) conj(flip(BBpar_data.').') zeros(nsym*1.5,guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        BBifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(BBpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        BBcyclic_add_data = [BBifft_data(:,((nfft - ncp +1):nfft)) BBifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        BBofdm_tx = reshape(BBcyclic_add_data.',nbitsym*nsym*1.5,1);

        %%
        % Channel Effect 1     
        BBt=(0:(nbitsym*nsym*1.5-1))*dt;
        BBtheta=(theta_0(k)./(d0+v*BBt+0.5*a*BBt.^2).^3.346).';
        BBch_data=BBofdm_tx.*BBtheta;
                
        
        %% 
        % Receiver        
        %parallalize received data
        BBpar_rec_data = reshape(BBch_data.', nbitsym, nsym*1.5).';

        % remove cyclic prefix
        BBcyclic_pre_rem=BBpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        BBfft_data = fftshift(fft(BBcyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        BBrem_pilot = BBfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        BBch_rem_pilot = reshape((BBrem_pilot.*alpha).',1,[])+noise3;
      
        % parallalize ch data
        BBpardata=reshape(BBch_rem_pilot, subcar*1.5,[]).';

        % separate data and training matrix
        BBrecpilot=abs(repmat(BBpardata(:,1:subcar/2),1,2));
        BBrecdata=BBpardata(:,subcar/2+1:end);
        
        % calculate SNR
        BBEsN0(n,k)=mean(abs(BBch_rem_pilot(symbols)).^2)/mean(abs(noise3).^2);

        % estimate qam_data_tx
        BBqam_data_rx=BBrecdata./BBrecpilot;

        % serialize qam data
        BBserdata=reshape(BBqam_data_rx.', subcar*nsym/2,[]).';

        % demodulate the qam data in the ofdm symbols
        BBbin_data_rx = reshape(qamdemod(BBserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        BBBER(n,k)=sum((data(BBbits)~=BBbin_data_rx(BBbits)))/length(BBbits);
        
        %% THIS IS FOR BLOCK TYPE CE ALGORITHM (3 data packet, 1 pilot packet)
        rng(n);
        % TX
        % parallalize bin data and do qam modulation
        BBBqam_data_tx = reshape([(rand(subcar/2,nsym/3)>0.5)*2-1; reshape(qam_data_tx,subcar*1.5,[])],1,[]);
        
        % parallalize qam data for collective ifft
        BBBpar_data = reshape(BBBqam_data_tx.', subcar/2, nsym*(1+1/3)).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        BBBpilot_ins_data=[zeros(nsym*(1+1/3),guard) BBBpar_data zeros(nsym*(1+1/3),1) conj(flip(BBBpar_data.').') zeros(nsym*(1+1/3),guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        BBBifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(BBBpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        BBBcyclic_add_data = [BBBifft_data(:,((nfft - ncp +1):nfft)) BBBifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        BBBofdm_tx = reshape(BBBcyclic_add_data.',nbitsym*nsym*(1+1/3),1);

        %%
        % Channel Effect 1     
        BBBt=(0:(nbitsym*nsym*(1+1/3)-1))*dt;
        BBBtheta=(theta_0(k)./(d0+v*BBBt+0.5*a*BBBt.^2).^3.346).';
        BBBch_data=BBBofdm_tx.*BBBtheta;
                
        
        %% 
        % Receiver        
        %parallalize received data
        BBBpar_rec_data = reshape(BBBch_data.', nbitsym, nsym*(1+1/3)).';

        % remove cyclic prefix
        BBBcyclic_pre_rem=BBBpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        BBBfft_data = fftshift(fft(BBBcyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        BBBrem_pilot = BBBfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        BBBch_rem_pilot = reshape((BBBrem_pilot.*alpha).',1,[])+noise4;
      
        % parallalize ch data
        BBBpardata=reshape(BBBch_rem_pilot, subcar*2,[]).';

        % separate data and training matrix
        BBBrecpilot=abs(repmat(BBBpardata(:,1:subcar/2),1,3));
        BBBrecdata=BBBpardata(:,subcar/2+1:end);
        
        % calculate SNR
        BBBEsN0(n,k)=mean(abs(BBBch_rem_pilot(symbols)).^2)/mean(abs(noise4).^2);

        % estimate qam_data_tx
        BBBqam_data_rx=BBBrecdata./BBBrecpilot;

        % serialize qam data
        BBBserdata=reshape(BBBqam_data_rx.', subcar*nsym/2,[]).';

        % demodulate the qam data in the ofdm symbols
        BBBbin_data_rx = reshape(qamdemod(BBBserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        BBBBER(n,k)=sum((data(BBBbits)~=BBBbin_data_rx(BBBbits)))/length(BBBbits);
        
        %% THIS IS FOR BLOCK TYPE CE ALGORITHM ("kb" data packet, 1 pilot packet)
        rng(n);
        % TX
        % parallalize bin data and do qam modulation
        BBBBqam_data_tx = reshape([(rand(subcar/2,nsym/kb)>0.5)*2-1; reshape(qam_data_tx,subcar*kb/2,[])],1,[]);
        
        % parallalize qam data for collective ifft
        BBBBpar_data = reshape(BBBBqam_data_tx.', subcar/2, round(nsym*(1+1/kb))).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        BBBBpilot_ins_data=[zeros(round(nsym*(1+1/kb)),guard) BBBBpar_data zeros(round(nsym*(1+1/kb)),1) conj(flip(BBBBpar_data.').') zeros(round(nsym*(1+1/kb)),guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        BBBBifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(BBBBpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        BBBBcyclic_add_data = [BBBBifft_data(:,((nfft - ncp +1):nfft)) BBBBifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        BBBBofdm_tx = reshape(BBBBcyclic_add_data.',round(nbitsym*nsym*(1+1/kb)),1);

        %%
        % Channel Effect 1     
        BBBBt=(0:(round(nbitsym*nsym*(1+1/kb))-1))*dt;
        BBBBtheta=(theta_0(k)./(d0+v*BBBBt+0.5*a*BBBBt.^2).^3.346).';
        BBBBch_data=BBBBofdm_tx.*BBBBtheta;
                
        
        %% 
        % Receiver        
        %parallalize received data
        BBBBpar_rec_data = reshape(BBBBch_data.', nbitsym, round(nsym*(1+1/kb))).';

        % remove cyclic prefix
        BBBBcyclic_pre_rem=BBBBpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        BBBBfft_data = fftshift(fft(BBBBcyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        BBBBrem_pilot = BBBBfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        BBBBch_rem_pilot = reshape((BBBBrem_pilot.*alpha).',1,[])+noisekb;
      
        % parallalize ch data
        BBBBpardata=reshape(BBBBch_rem_pilot, subcar*(kb+1)/2,[]).';

        % separate data and training matrix
        BBBBrecpilot=abs(repmat(BBBBpardata(:,1:subcar/2),1,kb));
        BBBBrecdata=BBBBpardata(:,subcar/2+1:end);
        
        % calculate SNR
        BBBBEsN0(n,k)=mean(abs(BBBBch_rem_pilot(symbols)).^2)/mean(abs(noisekb).^2);

        % estimate qam_data_tx
        BBBBqam_data_rx=BBBBrecdata./BBBBrecpilot;

        % serialize qam data
        BBBBserdata=reshape(BBBBqam_data_rx.', subcar*nsym/2,[]).';

        % demodulate the qam data in the ofdm symbols
        BBBBbin_data_rx = reshape(qamdemod(BBBBserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        BBBBBER(n,k)=sum((data(BBBBbits)~=BBBBbin_data_rx(BBBBbits)))/length(BBBBbits);
        
        %% THIS IS FOR COMB TYPE CE ALGORITHM (1 data symbol, 1 pilot symbol)
        rng(n);
        % TX
        % parallalize bin data and do qam modulation
        Cqam_data_tx = reshape([(rand(1,nsym*subcar/2)>0.5)*2-1; qam_data_tx],1,[]); 
        
        % parallalize qam data for collective ifft
        Cpar_data = reshape(Cqam_data_tx.', subcar/2, nsym*2).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        Cpilot_ins_data=[zeros(nsym*2,guard) Cpar_data zeros(nsym*2,1) conj(flip(Cpar_data.').') zeros(nsym*2,guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        Cifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(Cpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        Ccyclic_add_data = [Cifft_data(:,((nfft - ncp +1):nfft)) Cifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        Cofdm_tx = reshape(Ccyclic_add_data.',nbitsym*nsym*2,1);

        %%
        % Channel Effect 1  
        Ct=(0:(nbitsym*nsym*2-1))*dt;
        Ctheta=(theta_0(k)./(d0+v*Ct+0.5*a*Ct.^2).^3.346).';
        Cch_data=Cofdm_tx.*Ctheta;
        
        
        
        %% 
        % Receiver
        %parallalize received data
        Cpar_rec_data = reshape(Cch_data.', nbitsym, nsym*2).';

        % remove cyclic prefix
        Ccyclic_pre_rem=Cpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        Cfft_data = fftshift(fft(Ccyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        Crem_pilot = Cfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        Cch_rem_pilot = reshape((Crem_pilot.*alpha).',1,[])+noise2;
        
        % Parallelize data
        Cpardata=reshape(Cch_rem_pilot,subcar/2,[]).';
        
        % separate data and training matrix
        Crecpilot=abs(Cpardata(:,1:2:end));
        Crecdata=Cpardata(:,2:2:end);
        
        % calculate SNR
        CEsN0(n,k)=mean(abs(Cch_rem_pilot(symbols)).^2)/mean(abs(noise2).^2);

        % estimate qam_data_tx
        tmp= interp1(1:(subcar/2)/2,Crecpilot.',1:0.5:((subcar/2)+0.5),'linear','extrap');
        Cresp=tmp(2:2:subcar/2,:).';
        Cqam_data_rx=Crecdata./Cresp;

        % serialize qam data
        Cserdata=reshape(Cqam_data_rx.', 1,[]).';

        % demodulate the qam data in the ofdm symbols
        Cbin_data_rx = reshape(qamdemod(Cserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        CBER(n,k)=sum((data(Cbits)~=Cbin_data_rx(Cbits)))/length(Cbits);

%% THIS IS FOR COMB TYPE CE ALGORITHM (2 data symbol, 1 pilot symbol)
rng(n);
        % TX
        % parallalize bin data and do qam modulation
        CCqam_data_tx = reshape([(rand(1,nsym*subcar/4)>0.5)*2-1; reshape(qam_data_tx, 2,[])],1,[]); 
        
        % parallalize qam data for collective ifft
        CCpar_data = reshape(CCqam_data_tx.', subcar/2, nsym*(1+1/2)).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        CCpilot_ins_data=[zeros(nsym*(1+1/2),guard) CCpar_data zeros(nsym*(1+1/2),1) conj(flip(CCpar_data.').') zeros(nsym*(1+1/2),guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        CCifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(CCpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        CCcyclic_add_data = [CCifft_data(:,((nfft - ncp +1):nfft)) CCifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        CCofdm_tx = reshape(CCcyclic_add_data.',nbitsym*nsym*(1+1/2),1);

        %%
        % Channel Effect 1  
        CCt=(0:(nbitsym*nsym*(1+1/2)-1))*dt;
        CCtheta=(theta_0(k)./(d0+v*CCt+0.5*a*CCt.^2).^3.346).';
        CCch_data=CCofdm_tx.*CCtheta;
        
        
        
        %% 
        % Receiver
        %parallalize received data
        CCpar_rec_data = reshape(CCch_data.', nbitsym, nsym*(1+1/2)).';

        % remove cyclic prefix
        CCcyclic_pre_rem=CCpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        CCfft_data = fftshift(fft(CCcyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        CCrem_pilot = CCfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        CCch_rem_pilot = reshape((CCrem_pilot.*alpha).',1,[])+noise3;
        
        % Parallelize data
        CCpardata=reshape(CCch_rem_pilot,subcar/2,[]).';
        
        % separate data and training matrix
        CCrecpilot=abs(CCpardata(:,1:3:end));
        CCrecdata=CCpardata(:,sort([2:3:subcar/2 3:3:subcar/2]));
        
        % calculate SNR
        CCEsN0(n,k)=mean(abs(CCch_rem_pilot(symbols)).^2)/mean(abs(noise3).^2);

        % estimate qam_data_tx
        tmp= interp1(1:3:subcar/2,CCrecpilot.',1:subcar/2,'linear','extrap');
        CCresp=tmp(sort([2:3:subcar/2 3:3:subcar/2]),:).';
        CCqam_data_rx=CCrecdata./CCresp;

        % serialize qam data
        CCserdata=reshape(CCqam_data_rx.', 1,[]).';

        % demodulate the qam data in the ofdm symbols
        CCbin_data_rx = reshape(qamdemod(CCserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        CCBER(n,k)=sum((data(CCbits)~=CCbin_data_rx(CCbits)))/length(CCbits);

%% THIS IS FOR COMB TYPE CE ALGORITHM (3 data symbol, 1 pilot symbol)
rng(n);
        % TX
        % parallalize bin data and do qam modulation
        CCCqam_data_tx = reshape([(rand(1,nsym*subcar/6)>0.5)*2-1; reshape(qam_data_tx, 3,[])],1,[]); 
        
        % parallalize qam data for collective ifft
        CCCpar_data = reshape(CCCqam_data_tx.', subcar/2, nsym*(1+1/3)).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        CCCpilot_ins_data=[zeros(nsym*(1+1/3),guard) CCCpar_data zeros(nsym*(1+1/3),1) conj(flip(CCCpar_data.').') zeros(nsym*(1+1/3),guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        CCCifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(CCCpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        CCCcyclic_add_data = [CCCifft_data(:,((nfft - ncp +1):nfft)) CCCifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        CCCofdm_tx = reshape(CCCcyclic_add_data.',nbitsym*nsym*(1+1/3),1);

        %%
        % Channel Effect 1  
        CCCt=(0:(nbitsym*nsym*(1+1/3)-1))*dt;
        CCCtheta=(theta_0(k)./(d0+v*CCCt+0.5*a*CCCt.^2).^3.346).';
        CCCch_data=CCCofdm_tx.*CCCtheta;
        
        
        
        %% 
        % Receiver
        %parallalize received data
        CCCpar_rec_data = reshape(CCCch_data.', nbitsym, nsym*(1+1/3)).';

        % remove cyclic prefix
        CCCcyclic_pre_rem=CCCpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        CCCfft_data = fftshift(fft(CCCcyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        CCCrem_pilot = CCCfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        CCCch_rem_pilot = reshape((CCCrem_pilot.*alpha).',1,[])+noise4;
        
        % Parallelize data
        CCCpardata=reshape(CCCch_rem_pilot,subcar/2,[]).';
        
        % separate data and training matrix
        CCCrecpilot=abs(CCCpardata(:,1:4:end));
        CCCrecdata=CCCpardata(:,sort([2:4:subcar/2 3:4:subcar/2 4:4:subcar/2]));
        
        % calculate SNR
        CCCEsN0(n,k)=mean(abs(CCCch_rem_pilot(symbols)).^2)/mean(abs(noise4).^2);

        % estimate qam_data_tx
        tmp= interp1(1:4:subcar/2,CCCrecpilot.',1:subcar/2,'linear','extrap');
        CCCresp=tmp(sort([2:4:subcar/2 3:4:subcar/2 4:4:subcar/2]),:).';
        CCCqam_data_rx=CCCrecdata./CCCresp;

        % serialize qam data
        CCCserdata=reshape(CCCqam_data_rx.', 1,[]).';

        % demodulate the qam data in the ofdm symbols
        CCCbin_data_rx = reshape(qamdemod(CCCserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        CCCBER(n,k)=sum((data(CCCbits)~=CCCbin_data_rx(CCCbits)))/length(CCCbits);
        
%% THIS IS FOR COMB TYPE CE ALGORITHM ("kc" data symbol, 1 pilot symbol)
rng(n);
        % TX
        % parallalize bin data and do qam modulation
        CCCCqam_data_tx = reshape([(rand(1,nsym*subcar/(2*kc))>0.5)*2-1; reshape(qam_data_tx, kc,[])],1,[]); 
        
        % parallalize qam data for collective ifft
        CCCCpar_data = reshape(CCCCqam_data_tx.', subcar/2, nsym*(1+1/kc)).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        CCCCpilot_ins_data=[zeros(nsym*(1+1/kc),guard) CCCCpar_data zeros(nsym*(1+1/kc),1) conj(flip(CCCCpar_data.').') zeros(nsym*(1+1/kc),guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        CCCCifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(CCCCpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        CCCCcyclic_add_data = [CCCCifft_data(:,((nfft - ncp +1):nfft)) CCCCifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        CCCCofdm_tx = reshape(CCCCcyclic_add_data.',nbitsym*nsym*(1+1/kc),1);

        %%
        % Channel Effect 1  
        CCCCt=(0:(nbitsym*nsym*(1+1/kc)-1))*dt;
        CCCCtheta=(theta_0(k)./(d0+v*CCCCt+0.5*a*CCCCt.^2).^3.346).';
        CCCCch_data=CCCCofdm_tx.*CCCCtheta;
        
        
        
        %% 
        % Receiver
        %parallalize received data
        CCCCpar_rec_data = reshape(CCCCch_data.', nbitsym, nsym*(1+1/kc)).';

        % remove cyclic prefix
        CCCCcyclic_pre_rem=CCCCpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        CCCCfft_data = fftshift(fft(CCCCcyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the qam data
        CCCCrem_pilot = CCCCfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        CCCCch_rem_pilot = reshape((CCCCrem_pilot.*alpha).',1,[])+noisekc;
        
        % Parallelize data
        CCCCpardata=reshape(CCCCch_rem_pilot,subcar/2,[]).';
        
        % separate data and training matrix
        CCCCrecpilot=abs(CCCCpardata(:,1:(kc+1):end));
        CCCCrecdata=CCCCpardata(:,sort(setdiff(1:subcar/2,1:(kc+1):subcar/2)));
        
        % calculate SNR
        CCCCEsN0(n,k)=mean(abs(CCCCch_rem_pilot(symbols)).^2)/mean(abs(noisekc).^2);

        % estimate qam_data_tx
        tmp= interp1(1:(kc+1):subcar/2,CCCCrecpilot.',1:subcar/2,'linear','extrap');
        CCCCresp=tmp(sort(setdiff(1:subcar/2,1:(kc+1):subcar/2)),:).';
        CCCCqam_data_rx=CCCCrecdata./CCCCresp;

        % serialize qam data
        CCCCserdata=reshape(CCCCqam_data_rx.', 1,[]).';

        % demodulate the qam data in the ofdm symbols
        CCCCbin_data_rx = reshape(qamdemod(CCCCserdata, qam,'OutputType','bit','UnitAveragePower',true),1,[]);

        % calculate BER
        CCCCBER(n,k)=sum((data(CCCCbits)~=CCCCbin_data_rx(CCCCbits)))/length(CCCCbits);
        
    end
   
end

%% calculate SNR and BER
EbN0=mean(EsN0)/nbitqam;
EbN0_dB=10*log10(EbN0);
BER=mean(BER);
CEbN0=mean(CEsN0)/nbitqam;
CEbN0_dB=10*log10(CEbN0);
CBER=mean(CBER);
CCEbN0=mean(CCEsN0)/nbitqam;
CCEbN0_dB=10*log10(CCEbN0);
CCBER=mean(CCBER);
CCCEbN0=mean(CCCEsN0)/nbitqam;
CCCEbN0_dB=10*log10(CCCEbN0);
CCCBER=mean(CCCBER);
CCCCEbN0=mean(CCCCEsN0)/nbitqam;
CCCCEbN0_dB=10*log10(CCCCEbN0);
CCCCBER=mean(CCCCBER);
BEbN0=mean(BEsN0)/nbitqam;
BEbN0_dB=10*log10(BEbN0);
BBER=mean(BBER);
BBEbN0=mean(BBEsN0)/nbitqam;
BBEbN0_dB=10*log10(BBEbN0);
BBBER=mean(BBBER);
BBBEbN0=mean(BBBEsN0)/nbitqam;
BBBEbN0_dB=10*log10(BBBEbN0);
BBBBER=mean(BBBBER);
BBBBEbN0=mean(BBBBEsN0)/nbitqam;
BBBBEbN0_dB=10*log10(BBBBEbN0);
BBBBBER=mean(BBBBBER);
SEbN0=mean(SEsN0)/nbitqam;
SEbN0_dB=10*log10(SEbN0);
SBER=mean(SBER);

% %% Plot SNR vs. BER for block type CE
% figure
% semilogy(BEbN0_dB,BBER,'p-')
% hold on
% semilogy(BBEbN0_dB,BBBER,'s-') 
% hold on
% semilogy(BBBEbN0_dB,BBBBER,'d-') 
% hold on
% semilogy(BBBBEbN0_dB,BBBBBER,'h-') 
% hold on
% semilogy(CEbN0_dB,CBER,'+-')
% hold on
% semilogy(CCEbN0_dB,CCBER,'x-')
% hold on
% semilogy(CCCEbN0_dB,CCCBER,'*-')
% hold on
% semilogy(CCCCEbN0_dB,CCCCBER,'v-')
% hold on
% semilogy(SEbN0_dB,SBER,'o-')
% hold on
% semilogy(EbN0_dB,BER,'^-')
% 
% legend('Block type CE (1 data - 1 pilot packet)','Block type CE (2 data - 1 pilot packet)','Block type CE (3 data - 1 pilot packet)',...
%     ['Block type CE (' num2str(kb) ' data - 1 pilot packet)'],'Comb type CE (1 data - 1 pilot symbol)','Comb type CE (2 data - 1 pilot symbol)',...
%     'Comb type CE (3 data - 1 pilot symbol)',['Comb type CE (' num2str(kc) ' data - 1 pilot symbol)'],['ST-based CE, M=' num2str(SM) ' \alpha=' num2str(Salpha)],'Proposed CE')
% xlabel('Eb/N0 (dB)')
% ylabel('BER')
% title('CE Performance Comparison')

const = (nbitqam*subcar/2)/(nfft+ncp);

figure
semilogy(const*1/2*(1-BBER),'pb-')
hold on
semilogy(const*2/3*(1-BBBER),'sr-') 
hold on
semilogy(const*3/4*(1-BBBBER),'dy-')
hold on
semilogy(const*1/2*(1-CBER),'+k-')
hold on
semilogy(const*2/3*(1-CCBER),'xc-')
hold on
semilogy(const*3/4*(1-CCCBER),'*m-')
hold on
semilogy(const*(1-SBER),'ob-')
hold on
semilogy(const*(1-BER),'^g-')
set(gca,'xticklabel',ambient_light.')
legend('Block PS1','Block PS2','Block PS3','Comb PS1','Comb PS2','Comb PS3','ST','Proposed')
xlabel('Ambient Light')
ylabel('Spectral Efficiency (C [bits/s/Hz])')
title('8m 0deg: a:parking lot b: night lights off c: night lights on d: sunrise partial cloudy e: sunny f: cloudy g: partial cloudy h: sunset i: cloudy sunset j: shadow')
grid on


