%% DCO-OFDM, Only Half Data Subcarriers are used for psk Data for Hermitian Symmetry
clear
clc
close all

%% Local Variables: number of bits, samples, guard symbols, temp variables for interfunction calls
psk=16;
nfft=64; %fft size
subcar=48; %number of data subcarriers
nbitpsk = log2(psk);
guard = (nfft-subcar)/2;
nsym=3000; %number of ofdm symbols
nbitdata = subcar*nbitpsk*nsym/2;

ncp=4; %cyclic prefix size
nbitsym = nfft + ncp;

data=(rand(1,nbitdata)>0.5)+0;

% magnitude of the normalized CFR
alpha=[0.609981229934408,0.715459955559726,0.636965602683617,0.694127681700978,0.695098676789100,0.681699620219806,0.705628746895461,0.834320339106743,0.873378196662836,0.827668204581301,0.768804769140058,0.756148741911805,0.926978115212095,0.913129736817394,0.871793084001999,0.848536621661314,0.909045348037926,0.834381831138879,0.896913346571862,0.928182910840626,0.932323910122903,0.907741812468878,0.943650044868540,1];

symbols=1:nsym*subcar/2;
bits=1:nsym*nbitpsk*subcar/2;
Cbits=1:nsym*nbitpsk*subcar/2*(1/2);
CCbits=1:nsym*nbitpsk*subcar/2*(2/3);
CCCbits=1:nsym*nbitpsk*subcar/2*(3/4);
Bbits=1:nsym*nbitpsk*subcar/2*(1/2);
BBbits=1:nsym*nbitpsk*subcar/2*(2/3);
BBBbits=1:nsym*nbitpsk*subcar/2*(3/4);

%theta = logspace(-1,4,100); 

% theta for dynamic case
d0=1;
%theta_0 corresponds to the EIRP, and is swept to get different SNRs under noise
theta_0=logspace(-1,4,100);
v=0; %m/s2
a=0; %m/s2
dt=1/(2*1e6); %CFR is measured for 2MHz bandwidth

for n=1:10
    %Noise vectors are created with a normalized power. One can think that
    %the normaliation factor is rescaled with noise power.
    noise1=normrnd(0,1,[1,subcar*nsym/2])+i*normrnd(0,1,[1,subcar*nsym/2]);
    noise2=normrnd(0,1,[1,subcar*nsym])+i*normrnd(0,1,[1,subcar*nsym]);
    noise3=normrnd(0,1,[1,subcar*nsym*0.75])+i*normrnd(0,1,[1,subcar*nsym*0.75]);
    noise4=normrnd(0,1,[1,subcar*nsym/2*(1+1/3)])+i*normrnd(0,1,[1,subcar*nsym/2*(1+1/3)]);
    
    for k=1:length(theta_0)
        %% THIS IS FOR THE PROPOSED ALGORITHM
        % TX
        % parallalize bin data and do psk modulation
        psk_data_tx = pskmod(bi2de(reshape(data, nbitpsk,[]).').',psk);
        
        % parallalize psk data for collective ifft
        par_data = reshape(psk_data_tx.', subcar/2, nsym).';

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
        
        % remove the guard to get the psk data
        rem_pilot = fft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        ch_rem_pilot = reshape((rem_pilot.*alpha).',1,[])+noise1;
        
        % calculate SNR
        EsN0(n,k)=mean(abs(ch_rem_pilot(symbols)).^2)/mean(abs(noise1).^2);
        
        % parallalize ch data
        pardata=reshape(ch_rem_pilot, subcar/2,[]).';

        % estimate theta
        pskconst=pskmod(0:(psk-1),psk);
        est_theta=sum(((abs(real(pardata))+abs(imag(pardata)))./alpha).')./(mean(abs(real(pskconst))+abs(imag(pskconst)))*subcar/2);

        % calculate theta*data
        td=(pardata./alpha).';

        % estimate psk_data_tx
        psk_data_rx=td./est_theta;

        % serialize psk data
        serdata=reshape(psk_data_rx, subcar*nsym/2,[]).';

        % demodulate the psk data in the ofdm symbols
        bin_data_rx = reshape(de2bi(reshape(pskdemod(serdata, psk),1,[])).',1,[]);

        % calculate BER
        BER(n,k)=sum((data(bits)~=bin_data_rx(bits)))/length(bits);
        
        
        %% THIS IS FOR BLOCK TYPE CE ALGORITHM (1 data packet, 1 pilot packet)
        % TX
        % parallalize bin data and do psk modulation
        Bpsk_data_tx = reshape([(rand(subcar/2,nsym)>0.5)*2-1; reshape(psk_data_tx,subcar/2,[])],1,[]); 
        
        % parallalize psk data for collective ifft
        Bpar_data = reshape(Bpsk_data_tx.', subcar/2, nsym*2).';

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
        
        % remove the guard to get the psk data
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

        % estimate psk_data_tx
        Bpsk_data_rx=Brecdata./Brecpilot;

        % serialize psk data
        Bserdata=reshape(Bpsk_data_rx.', subcar*nsym/2,[]).';

        % demodulate the psk data in the ofdm symbols
        Bbin_data_rx = reshape(de2bi(reshape(pskdemod(Bserdata, psk),1,[])).',1,[]);

        % calculate BER
        BBER(n,k)=sum((data(Bbits)~=Bbin_data_rx(Bbits)))/length(Bbits);
       

        %% THIS IS FOR COMB TYPE CE ALGORITHM (1 data symbol, 1 pilot symbol)
        % TX
        % parallalize bin data and do psk modulation
        Cpsk_data_tx = reshape([(rand(1,nsym*subcar/2)>0.5)*2-1; psk_data_tx],1,[]); 
        
        % parallalize psk data for collective ifft
        Cpar_data = reshape(Cpsk_data_tx.', subcar/2, nsym*2).';

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
        
        % remove the guard to get the psk data
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

        % estimate psk_data_tx
        tmp= interp1(1:(subcar/2)/2,Crecpilot.',1:0.5:((subcar/2)+0.5),'linear','extrap');
        Cresp=tmp(2:2:subcar/2,:).';
        Cqam_data_rx=Crecdata./Cresp;

        % serialize psk data
        Cserdata=reshape(Cqam_data_rx.', 1,[]).';

        % demodulate the psk data in the ofdm symbols
        Cbin_data_rx = reshape(de2bi(reshape(pskdemod(Cserdata, psk),1,[])).',1,[]);

        % calculate BER
        CBER(n,k)=sum((data(Cbits)~=Cbin_data_rx(Cbits)))/length(Cbits);

        %% THIS IS FOR BLOCK TYPE CE ALGORITHM (2 data packet, 1 pilot packet)
        % TX
        % parallalize bin data and do psk modulation
        BBpsk_data_tx = reshape([(rand(subcar/2,nsym/2)>0.5)*2-1; reshape(psk_data_tx,subcar,[])],1,[]); % Pilot power = Mean Absolute Power of psk_data_tx
        
        % parallalize psk data for collective ifft
        BBpar_data = reshape(BBpsk_data_tx.', subcar/2, nsym*1.5).';

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
        
        % remove the guard to get the psk data
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

        % estimate psk_data_tx
        BBpsk_data_rx=BBrecdata./BBrecpilot;

        % serialize psk data
        BBserdata=reshape(BBpsk_data_rx.', subcar*nsym/2,[]).';

        % demodulate the psk data in the ofdm symbols
        BBbin_data_rx = reshape(de2bi(reshape(pskdemod(BBserdata, psk),1,[])).',1,[]);

        % calculate BER
        BBBER(n,k)=sum((data(BBbits)~=BBbin_data_rx(BBbits)))/length(BBbits);


        %% THIS IS FOR COMB TYPE CE ALGORITHM (2 data symbol, 1 pilot symbol)
        % TX
        % parallalize bin data and do psk modulation
        CCpsk_data_tx = reshape([(rand(1,nsym*subcar/4)>0.5)*2-1; reshape(psk_data_tx, 2,[])],1,[]); 
        
        % parallalize psk data for collective ifft
        CCpar_data = reshape(CCpsk_data_tx.', subcar/2, nsym*1.5).';

        % generate hermitian symmetric data and insert guard at non-data subcarrier locations
        CCpilot_ins_data=[zeros(nsym*1.5,guard) CCpar_data zeros(nsym*1.5,1) conj(flip(CCpar_data.').') zeros(nsym*1.5,guard-1)] ;

        % take ifft of all the ofdm symbols in the frame
        CCifft_data = sqrt(nfft*nfft/subcar)*ifft(ifftshift(CCpilot_ins_data.',1)).';

        % add cyclic prefix to each ofdm symbol
        CCcyclic_add_data = [CCifft_data(:,((nfft - ncp +1):nfft)) CCifft_data];

        % arrange symbols row wise and searialize to transmit time domain signal
        CCofdm_tx = reshape(CCcyclic_add_data.',nbitsym*nsym*1.5,1);

        %%
        % Channel Effect 1  
        CCt=(0:(nbitsym*nsym*(1+1/2)-1))*dt;
        CCtheta=(theta_0(k)./(d0+v*CCt+0.5*a*CCt.^2).^3.346).';
        CCch_data=CCofdm_tx.*CCtheta;
                        
        
        %% 
        % Receiver        
        %parallalize received data
        CCpar_rec_data = reshape(CCch_data.', nbitsym, nsym*1.5).';

        % remove cyclic prefix
        CCcyclic_pre_rem=CCpar_rec_data(:,ncp+1:end);

        % take fft of time domain ofdm symbols
        CCfft_data = fftshift(fft(CCcyclic_pre_rem.'),1).'/sqrt(nfft*nfft/subcar);
        
        % remove the guard to get the psk data
        CCrem_pilot = CCfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        CCch_rem_pilot = reshape((CCrem_pilot.*alpha).',1,[])+noise3;
     
        % parallalize ch data
        CCpardata=reshape(CCch_rem_pilot,subcar/2,[]).';

       % separate data and training matrix
        CCrecpilot=abs(CCpardata(:,1:3:end));
        CCrecdata=CCpardata(:,sort([2:3:subcar/2 3:3:subcar/2]));
        
        % calculate SNR
        CCEsN0(n,k)=mean(abs(CCch_rem_pilot(symbols)).^2)/mean(abs(noise3).^2);

        % estimate psk_data_tx
        tmp= interp1(1:3:subcar/2,CCrecpilot.',1:subcar/2,'linear','extrap');
        CCresp=tmp(sort([2:3:subcar/2 3:3:subcar/2]),:).';
        CCqam_data_rx=CCrecdata./CCresp;

        % serialize psk data
        CCserdata=reshape(CCqam_data_rx.', 1,[]).';

        % demodulate the psk data in the ofdm symbols
        CCbin_data_rx = reshape(de2bi(reshape(pskdemod(CCserdata, psk),1,[])).',1,[]);

        % calculate BER
        CCBER(n,k)=sum((data(CCbits)~=CCbin_data_rx(CCbits)))/length(CCbits);

        %% THIS IS FOR BLOCK TYPE CE ALGORITHM (3 data packet, 1 pilot packet)
        % TX
        % parallalize bin data and do psk modulation
        BBBpsk_data_tx = reshape([(rand(subcar/2,nsym/3)>0.5)*2-1; reshape(psk_data_tx,subcar*1.5,[])],1,[]);
        
        % parallalize psk data for collective ifft
        BBBpar_data = reshape(BBBpsk_data_tx.', subcar/2, nsym*(1+1/3)).';

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
        
        % remove the guard to get the psk data
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

        % estimate psk_data_tx
        BBBpsk_data_rx=BBBrecdata./BBBrecpilot;

        % serialize psk data
        BBBserdata=reshape(BBBpsk_data_rx.', subcar*nsym/2,[]).';

        % demodulate the psk data in the ofdm symbols
        BBBbin_data_rx = reshape(de2bi(reshape(pskdemod(BBBserdata, psk),1,[])).',1,[]);

        % calculate BER
        BBBBER(n,k)=sum((data(BBBbits)~=BBBbin_data_rx(BBBbits)))/length(BBBbits);

%% THIS IS FOR COMB TYPE CE ALGORITHM (3 data symbol, 1 pilot symbol)
        % TX
        % parallalize bin data and do psk modulation
        CCCpsk_data_tx = reshape([(rand(1,nsym*subcar/6)>0.5)*2-1; reshape(psk_data_tx, 3,[])],1,[]); 
        
        % parallalize psk data for collective ifft
        CCCpar_data = reshape(CCCpsk_data_tx.', subcar/2, nsym*(1+1/3)).';

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
        
        % remove the guard to get the psk data
        CCCrem_pilot = CCCfft_data(:,guard+(1:subcar/2));
    
        % Channel Effect 2
        CCCch_rem_pilot = reshape((CCCrem_pilot.*alpha).',1,[])+noise4;
     
        % parallalize ch data
        CCCpardata=reshape(CCCch_rem_pilot,subcar/2,[]).';

        % separate data and training matrix
        CCCrecpilot=abs(CCCpardata(:,1:4:end));
        CCCrecdata=CCCpardata(:,sort([2:4:subcar/2 3:4:subcar/2 4:4:subcar/2]));
        
        % calculate SNR
        CCCEsN0(n,k)=mean(abs(CCCch_rem_pilot(symbols)).^2)/mean(abs(noise4).^2);

        % estimate psk_data_tx
        tmp= interp1(1:4:subcar/2,CCCrecpilot.',1:subcar/2,'linear','extrap');
        CCCresp=tmp(sort([2:4:subcar/2 3:4:subcar/2 4:4:subcar/2]),:).';
        CCCqam_data_rx=CCCrecdata./CCCresp;

        % serialize psk data
        CCCserdata=reshape(CCCqam_data_rx.', 1,[]).';

        % demodulate the psk data in the ofdm symbols
        CCCbin_data_rx = reshape(de2bi(reshape(pskdemod(CCCserdata, psk),1,[])).',1,[]);

        % calculate BER
        CCCBER(n,k)=sum((data(CCCbits)~=CCCbin_data_rx(CCCbits)))/length(CCCbits);

    end
   
end

%% calculate SNR and BER
EbN0=mean(EsN0)/nbitpsk;
EbN0_dB=10*log10(EbN0);
BEbN0=mean(BEsN0)/nbitpsk;
BEbN0_dB=10*log10(BEbN0);
BBEbN0=mean(BBEsN0)/nbitpsk;
BBEbN0_dB=10*log10(BBEbN0);
BBBEbN0=mean(BBBEsN0)/nbitpsk;
BBBEbN0_dB=10*log10(BBBEbN0);
CEbN0=mean(CEsN0)/nbitpsk;
CEbN0_dB=10*log10(CEbN0);
CCEbN0=mean(CCEsN0)/nbitpsk;
CCEbN0_dB=10*log10(CCEbN0);
CCCEbN0=mean(CCCEsN0)/nbitpsk;
CCCEbN0_dB=10*log10(CCCEbN0);

BER=mean(BER);
BBER=mean(BBER);
BBBER=mean(BBBER);
BBBBER=mean(BBBBER);
CBER=mean(CBER);
CCBER=mean(CCBER);
CCCBER=mean(CCCBER);

%% Plot SNR vs. BER for block type CE
figure
semilogy(BEbN0_dB,BBER,'p-')
hold on
semilogy(BBEbN0_dB,BBBER,'s-') 
hold on
semilogy(BBBEbN0_dB,BBBBER,'d-') 
hold on
semilogy(CEbN0_dB,CBER,'+-')
hold on
semilogy(CCEbN0_dB,CCBER,'x-')
hold on
semilogy(CCCEbN0_dB,CCCBER,'*-')
hold on
semilogy(EbN0_dB,BER,'^-')

legend('Block type CE (1 data - 1 pilot packet)','Block type CE (2 data - 1 pilot packet)','Block type CE (3 data - 1 pilot packet)',...
    'Comb type CE (1 data - 1 pilot symbol)','Comb type CE (2 data - 1 pilot symbol)','Comb type CE (3 data - 1 pilot symbol)','Proposed CE')
xlabel('Eb/N0 (dB)')
ylabel('BER')
title('CE Performance Comparison')

figure
semilogy(BEbN0_dB,nbitpsk*1/2*(1-BBER),'pb-')
hold on
semilogy(BBEbN0_dB,nbitpsk*2/3*(1-BBBER),'sr-') 
hold on
semilogy(BBBEbN0_dB,nbitpsk*3/4*(1-BBBBER),'dy-')
hold on
semilogy(CEbN0_dB,nbitpsk*1/2*(1-CBER),'+k-')
hold on
semilogy(CCEbN0_dB,nbitpsk*2/3*(1-CCBER),'xc-')
hold on
semilogy(CCCEbN0_dB,nbitpsk*3/4*(1-CCCBER),'*m-')
hold on
semilogy(EbN0_dB,nbitpsk*(1-BER),'^g-')
legend('Block PS1','Block PS2','Block PS3','Comb PS1','Comb PS2','Comb PS3','Proposed')
xlabel('Eb/N0 (dB)')
ylabel('Average Throughput (per subcarrier)')