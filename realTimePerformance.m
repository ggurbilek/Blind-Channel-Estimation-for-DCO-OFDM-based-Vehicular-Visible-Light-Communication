clear all
close all
clc

load('sumodata.mat')

figure
plot(vtime,dist)
xlabel('Time (s)')
ylabel('Distance (m)')

t=0:0.0000005:vtime(end);
d=interp1(vtime,dist,t);

const=[4 16 64];
for p=1:3
    %% Local Variables: number of bits, samples, guard symbols, temp variables for interfunction calls
    qam=const(p);
    nfft=64; %fft size
    subcar=48; %number of data subcarriers
    nbitqam = log2(qam);
    guard = (nfft-subcar)/2;
    nsym=10000; %number of ofdm symbols
    nbitdata = subcar*nbitqam*nsym/2;

    ncp=4; %cyclic prefix size
    nbitsym = nfft + ncp;

    data=(rand(1,nbitdata)>0.5)+0;

    % magnitude of the normalized CFR
    alpha=[0.609981229934408,0.715459955559726,0.636965602683617,0.694127681700978,0.695098676789100,0.681699620219806,0.705628746895461,0.834320339106743,0.873378196662836,0.827668204581301,0.768804769140058,0.756148741911805,0.926978115212095,0.913129736817394,0.871793084001999,0.848536621661314,0.909045348037926,0.834381831138879,0.896913346571862,0.928182910840626,0.932323910122903,0.907741812468878,0.943650044868540,1];

    symbols=1:nsym*subcar/2;
    bits=1:nsym*nbitqam*subcar/2;

    % theta
    theta_0=db2mag(8.114+30); %dBm to dBW

    % noise power (receiver sensitivity)
    np_dbm=-95;
    np_dbW=np_dbm+30;
    np=db2mag(np_dbW);
    
    for n=1:10
        
        noise1=np*(normrnd(0,1,[1,subcar*nsym/2])+i*normrnd(0,1,[1,subcar*nsym/2]));

        for k=1:floor(length(t)/(nbitsym*nsym))
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
            tt=t((nbitsym*nsym)*(k-1)+1:(k*nbitsym*nsym));
            tt0(k)=t((nbitsym*nsym)*(k-1)+1);
            dd=d((nbitsym*nsym)*(k-1)+1:(k*nbitsym*nsym));
            theta=(theta_0./dd.^3.346).';
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

            % calculate SNR and BER
            BER(n,k)=sum((data(bits)~=bin_data_rx(bits)))/length(bits);


            %% THIS IS FOR THE PROPOSED ALGORITHM (PSK)
            % TX
            % parallalize bin data and do qam modulation
            psk_data_tx = pskmod(bi2de(reshape(data, nbitqam,[]).').',qam);

            % parallalize qam data for collective ifft
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
            tt=t((nbitsym*nsym)*(k-1)+1:(k*nbitsym*nsym));
            tt0(k)=t((nbitsym*nsym)*(k-1)+1);
            dd=d((nbitsym*nsym)*(k-1)+1:(k*nbitsym*nsym));
            theta=(theta_0./dd.^3.346).';
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
            PEsN0(n,k)=mean(abs(ch_rem_pilot(symbols)).^2)/mean(abs(noise1).^2);

            % parallalize ch data
            pardata=reshape(ch_rem_pilot, subcar/2,[]).';

            % estimate theta
            pskconst=pskmod(0:(qam-1),qam);
            est_theta=sum(((abs(real(pardata))+abs(imag(pardata)))./alpha).')./(mean(abs(real(pskconst))+abs(imag(pskconst)))*subcar/2);

            % calculate theta*data
            td=(pardata./alpha).';

            % estimate qam_data_tx
            psk_data_rx=td./est_theta;

            % serialize qam data
            serdata=reshape(psk_data_rx, subcar*nsym/2,[]).';

            % demodulate the qam data in the ofdm symbols
            bin_data_rx = reshape(de2bi(reshape(pskdemod(serdata, qam),1,[])).',1,[]);

            % calculate SNR and BER
            PBER(n,k)=sum((data(bits)~=bin_data_rx(bits)))/length(bits);



        end

    end

    %% calculate SNR and BER
    EbN0=mean(EsN0)/nbitqam;
    EbN0_dB=10*log10(EbN0);
    PEbN0=mean(PEsN0)/nbitqam;
    PEbN0_dB=10*log10(PEbN0);
    BER=mean(BER);
    PBER=mean(PBER);

qsnr(p,:)=EbN0_dB;
qber(p,:)=BER;
psnr(p,:)=PEbN0_dB;
pber(p,:)=PBER;
end
%% Plot SNR vs. BER for block type CE
figure
semilogy(qsnr(2,:),qber(2,:),'d')
hold on
semilogy(qsnr(3,:),qber(3,:),'s')
hold on
semilogy(psnr(1,:),pber(1,:),'p')
hold on
semilogy(psnr(2,:),pber(2,:),'h')
hold on
semilogy(psnr(3,:),pber(3,:),'^')
xlabel('Eb/N0 (dB)')
ylabel('BER')
legend('16-QAM','64-QAM','QPSK','16-PSK','64-PSK')
grid on

figure
semilogy(tt0,2*(1-pber(1,:)))
hold on
semilogy(tt0,4*(1-pber(2,:)))
hold on
semilogy(tt0,6*(1-pber(3,:)))
hold on
semilogy(tt0,4*(1-qber(2,:)))
hold on
semilogy(tt0,6*(1-qber(3,:)))
xlabel('Time (s)')
ylabel('Average Throughput (bits per subcarrier)')
legend('QPSK/4-PSK','16-PSK','64-PSK','16-QAM','64-QAM')
ylim([0 6.5])
grid on




