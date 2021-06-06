clear all
%close all
clc

load('qamsnrber.mat')

for g=1:7

    v=(g-1)*5; %m/s2

    %qam=2^f;
    nfft=64; %fft size
    subcar=48; %number of data subcarriers

    ncp=4; %cyclic prefix size

    %speed=(g-1)*5;
    
    for f=2:6
        p_se(f-1) = (f*subcar/2)/(nfft+ncp)*(1-qamsnrber{f,g}{1,1}(length(find(qamsnrber{f,g}{2,1}<10))+1));
        b1_se(f-1) = (f*subcar/2)/(nfft+ncp)*1/2*(1-qamsnrber{f,g}{1,2}(length(find(qamsnrber{f,g}{2,2}<10))+1));
        b2_se(f-1) = (f*subcar/2)/(nfft+ncp)*2/3*(1-qamsnrber{f,g}{1,3}(length(find(qamsnrber{f,g}{2,3}<10))+1));
        b3_se(f-1) = (f*subcar/2)/(nfft+ncp)*3/4*(1-qamsnrber{f,g}{1,4}(length(find(qamsnrber{f,g}{2,4}<10))+1));
        c1_se(f-1) = (f*subcar/2)/(nfft+ncp)*1/2*(1-qamsnrber{f,g}{1,5}(length(find(qamsnrber{f,g}{2,5}<10))+1));
        c2_se(f-1) = (f*subcar/2)/(nfft+ncp)*2/3*(1-qamsnrber{f,g}{1,6}(length(find(qamsnrber{f,g}{2,6}<10))+1));
        c3_se(f-1) = (f*subcar/2)/(nfft+ncp)*3/4*(1-qamsnrber{f,g}{1,7}(length(find(qamsnrber{f,g}{2,7}<10))+1));
        st_se(f-1) = (f*subcar/2)/(nfft+ncp)*(1-qamsnrber{f,g}{1,8}(length(find(qamsnrber{f,g}{2,8}<10))+1));
    end
    
    figure
    semilogy(2:6,b1_se,'pb-')
    hold on
    semilogy(2:6,b2_se,'sr-') 
    hold on
    semilogy(2:6,b3_se,'dy-')
    hold on
    semilogy(2:6,c1_se,'+k-')
    hold on
    semilogy(2:6,c2_se,'xc-')
    hold on
    semilogy(2:6,c3_se,'*m-')
    hold on
    semilogy(2:6,st_se,'ob-')
    hold on
    semilogy(2:6,p_se,'^g-')
    legend('Block PS1','Block PS2','Block PS3','Comb PS1','Comb PS2','Comb PS3','ST','Proposed')
    xlabel('log_2M of M-QAM')
    ylabel('Spectral Efficiency (C [bits/s/Hz])')
    title(['Speed ' num2str(v) 'm/s and 10 dB SNR'])
    grid on

end

