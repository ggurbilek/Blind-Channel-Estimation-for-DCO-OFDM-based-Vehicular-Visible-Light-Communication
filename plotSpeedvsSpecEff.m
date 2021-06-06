clear all
%close all
clc

load('qamsnrber.mat')

for f=2:6

    g=1:7;
    v=(g-1)*5; %m/s2

    qam=2^f;
    nfft=64; %fft size
    subcar=48; %number of data subcarriers
    nbitqam = log2(qam);

    ncp=4; %cyclic prefix size
    nbitsym = nfft + ncp;

    const = (nbitqam*subcar/2)/(nfft+ncp);

    speed=(g-1)*5;
    
    for a=1:7
        p_se(a) = qamsnrber{f,a}{1,1}(length(find(qamsnrber{f,a}{2,1}<10))+1);
        b1_se(a) = qamsnrber{f,a}{1,2}(length(find(qamsnrber{f,a}{2,2}<10))+1);
        b2_se(a) = qamsnrber{f,a}{1,3}(length(find(qamsnrber{f,a}{2,3}<10))+1);
        b3_se(a) = qamsnrber{f,a}{1,4}(length(find(qamsnrber{f,a}{2,4}<10))+1);
        c1_se(a) = qamsnrber{f,a}{1,5}(length(find(qamsnrber{f,a}{2,5}<10))+1);
        c2_se(a) = qamsnrber{f,a}{1,6}(length(find(qamsnrber{f,a}{2,6}<10))+1);
        c3_se(a) = qamsnrber{f,a}{1,7}(length(find(qamsnrber{f,a}{2,7}<10))+1);
        st_se(a) = qamsnrber{f,a}{1,8}(length(find(qamsnrber{f,a}{2,8}<10))+1);
    end
    
    figure
    semilogy(speed,const*1/2*(1-b1_se),'pb-')
    hold on
    semilogy(speed,const*2/3*(1-b2_se),'sr-') 
    hold on
    semilogy(speed,const*3/4*(1-b3_se),'dy-')
    hold on
    semilogy(speed,const*1/2*(1-c1_se),'+k-')
    hold on
    semilogy(speed,const*2/3*(1-c2_se),'xc-')
    hold on
    semilogy(speed,const*3/4*(1-c3_se),'*m-')
    hold on
    semilogy(speed,const*(1-st_se),'ob-')
    hold on
    semilogy(speed,const*(1-p_se),'^g-')
    legend('Block PS1','Block PS2','Block PS3','Comb PS1','Comb PS2','Comb PS3','ST','Proposed')
    xlabel('Relative Speed (m/s)')
    ylabel('Spectral Efficiency (C [bits/s/Hz])')
    title([num2str(qam) ' QAM'])
    grid on

end

