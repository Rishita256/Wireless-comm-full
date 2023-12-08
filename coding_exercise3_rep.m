clear;
clear all;
clc;

rng(2023);

%GENERATE SYMBOLS( or BITS FOR BPSK)
for numRep = [2 4]; %number of repetitions
numBits = 1024*2/numRep;% number of bits
numIter = 1000; %number of iterations
txBits = rand(1,numBits)<0.5;% generate bits with rand function


txBits_rep = repelem(txBits, numRep);

%QPSK MODULATION
modulatedSignal = myQPSKmod(txBits_rep);


%time domain signal without cp
txSym_without_cp = ifft(modulatedSignal.',1024);

%OFDM MODULATION
txSym_with_cp = [];
%cp addition
txSym_with_cp = [txSym_without_cp(769:1024) ; txSym_without_cp];

txWaveform = txSym_with_cp(:); %time domain waveform



%-------------------------------------------------------------------------

%SNR LOOP

esnodb = -60:4:5;
for snrIdx = 1:numel(esnodb)

    esno(snrIdx) = 10^(esnodb(snrIdx)/10);

    snr = sqrt(esno(snrIdx));

    for iter = 1:numIter

        % generating additive white gaussian noise
        noise = (1/sqrt(snr*2*1024))*(randn(length(txWaveform),1) + 1i*randn(length(txWaveform),1));

        %modelling N tap channel
        N_taps = 2; % number of taps
        h = (1/sqrt(2))*(randn(1,N_taps) + 1i*randn(1,N_taps));
        H = fft(h,1024);% estimation in time domain

        % Model to the received time domain waveform
        rxWaveform1 = conv(h,txWaveform) ; %(multiplication in freq domain)
        rxWaveform2 = rxWaveform1(1:1280) +  noise;
        rxWaveform = rxWaveform2(1:1280);


        %OFDM DEMODULATION
        %cp removal
        td_signal_without_cp = rxWaveform(257:1152);
        %conversion to frequency domain
        y = fft(td_signal_without_cp,1024);

        y_tilda = (conj(H).*y.')./abs(H);


        %QPSK DEMODULATION
        demodBits = myQPSKdemod(y_tilda, numRep);

        
        %PROBABILITY OF ERROR CALCULATION

        BER = round(abs(txBits-demodBits),2);% bit error rate

        numErrors(snrIdx,iter) = sum(BER);% percentage of BER


    end
    BER_percent_avg(snrIdx) = sum(numErrors(snrIdx,:))/(numBits*numIter);




end


plot(esnodb,BER_percent_avg(:));
hold on
title('BER percentage vs SNR');

end

legend('Code rate = 1/2', 'Code rate = 1/4')

%% 


function modulatedSignalf = myQPSKmod(txBits_repf)
%QPSK MODULATION

constellation = (1/sqrt(2))*[1+1i 1-1i -1+1i -1-1i];% qpsk constellation

reshaped_txBits = transpose(reshape(txBits_repf',[2,length(txBits_repf)/2]));
%we reshape the vector of bits to obtain a vector of symbols of size M each

reshaped_txBits_str = num2str(reshaped_txBits);
%convert from number to string

symNumVec = bin2dec(reshaped_txBits_str);
%maps bits representing symbols to integers

modulatedSignalf = constellation(symNumVec+1);
%symbol numbers mapped to constellation
end


function demodBitsf = myQPSKdemod(y_tildaf, numRep)
y1(1:2:2048) = real(y_tildaf);
        y1(2:2:2048) = imag(y_tildaf);

        fd_rx2  = (reshape(y1.',[numRep],[])).';
        fd_rx3 = mean(fd_rx2,2); % soft decoding
       


        demodBitsf = [];
        for ii = 1: length(fd_rx3)
            if(fd_rx3(ii)>0)
                demodBitsf = [demodBitsf 0];
            else
                demodBitsf = [demodBitsf 1];
            end
        end
end
      
