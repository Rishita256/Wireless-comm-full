clear;
clear all;
clc;

rng(2023);
%GENERATE SYMBOLS( or BITS FOR BPSK)
numBits = 1000;% number of bits
numIter = 1000;% number of iterations
txBits = (rand(1,numBits)<0.5)';%generate bits with rand function

%QPSK MODULATION

modulatedSignal = myQPSKmod(txBits);


alamouti_signal = [];
%ALAMOUTI
for antenna = [1]

    for idx = 1: length(modulatedSignal)
        if(rem(idx,2) == 1)
            alamouti_signal(idx,antenna) = modulatedSignal(idx);
        else
            alamouti_signal(idx, antenna) = -conj(modulatedSignal(idx));
        end
    end
end

for antenna = [2]

    for idx = 1: length(modulatedSignal)
        if(rem(idx,2) == 1)
            alamouti_signal(idx,antenna) = modulatedSignal(idx+1);
        else
            alamouti_signal(idx, antenna) = conj(modulatedSignal(idx-1));
        end
    end
end

txWaveform = alamouti_signal;

%-------------------------------------------------------------------------

%SNR LOOP
numIter = 1000; %number of iterations for each snr value
esnodb = -20:5:20;
for snrIdx = 1:numel(esnodb)

    esno(snrIdx) = 10^(esnodb(snrIdx)/10);

    snr = sqrt(esno(snrIdx));

    for iter = 1:numIter

        noise = (1/sqrt(2))*(randn(length(txWaveform),1) + 1i*randn(length(txWaveform),1));
        % generating additive white gaussian noise
        noise = (1/snr).*noise;

        for antenna = [1 2]
            %modelling channel
            h(antenna) = (1/sqrt(2.0))*complex(randn(1,1), randn(1,1));

            % Add AWGN to the received time domain waveform
            rxWaveform(:,antenna) = h(antenna)*txWaveform(:,antenna) ;
        end
rxWaveform_final = sum(rxWaveform,2)+noise;

        %ALAMOUTI DECODING
for y_idx = 1:length(rxWaveform_final)
    if(rem(y_idx,2)==1)
        u(y_idx) = [conj(h(1)) h(2)]/sqrt(abs(h(1)^2+h(2)^2))*[rxWaveform_final(y_idx) ; conj(rxWaveform_final(y_idx+1))];
    else
        u(y_idx) = [conj(h(2)) -h(1)]/sqrt(abs(h(1)^2+h(2)^2))*[rxWaveform_final(y_idx-1); conj(rxWaveform_final(y_idx))];

    end
end
    rx_decoded = transpose(u);
         
        %QPSK DEMODULATION

       constellation = [1+1i 1-1i -1+1i -1-1i];% qpsk constellation
        for ii = 1: length(rx_decoded)
            distances = abs(constellation-rx_decoded);
            % distances is the norm of received symbols and each constellation
            [~ , index] = min(distances,[],2)   ;
            %min of distances gives the index of the most appropriate constellation
            %detModSym = constellation(index);

            %detModBit(N*i-N+1 :N*i) = dec2bin(index-1, N);
            detModBit = double(dec2bin(index-1, 2))-48;

        end
       
        %demodulatedSymbols = detModSym;
        demodBits = reshape(detModBit', [], 1);

    

        %PROBABILITY OF ERROR CALCULATION

        BER = round(abs(txBits-demodBits),2);% bit error rate
   
        numErrors(iter,snrIdx) = sum(BER);% percentage of BER

    
   
    BER_percent_avg(snrIdx) = mean(numErrors(:,snrIdx))/numBits;
    end
end

plot(esnodb,BER_percent_avg(:));
title('BER percentage vs SNR');

%% 


function modulatedSignalf = myQPSKmod(txBitsf)
constellation = [1+1i 1-1i -1+1i -1-1i];% qpsk constellation

reshaped_txBits = transpose(reshape(txBitsf,[2,length(txBitsf)/2]));
%we reshape the vector of bits to obtain a vector of symbols of size M each

reshaped_txBits_str = num2str(reshaped_txBits);
%convert from number to string

symNumVec = bin2dec(reshaped_txBits_str)+1;
%maps bits representing symbols to integers

modulatedSignalf = transpose(constellation(symNumVec));
%symbol numbers mapped to constellation
end
