clear;
clear all;
clc;

rng(2023);

%GENERATE SYMBOLS( or BITS FOR BPSK)
numBits = 1024;% number of bits
numIter = 1000;
txBits = rand(1,numBits)<0.5;%generate bits with rand function

trellis = poly2trellis(3,[6 7]);
tbl = 32;
rate = 1/2;

txBits_conv=convenc(txBits,trellis);
%txBits_conv = txBits;

%QPSK MODULATION

constellation = (1/sqrt(2))*[1+1i 1-1i -1+1i -1-1i];% qpsk constellation

reshaped_txBits = transpose(reshape(txBits_conv',[2,length(txBits_conv)/2]));
%we reshape the vector of bits to obtain a vector of symbols of size M each

reshaped_txBits_str = num2str(reshaped_txBits);
%convert from number to string

symNumVec = bin2dec(reshaped_txBits_str);
%maps bits representing symbols to integers

modulatedSignal = constellation(symNumVec+1);
%symbol numbers mapped to constellation


txSym_without_cp = ifft(modulatedSignal.',1024);

%OFDM MODULATION
txSym_with_cp = [];

txSym_with_cp = [txSym_without_cp(769:1024) ; txSym_without_cp];

txWaveform = txSym_with_cp(:); %time domain waveform



%-------------------------------------------------------------------------

%SNR LOOP

esnodb = -60:4:5;
for snrIdx = 1:numel(esnodb)

    esno(snrIdx) = 10^(esnodb(snrIdx)/10);

    snr = sqrt(esno(snrIdx));

    for iter = 1:numIter

        %no = 1/(sqrt(2*snr*1024));
        % generating additive white gaussian noise
        noise = (1/(snr*sqrt(2)*1024))*(randn(length(txWaveform),1) + 1i*randn(length(txWaveform),1));

        %modelling N tap channel
        N_taps = 1; % number of taps
        h = (1/sqrt(2))*(randn(1,N_taps) + 1i*randn(1,N_taps));
        H = fft(h,1024);% estimation in time domain

        % Model to the received time domain waveform
        
        rxWaveform1 = conv(h,txWaveform) ;
        rxWaveform2 = rxWaveform1(1:1280) +  noise;
        rxWaveform = rxWaveform2(1:1280);




        %OFDM DEMODULATION

        td_signal_without_cp = rxWaveform(257:1152);

        y = fft(td_signal_without_cp,1024);

        y_tilda = (conj(H).*y.')./abs(H);


        %QPSK DEMODULATION
        y1(1:2:2048) = real(y_tilda);
        y1(2:2:2048) = imag(y_tilda);

        % fd_rx2  = (reshape(y1.',[numRep],[])).';
        % fd_rx3 = mean(fd_rx2,2); % soft decoding
       


        demodBits = [];
        for ii = 1: length(y1)
            if(y1(ii)>0)
                demodBits = [demodBits 0];
            else
                demodBits = [demodBits 1];
            end
        end


         %decoded = demodBits;
         decoded = vitdec(demodBits,trellis,tbl,'trunc','hard');
        % decoded = vitdec((round(y1)+64),trellis,tbl,'trunc','soft',127 );

        %     distances = abs(real(constellation)-real(fd_rx5(ii)))+abs(imag(constellation)-imag(fd_rx5(ii)));
        %     % distances is the norm of received symbols and each constellation
        %     [~ , index] = min(distances,[],2)   ;
        %     %min of distances gives the index of the most appropriate constellation
        %     detModSym(ii) = constellation(index);
        % 
        % 
        % end
        % 
        % demodulatedSymbols = detModSym.';
        % demodBits = [];
        % for ii = 1: length(demodulatedSymbols)
        % 
        %     if(demodulatedSymbols(ii) == constellation(1))
        %         demodBits = [demodBits 0 0];
        %     elseif(demodulatedSymbols(ii) == constellation(2))
        %         demodBits = [demodBits 0 1];
        %     elseif(demodulatedSymbols(ii) == constellation(3))
        %         demodBits = [demodBits 1 0];
        %     elseif(demodulatedSymbols(ii) == constellation(4))
        %         demodBits = [demodBits 1 1];
        %     end
        % end

        
      

        %PROBABILITY OF ERROR CALCULATION

        BER = round(abs(txBits-decoded),2);% bit error rate

        numErrors(snrIdx,iter) = sum(BER);% percentage of BER


    end
    BER_percent_avg(snrIdx) = sum(numErrors(snrIdx,:))/(numBits*numIter);




end


plot(esnodb,BER_percent_avg(:));
title('BER percentage vs SNR');
