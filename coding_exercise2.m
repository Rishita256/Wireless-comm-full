clear;
clear all;
clc;


%GENERATE SYMBOLS( or BITS FOR BPSK)
numBits = 10^3;% number of bits
numIter = 1000;
txBits = (rand(1,numBits)<0.5)';%generate bits with rand function

 numRepetitions = [1:5];

 for numIdx = 1:numel(numRepetitions)
     numRep = numRepetitions(numIdx);
    txBits_encoded = repelem(txBits, numRep);


    %BPSK MODULATION
    modulatedSignal = myModulator(txBits_encoded);
    modulatedSignalRef = myModulator(txBits);


    %GENERATING NOISE
    esnodb = 0:1:10;
    for snrIdx = 1:numel(esnodb)

        esno(snrIdx) = 10^(esnodb(snrIdx)/10);

        snr = sqrt(esno(snrIdx));

        for iter = 1:numIter

            noise = (1/sqrt(2))*(randn(length(modulatedSignal),1) + 1i*randn(length(modulatedSignal),1));
            % generating additive white gaussian noise
            noise = (1/snr).*noise;
           

            %MODELLING THE RECEIVED SIGNAL

           % for transmitterIdx = 1:numTransmitters
                h_tilda = (1/sqrt(2.0))*complex(randn(length(modulatedSignal),1), randn(length(modulatedSignal),1));
                h =  h_tilda;
                rxModSymbols1 = h.*modulatedSignal ;
           % end
            %OBTAINING Y_TILDA
            rxModSymbols2 =rxModSymbols1 + noise; %noise modelled at receiver

            rxModSymbols3 = transpose(reshape(rxModSymbols2, [numRep,length(txBits)]));
            h3 = transpose(reshape(h, [numRep,length(txBits)]));
            
            norm_h = sqrt(sum(abs(h3).^2,2));
            
            for ii = 1:1000
                fac = conj(h3(ii,:))./norm_h(ii);
                rxModSymbols(ii) = fac*(rxModSymbols3(ii,:).');
            end
            
           

            % h2 = sum(h);% channel estimate 
            % norm_h = sqrt(sum(abs(h).^2));
            % 
            % rxModSymbols = (conj(h2)/norm_h) *rxModSymbols2;
            % 
            % 

            %DEMODULATION OF BPSK SIGNAL

            [demodulatedSymbols] = myDemodulator(rxModSymbols);
            % demodulating signal by mapping symbols to constellation


            %CALCULATION OF SYMBOL ERROR RATE
            SER =round(abs(modulatedSignalRef-demodulatedSymbols),2)>0;% symbol error rate


            numErrors(snrIdx,iter) = sum(SER);% percentage of SER

        end
        SER_percent_avg(snrIdx, numRep) = sum(numErrors(snrIdx,:))/(numBits*numIter);
        %       SER_percent_avg1 = squeeze(SER_percent(settings,:,:));
        %       SER_percent_avg(settings) = mean(SER_percent_avg1);
        %      plot(SNRdB,  mean(ProbOfSuccess1));
        %     hold on;


    end


    hold on
    plot(esnodb,SER_percent_avg(:,numRep));



 end

legend('SER percentage: L = 1','SER percentage: L = 3','SER percentage: L = 6')



%% Functions

function mod = myModulator(txBits)
%bits are input to the function and they need to be mapped to the provided
%constellation
numBits = length(txBits);%number of bits in the binary vector
constellation = [-1 1];% bpsk constellation

modSymbols = constellation(txBits+1);
%symbol numbers mapped to constellation
mod = modSymbols';
end


function mod = orthoModulator(txBits1)

numBits1 =length(txBits1);% number of bits
txBits1 = rand(1,numBits1)<0.5;%generate bits with rand function

for ii = 1:numBits1
    if(txBits1(ii)==0)
        mod(ii,1) = 1;
        mod(ii,2) = 0;
        %mod(ii,:) = [1 0]
    elseif(txBits1(ii)==1)
        mod(ii,1) = 0;
        mod(ii,2) = 1;
        %mod(ii,:) = [0 1]
    end
end

end

function [demodSym, demodBits] = myDemodulator(rxModSymbols)
%symbols are input to the function and they need to be mapped to the provided
%constellation to give bits
constellation = [-1 1];% bpsk constellation


for ii = 1: length(rxModSymbols)
%     distances = abs(constellation-real(rxModSymbols(ii)));
%     % distances is the norm of received symbols and each constellation
%     [~ , index] = min(distances)   ;
%     %min of distances gives the index of the most appropriate constellation

    if(real(rxModSymbols(ii))>0)
        index = 2;
    else
        index = 1;
    end

    detModSym(ii) = constellation(index);
end

demodSym = detModSym';

end
