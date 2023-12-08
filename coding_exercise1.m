clear;
clear all;
clc;


%GENERATE SYMBOLS( or BITS FOR BPSK)
numBits =10^5;% number of bits
numIter = 100;
txBits = rand(1,numBits)<0.5;%generate bits with rand function


%BPSK MODULATION
modulatedSignal = myModulator(txBits);


%ORTHOGONAL MODULATION
for i = 1:numBits
    if(txBits(i)==0)
        mod(i,1) = sqrt(2);
        mod(i,2) = 0;
    elseif(txBits(i)==1)
        mod(i,1) = 0;
        mod(i,2) = sqrt(2);
    end
end
orthoSignal = mod;


%orthoSignal = orthoModulator(txBits);
% modulating signal by mapping bits to constellation



for settings = [1 2 3]
    %setting 1 for bpsk
    %setting 2 for coherent detection of bpsk
    %setting 3 for non coherent detection of orthogonal signals

    %GENERATING NOISE

    esnodb = 0:2:20;
    for snrIdx = 1:numel(esnodb)

        esno(snrIdx) = 10^(esnodb(snrIdx)/10);

        snr = sqrt(esno(snrIdx));

        for iter = 1:numIter

            if(settings==1||settings==2)
                %rxModSymbols = awgn( modulatedSignal , snr );% adding additive white gaussian noise
                noise = (1/sqrt(2))*(randn(length(modulatedSignal),1) + 1i*randn(length(modulatedSignal),1));
                % generating additive white gaussian noise
                noise = (1/snr).*noise;




            elseif(settings==3 )
                [len, wid] = size(orthoSignal);
                noise = (1/sqrt(2))*(randn(len,wid) + 1i*randn(len,wid));
                % generating additive white gaussian noise
                noise = (1/snr).*noise;
            end

            %MODELLING THE RECEIVED SIGNAL
            if(settings==1)
                h =  1;
                rxModSymbols = h*modulatedSignal + noise;
                % adding noise to modulated signal

            elseif(settings==2)

                h_tilda = (1/sqrt(2.0))*complex(randn(1,1), randn(1,1));
                h =  h_tilda;
                rxModSymbols1= h*modulatedSignal + noise;



                rxModSymbols = (conj(h)/sqrt(abs(h).^2)) *rxModSymbols1;


            elseif(settings==3)
                h_tilda = (1/sqrt(2.0))*complex(randn(1,1), randn(1,1));
                h =  h_tilda;
                rxModSymbols2 = h*orthoSignal + noise;

                %DEMODULATION OF ORTHOGONAL SIGNAL
                for j = 1:len
                    if(abs(rxModSymbols2(j,1))>abs(rxModSymbols2(j,2)))
                        demodulatedSymbols1(j) = 0;
                    else
                        demodulatedSymbols1(j) = 1;

                    end

                end
                demodulatedSymbols = transpose(demodulatedSymbols1);
            end
            %DEMODULATION OF BPSK SIGNAL
            if(settings==1||settings==2)
                [demodulatedSymbols] = myDemodulator(rxModSymbols);
                % demodulating signal by mapping symbols to constellation

                %CALCULATION OF SYMBOL ERROR RATE
                SER =round(abs(modulatedSignal-demodulatedSymbols),2)>0;% symbol error rate
            elseif(settings==3)
                SER =round(abs(txBits'-demodulatedSymbols),2)>0;% symbol error rate
            end
            % BER = round(abs(txBits-demodulatedBits),2);% bit error rate

            numErrors(settings,snrIdx,iter) = sum(SER);% percentage of SER

        end
        SER_percent_avg(settings,snrIdx) = sum(numErrors(settings,snrIdx,:))/(numBits*numIter);
        %       SER_percent_avg1 = squeeze(SER_percent(settings,:,:));
        %       SER_percent_avg(settings) = mean(SER_percent_avg1);
        %      plot(SNRdB,  mean(ProbOfSuccess1));
        %     hold on;


        %THEORITICAL BEHAVIOUR OF EACH SCENARIO
        if(settings==1)
            theoriticalPE(settings,snrIdx) = qfunc(sqrt(2*esno(snrIdx)));
        elseif(settings==2)
            theoriticalPE(settings,snrIdx) = 1/2*(1- sqrt( esno(snrIdx) /(1+ esno(snrIdx) )) );
        elseif(settings==3)
            theoriticalPE(settings,snrIdx) = 1/(2* esno(snrIdx) +2);
        end
    end


    hold on
    semilogy(esnodb,SER_percent_avg(settings,:));

    semilogy(esnodb,theoriticalPE(settings,:));

end

legend('SER percentage: bpsk','TheoriticalPe: bpsk','SER percentage: coherent bpsk',...
    'TheoriticalPe: coherent bpsk','SER percentage: non coherent orthogonal',...
    'TheoriticalPe: non coherent orthogonal')

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
    elseif(txBits1(ii)==1)
        mod(ii,1) = 0;
        mod(ii,2) = 1;
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
