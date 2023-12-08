clear;
clear all;
clc;

rng(2023);

Nc = 128;
power = 1;
total_power = power*Nc;
threshold = 0.01;

%set lambda
lambda = 4;

%generate h and sort
h_tilda = (1/sqrt(2.0))*complex(randn(Nc,1), randn(Nc,1));
h_tilda_sorted = sort(h_tilda);

SNRdB = 5;


for snrIdx = 1:numel(SNRdB)

    SNR(snrIdx) = 10^(SNRdB(snrIdx)/10);

    %generate No
    No = 1/(sqrt(2*SNR(snrIdx))).*ones(Nc, 1);
    % % generating additive white gaussian noise
    % noise = no*(randn(length(txWaveform),1) + 1i*randn(length(txWaveform),1));

    %compute No/|h|^2
    var =  No./abs(h_tilda_sorted).^2;

    %find capacity
i =0;
    while i == 0
        %compute total power
        power_iter = max(((1/lambda) - var),0);
        power_iter_sum = sum(power_iter);

        %reset lambda
        if( (total_power-power_iter_sum) < threshold)
            i=1;
        else
            i=0;
            lambda = lambda-threshold/10;
        end

    end
end

if ((total_power-power_iter_sum) < 0)
    lambda = lambda+threshold/10;
     power_iter_final = max(((1/lambda) - var),0);
     power_iter_sum_final = sum(power_iter_final);
end



capacity = sum(log(1+power_iter_final./var));

capacity_without_waterfilling = sum(log(1+power./var));
