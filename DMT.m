clear all;
close all;

format long;

%   ---------------
%   A: Setting Parameters
%   ---------------
M = 16;                          %   Qam signal constellation
% Modulation order
k = log2(M);     % Bits per symbol
EbNoVec = (5:15)';      % Eb/No values (dB)
numSymPerFrame = 100;   % Number of QAM symbols per frame

no_of_data_points = 144;        %   have 128 data points
block_size = 8;                 %   size of each ofdm block
cp_len = ceil(0.1*block_size);  %   length of cyclic prefix
no_of_ifft_points = block_size;           %   128 points for the FFT/IFFT
no_of_fft_points = block_size;


n = 3e4; % Number of bits to process
nsamp = 1; % Oversampling rate


berEst = zeros(size(EbNoVec));


%   ---------------------------------------------
%   B:  %   +++++   TRANSMITTER    +++++
%   ---------------------------------------------
%   1.  Generate 1 x 128 vector of random data points
data_source = randsrc(1, no_of_data_points, 0:M-1);

%data_source = randi([0 1],1000*k,1);
%data_source=randi(9600,1);
%data_source=randi([0 M-1],1000,1);
figure(1)
stem(data_source); grid on; xlabel('Data Points'); ylabel('transmitted data phase representation')
title('Transmitted Data "O"')


%   2.  Perform QAM modulation
qam_modulated_data = qammod(data_source,M);
scatterplot(qam_modulated_data);title('MODULATED TRANSMITTED DATA');

%   3.  Do IFFT on each block
%   Make the serial stream a matrix where each column represents a pre-OFDM
%   block (w/o cyclic prefixing)
%   First: Find out the number of colums that will exist after reshaping
num_cols=(length(qam_modulated_data)/block_size);
data_matrix = reshape(qam_modulated_data, block_size, num_cols);
%   Second: Create empty matix to put the IFFT'd data
cp_start = block_size-cp_len;
cp_end = block_size;
%   Third: Operate columnwise & do CP
for i=1:num_cols
    ifft_data_matrix(:,i) = ifft((data_matrix(:,i)),no_of_ifft_points);
    %   Compute and append Cyclic Prefix
    for j=1:cp_len
       actual_cp(j,i) = ifft_data_matrix(j+cp_start,i);
    end
    %   Append the CP to the existing block to create the actual OFDM block
    ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_data_matrix(:,i));
end


%   4.  Convert to serial stream for transmission

[rows_ifft_data, cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;
%   Actual OFDM signal to be transmitted
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
figure(3)
plot(real(ofdm_signal)); xlabel('Time'); ylabel('Amplitude');
title('OFDM Signal');grid on;


%ofdmDims = info(ofdm_signal);
%numDC =ofdm_signal.DataInputSize(1);



%   ---------------------------------------------------------------
%     +++++   clipping as a PAPR reduction method    +++++
%   ---------------------------------------------------------------
%  avg=0.9;
%  clipped=ofdm_signal;
%  for i=1:length(clipped)
%  	if clipped(i) > avg
%  		clipped(i) = avg;
%      end
%      if clipped(i) < -avg
%  		clipped(i) = -avg;
%      end
%  end
%  figure(4)
%  plot(real(clipped)); xlabel('Time'); ylabel('Amplitude');
%  title('clipped Signal');grid on;

 avg=0;
 clipped=ofdm_signal;
 for i=1:length(clipped)
 	if clipped(i) > avg
 		clipped(i) = clipped(i);
     
    elseif clipped(i) < -avg

 		clipped(i) = 0;
  
    end
   
%     figure(4)
% plot(real(clipped)); xlabel('Time'); ylabel('Amplitude');
%  title('clipped Signal');grid on;
% axis([0 180 -1.5 1.5]);
 end
 %figure(4)
% plot(real(clipped)); xlabel('Time'); ylabel('Amplitude');
 %title('clipped Signal');grid on;
figure(4)
plot(real(clipped)); xlabel('Time'); ylabel('Amplitude');
 title('clipped Signal');grid on;
axis([0 180 -1.5 1.5]);




%   --------------------------------
%     +++++   CHANNEL    +++++
%   --------------------------------

%   Create a complex multipath channel
channel = randn(1,block_size) + sqrt(-1)*randn(1,block_size);



%   ------------------------------------------
%     +++++   RECEIVER    +++++
%   ------------------------------------------

%   1.  Pass the ofdm signal through the channel

after_channel = filter(channel, 1, ofdm_signal);


%SNR with respect to ber
EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(nsamp);
receivedSignal = awgn(qam_modulated_data,snr,'measured');
sPlotFig = scatterplot(receivedSignal,1,0,'g.');hold on
sPlotFig1 = scatterplot(after_channel,1,0,'r.');hold on

figure(5)
scatterplot(qam_modulated_data,1,0,'k*',sPlotFig)
scatterplot(qam_modulated_data,1,0,'k*',sPlotFig1)


%   2.   Add Noise

awgn_noise = awgn(zeros(1,length(after_channel)),0);

%   3.  Add noise to signal...

recvd_signal = awgn_noise+after_channel;

%   4.  Convert Data back to "parallel" form to perform FFT

recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);

%   5.  Remove CP

recvd_signal_matrix(1:cp_len,:)=[];

%   6.  Perform FFT

for i=1:cols_ifft_data
    
    %   FFT
    
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);

end



%   7.  Convert to serial stream

recvd_serial_data = reshape(fft_data_matrix,1,(block_size*num_cols));

%   8.  Demodulate the data

qam_demodulated_data = qamdemod(recvd_serial_data,M);

figure(7)

stem(qam_demodulated_data,'rx');

grid on;xlabel('Data Points');ylabel('received data phase representation');title('Received Data "X"')   

dataOutMatrix = de2bi(qam_demodulated_data,k);
dataOut = dataOutMatrix(:);

[numErrors,ber] = biterr(data_source,dataOut);
fprintf('\nThe binary coding bit error rate = %5.2e, based on %d errors\n',ber,numErrors)


%   ----------------------------------------------------
%   F:  %   +++++   RECEIVER of clipped signal    +++++
%   ----------------------------------------------------

%   1.  Pass the ofdm signal through the channel

after_channel = filter(channel, 1, clipped);

%   2.   Add Noise

awgn_noise = awgn(zeros(1,length(after_channel)),0);

%   3.  Add noise to signal...

recvd_signal = awgn_noise+after_channel;

%   4.  Convert Data back to "parallel" form to perform FFT

recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);

%   5.  Remove CP

recvd_signal_matrix(1:cp_len,:)=[];

%   6.  Perform FFT

for i=1:cols_ifft_data

    %   FFT
    
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);

end

%   7.  Convert to serial stream

recvd_serial_data = reshape(fft_data_matrix, 1,(block_size*num_cols));

%   8.  Demodulate the data

qam_demodulated_data = qamdemod(recvd_serial_data,M);

figure(8)

stem(qam_demodulated_data,'rx');

%BER
 noe=0; 
 for ii= 1:1:length(recvd_serial_data)
   if data_source(ii)~=recvd_serial_data(ii)
       
    noe=noe+1;
   else
       noe=noe;
   
   end 
 end 
% ber= noe/length(ofdm_signal); 
 
ber=( noe/length(data_source)); 
 
% 
% for n = 1:length(EbNoVec)
%     % Convert Eb/No to SNR
%     snrdB = EbNoVec(n) + 10*log10(k);
%     % Reset the error and bit counters
%     numErrs = 0;
%     numBits = 0;
%     
%     while numErrs < 200 && numBits < 1e9
%         % Generate binary data and convert to symbols
%         dataIn = randi([0 1],numSymPerFrame,k);
%         dataSym = bi2de(dataIn);
%         
%         % QAM modulate using 'Gray' symbol mapping
%         txSig = qammod(dataSym,M);
%         
%         % Pass through AWGN channel
%         rxSig = awgn(txSig,snrdB,'measured');
%         
%         % Demodulate the noisy signal
%         rxSym = qamdemod(rxSig,M);
%         % Convert received symbols to bits
%         dataOut = de2bi(rxSym,k);
%         
%         % Calculate the number of bit errors
%         nErrors = biterr(dataIn,dataOut);
%         
%         % Increment the error and bit counters
%         numErrs = numErrs + nErrors;
%         numBits = numBits + numSymPerFrame*k;
%     end
%     
%     % Estimate the BER
%     berEst(n) = numErrs/numBits;
% end
% berTheory = berawgn(EbNoVec,'qam',M);
% 
% semilogy(EbNoVec,berEst,'*')
% hold on
% semilogy(EbNoVec,berTheory)
% grid
% legend('Estimated BER','Theoretical BER')
% xlabel('Eb/No (dB)')
% ylabel('Bit Error Rate')

%http://www.raymaps.com/index.php/ber-64-qam-awgn/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE 64-QAM BER USING FORMULA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbNodB=0:2:16;
EbNo=10.^(EbNodB/10);
x=sqrt(3*k*EbNo/(M-1));
Pb=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));
figure(9) 
semilogy(EbNodB,Pb,'bs-')
title('QAM bit error rate');
xlabel('EbNo');
ylabel('Pb');

sinad(ofdm_signal,after_channel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%https://www.edaboard.com/showthread.php?143238-matlab-code-for-SNR-vs-BER-plot

% %% Signal Source
% % Create a binary data stream as a column vector.
% x = data_source; % Random binary data stream
% 
% % Plot first 40 bits in a stem plot.
% % stem(x(1:70),'filled');
% % title('Random Bits');
% % xlabel('Bit Index'); ylabel('Binary Value');
% %% Bit-to-Symbol Mapping
% % Convert the bits in x into k-bit symbols.
% xsym = bi2de(reshape(x,k,length(x)/k).','left-msb');
% 
% %% Stem Plot of Symbols 
% % Plot first 10 symbols in a stem plot.
% figure; % Create new figure window.
% stem(xsym(1:10));
% title('Random Symbols');
% xlabel('Symbol Index'); ylabel('Integer Value');
% 
% %% Modulation
% % Modulate using 16-QAM.
% y = qam_modulated_data;
% %% Transmitted Signal
% ytx = y;
% 
% %% Channel
% % Send signal over an AWGN channel.
% EbNo = 10; % In dB
% snr = EbNo + 10*log10(k) - 10*log10(nsamp);
% ynoisy = awgn(ytx,snr,'measured');
% 
% %% Received Signal
% yrx = ynoisy;
% 
% %% Scatter Plot
% % Create scatter plot of noisy signal and transmitted
% % signal on the same axes.
% h = scatterplot(yrx(1:nsamp*5e1),nsamp,3,'g.');
% hold on;
% scatterplot(ytx(1:5e1),1,0,'k*',h);
%  title('Received Signal');
% legend('Received Signal','Signal Constellation');
%  axis([-5 5 -5 5]); % Set axis ranges.
%  hold off
% % Demodulation
% % Demodulate signal using 16-QAM.
% zsym =qam_demodulated_data;
% %% Symbol-to-Bit Mapping
% % Undo the bit-to-symbol mapping performed earlier.
% z = de2bi(zsym,'left-msb'); % Convert integers to bits.
% % Convert z from a matrix to a vector.
% z = reshape(z.',numel(z),1);
% %% BER Computation
% % Compare x and z to obtain the number of errors and
% % the bit error rate.
% [errors,errorrate] = biterr(x,z);
% hold on
% semilogy(snr,errorrate);
% legend('snr','errorrate');
% title(' BER OF QAM');
% hold off


  %https://dsp.stackexchange.com/questions/16240/deriving-ser-ber-for-4qam-16qam-and-32qam
  
 % https://www.mathworks.com/matlabcentral/fileexchange/39011-ber-comparison-of-m-ary-qam
  
  %https://dsp.stackexchange.com/questions/42857/converting-e-b-n-0-to-snr