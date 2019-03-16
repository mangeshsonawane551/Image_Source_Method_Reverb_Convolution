%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%Author: Mangesh Sonawane 
%Program Details: An implementation of convolution operation that 
%convolves an input .wav file  with the impulse response simulated
%using the ?ImageSource...m? script, in order to make the sound in the 
%input .wav file appear to have been created and recorded in a virtualroom.
%Convolution has been implemented using time domain as well as frequency
%domain( fast convolution theorem).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


clc;
clear all;



%-------------------------------------------------------------------------%
                    %Reads input dry audio file
%-------------------------------------------------------------------------%
[x,Fs] = audioread('castanetcut_44.wav');

% Checks if stereo and converts to mono
if size(x,2) == 2                            %detects if stereo or mono
    x = sum(x,2) / size(x,2);                %convert to mono
end


%-------------------------------------------------------------------------%
                    %Reads input impulse response wave file
%-------------------------------------------------------------------------%
[impulse_response,Fs] = audioread('IR_10.8X13.5X20.3_WET_S1889125_Sonawane.wav');

% Extracts only left channel audio
if size(x,2) > 1
    impulse_response = impulse_response(:,1); 
end

% Length of dry audio 
Lx = length(x);

%Length of impulse response sound
Li = length(impulse_response);

%Summation of total lenght
N = Lx + Li;

%Initialising output vector
y = zeros(N,1);

%-------------------------------------------------------------------------%
          %Frequency domain convolution usinf fast convolution theorem
%-------------------------------------------------------------------------%
%Calculating FFT of input audio
fftx = fft(x,N);

%Calculating FFT of Impulse response
ffth = fft(impulse_response,N);

%Output 'y'is the inverse fft of dot multiplication on input audio and 
%impulse response
y = ifft(fftx.*ffth);


%-------------------------------------------------------------------------%
                        %Time domain convolution (UNCOMMENT TO USE)
%-------------------------------------------------------------------------%
% %Zero padding initial and end on input audio vector
% x = [zeros(Li,1); x]; 
% x(end+Li) = 0; 
% 
% %Transpose of input audio vector
% x = transpose(x);
% 
% %Transpose and flip array from left to right
% impulse_response = transpose(fliplr(impulse_response'));
% 
% %Start stop watch
% tic
% for n = 1:N-1
%     y(n) = x(n+1:n+Li)*impulse_response(1:Li);
% end
% 
% %End stop watch
% toc

%-------------------------------------------------------------------------%

%Normalise
y = y/max(abs(y));

%Output convolved audio
soundsc(y,Fs);



