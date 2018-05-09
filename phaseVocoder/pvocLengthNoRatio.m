function [y, lastph] = pvocLengthNoRatio(x, targetL, n, ph)
% y = pvoc(x, r, n)  Time-scale a signal to r times faster with phase vocoder
%      x is an input sound. n is the FFT size, defaults to 1024.  
%      Calculate the 25%-overlapped STFT, squeeze it by a factor of r, 
%      inverse spegram.
% 2000-12-05, 2002-02-13 dpwe@ee.columbia.edu.  Uses pvsample, stft, istft
% $Header: /home/empire6/dpwe/public_html/resources/matlab/pvoc/RCS/pvoc.m,v 1.3 2011/02/08 21:08:39 dpwe Exp $

if nargin < 3
  n = 1024;
end

% Effect of hanns at both ends is a cumulated cos^2 window (for
% r = 1 anyway); need to scale magnitudes by 2/3 for
% identity input/output
%scf = 2/3;
% 2011-02-07: this factor is now included in istft.m
scf = 1.0;
hop=n/4;
% Calculate the basic STFT, magnitude scaled
X = scf * stft(x', n, n, hop);


% Calculate the new timebase samples
[rows, cols] = size(X);
% With hann windowing on both input and output, 
% we need 25% window overlap for smooth reconstruction

lengthT=ceil((targetL-n)/hop +1);
t=(0:1/(lengthT-1):1)*(cols-2);
% Have to stay two cols off end because (a) counting from zero, and 
% (b) need col n AND col n+1 to interpolate

%length(y) = hop*(length(t)-1)+n
% (targetL-n)/hop +1 = length(t)
%difference=targetL-(hop*(length(t)-1)+n);
hopV= round((0: 1/(length(t)-1):1)*(targetL-n));
hopV(1)=1;
hopV(end)=targetL-n;

% Generate the new spectrogram
if (exist('ph','var'))
    X2 = pvsampleHopV(X, t, hopV, ph);
else
    X2 = pvsampleHopV(X, t, hopV);
end

lastph=angle(X2(:,end));

% Invert to a waveform
y = istftHopV(X2, n, n, hopV)';
