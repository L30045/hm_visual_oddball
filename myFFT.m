function [spec, freq] = myFFT(data, srate, varargin)
if isempty(varargin)
    binwidth = 1;
else
    binwidth = varargin{1};
end

[nCh,L] = size(data);
Y = fft(data,[],2);
P2 = abs(Y/L);
P1 = P2(:,1:L/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
f = srate*(0:(L/2))/L;
tmp_freq = 0:binwidth:(srate/2-binwidth);
bins = discretize(f,tmp_freq+binwidth/2);
spec = zeros(nCh,max(bins));

for b_i = 1:max(bins)
    spec(:,b_i) = sum(P1(:,bins==b_i),2);
end

freq = tmp_freq(2:end);


end