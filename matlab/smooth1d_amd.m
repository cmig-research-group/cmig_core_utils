function [vec_out] = smooth1d_amd(vec_in,sig1)

idim = length(vec_in);
ivec = [-floor(idim/2):ceil(idim/2)-1];
ivec = ifftshift(ivec);
sfiltvec = exp(-1/2*(ivec/sig1).^2); sfiltvec = sfiltvec/sum(sfiltvec);
filtvec = fft(sfiltvec);
if (size(filtvec,1) ~= size(vec_in,1))
  filtvec = filtvec.';
end
kvec = fft(ifftshift(vec_in));
vec_out = fftshift(ifft(kvec.*filtvec));

% ToDo
%   Default to old version of smooth1d, allow for additional argument to turn on correct scaling of sig1
%   Should use analytic form of Fourier transform of Gaussian

%keyboard

% Should permit zero-padding


% yvec = zeros(1,512); yvec(1)=1; xvec = 0:length(yvec)-1; sig = 100; yvec = smooth1d_amd(yvec,sig); yvec = yvec/yvec(1); 
% figure(1000); plot(xvec,yvec,xvec,exp(-4*(xvec/sig).^2))
