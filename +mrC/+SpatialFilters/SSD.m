function [OutAxx,W,A,D] = SSD(InAxx,signal_freqs,varargin)

% This function calculates finds a spatial filter based on the the SSD of
% InAxx. The signal is assumed in freqs, noise is assumed in the
% neighboring bins

% INPUT:
    % InAxx: EEG data in Axx format
    % freqs: Frequnecies of signal considered
%   <options>:
    % do_whitening: Apply whitening and dimensionality deflation prior to
    %               estimation of spatial filter. Default: true.
    % rank_ratio:   Threshold for dimensionality reduction based on 
    %               eigenvalue spectrum. Default: 10^-4
% OUTPUT:
    % OutAxx: Data in component space in Axx format
    % W: Spatial filter
    % A: Activation pattern
% 
% Written by Sebastian Bosse, 3.8.2018

opt	= ParseArgs(varargin,...
    'do_whitening', true, ...
    'rank_ratio', 10^-4 ...             
    );


freqs = [0:InAxx.nFr-2,InAxx.nFr-1:-1:1]*InAxx.dFHz;
signal_locs = find(ismember(freqs,signal_freqs)) ;
noise_locs = reshape(repmat(signal_locs,2,1)+[-1;+1],1,[]);

cmplx_full = [InAxx.Cos;InAxx.Cos(end-1:-1:2,:,:)] + 1i*[InAxx.Sin;-InAxx.Sin(end-1:-1:2,:,:)];
cmplx_signal =  cmplx_full(signal_locs,:,:);
cmplx_noise =  cmplx_full(noise_locs,:,:);

C_s =reshape(permute(cmplx_signal,[2,1,3]),size(cmplx_signal,2),[])*reshape(permute(cmplx_signal,[2,1,3]),size(cmplx_signal,2),[])';
C_n =reshape(permute(cmplx_noise,[2,1,3]),size(cmplx_noise,2),[])*reshape(permute(cmplx_noise,[2,1,3]),size(cmplx_noise,2),[])';

if sum(abs(imag(C_n(:))))/sum(abs(real(C_n(:))))>10^-10
    error('SSD: Covariance matrix is complex')
end
if sum(abs(imag(C_s(:))))/sum(abs(real(C_s(:))))>10^-10
    error('SSD: Covariance matrix is complex')
end
C_n = real(C_n);
C_s = real(C_s);



if opt.do_whitening % and deflate matrix dimensionality
    [V, D] = eig(C_s);
    [ev, desc_idxs] = sort(diag(D), 'descend');
    V = V(:,desc_idxs);
    
    dims = sum((ev/sum(ev))>opt.rank_ratio) ;
    P = V(:,1:dims)*diag(1./sqrt(ev(1:dims)));
else
    dims = size(C_s,1);
    P = eye(dims) ;
end

C_s_w = P' * C_s * P;
C_n_w = P' * C_n * P;


if max(max(abs(C_s_w-C_s_w')))>10^-10
    error('SSD: Whitened signal covariance matrix is not symmetric')
end

if max(max(abs(C_n_w-C_n_w')))>10^-10
    error('SSD: Whitened noise covariance matrix is not symmetric')
end

C_s_w = (C_s_w+C_s_w')/2 ;
C_n_w = (C_n_w+C_n_w')/2 ;

[W,D] =eig(C_s_w,C_n_w);
[D,sorted_idx] = sort(diag(D),'descend') ;
W = W(:,sorted_idx);
W = P * W;


A = C_s * W * pinv(W'*C_s*W);

if max(imag(A(:)))>0
    error('should not be complex!!!')
end

% project to ssd domain
OutAxx = InAxx ;

temp = W'*reshape(permute(InAxx.Cos,[2,1,3]),size(InAxx.Cos,2),[]);
OutAxx.Cos = permute(reshape(temp,dims,size(InAxx.Cos,1),size(InAxx.Cos,3)),[2,1,3]);

temp = W'*reshape(permute(InAxx.Sin,[2,1,3]),size(InAxx.Sin,2),[]);
OutAxx.Sin = permute(reshape(temp,dims,size(InAxx.Sin,1),size(InAxx.Sin,3)),[2,1,3]);

OutAxx.Amp = abs(OutAxx.Cos +1i *OutAxx.Sin);

temp = W'*reshape(permute(InAxx.Wave,[2,1,3]),size(InAxx.Wave,2),[]);
OutAxx.Wave = permute(reshape(temp,dims,size(InAxx.Wave,1),size(InAxx.Wave,3)),[2,1,3]);

    
    
