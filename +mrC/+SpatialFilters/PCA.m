function [OutAxx,W,A,D] = PCA(InAxx,varargin)

% This function calculates finds a spatial filter based on the the PCA of
% InAxx. 

% INPUT:
    % InAxx: EEG data in Axx format
%   <options>:
    % freq_range: Frequnecies of signal considered. 
    %             Default: all, except 0 Hz
    
% OUTPUT:
    % OutAxx: Data in component space in Axx format
    % W: Spatial filter
    % A: Activation pattern
    
% Written by Sebastian Bosse, 3.8.2018

opt	= ParseArgs(varargin,...
    'freq_range', InAxx.dFHz*[1:(InAxx.nFr-1)], ...
    'model_type','complex' ...
    );


freqs = [0:InAxx.nFr-2,InAxx.nFr-1:-1:1]*InAxx.dFHz;
locs = ismember(freqs,opt.freq_range) ;

% build full cmplx spectrum 
cmplx_signal = [InAxx.Cos;InAxx.Cos(end-1:-1:2,:,:)] + 1i*[InAxx.Sin;-InAxx.Sin(end-1:-1:2,:,:)];
% reduce cmplx spectrum according to desired freq_range
cmplx_signal = cmplx_signal(locs,:,:) ;

if strcmpi(opt.model_type,'cartesian')
    cmplx_signal = [real(cmplx_signal);imag(cmplx_signal)] ;
end

C =reshape(permute(cmplx_signal,[2,1,3]),size(cmplx_signal,2),[])*(reshape(permute(cmplx_signal,[2,1,3]),size(cmplx_signal,2),[]))';


if sum(abs(imag(C(:))))/sum(abs(real(C(:))))>10^-10
    error('PCA: Covariance matrix is complex')
end
C = real(C);

[W,D] = eig(C);
[D,sorted_idx] = sort(diag(D),'descend') ;
W = W(:,sorted_idx);
A = C * W * pinv(W'*C*W);
if max(imag(A(:)))>0
    error('should not be complex!!!')
end

% project to pca domain
OutAxx = InAxx ;
temp = W'*reshape(permute(InAxx.Cos,[2,1,3]),size(InAxx.Cos,2),[]);
OutAxx.Cos = permute(reshape(temp,size(InAxx.Cos,2),size(InAxx.Cos,1),size(InAxx.Cos,3)),[2,1,3]);

temp = W'*reshape(permute(InAxx.Sin,[2,1,3]),size(InAxx.Sin,2),[]);
OutAxx.Sin = permute(reshape(temp,size(InAxx.Sin,2),size(InAxx.Sin,1),size(InAxx.Sin,3)),[2,1,3]);

OutAxx.Amp = abs(OutAxx.Cos +1i *OutAxx.Sin);

temp = W'*reshape(permute(InAxx.Wave,[2,1,3]),size(InAxx.Wave,2),[]);
OutAxx.Wave = permute(reshape(temp,size(InAxx.Wave,2),size(InAxx.Wave,1),size(InAxx.Wave,3)),[2,1,3]);


    
    
