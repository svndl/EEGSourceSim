function [OutAxxs,W,A,D] = CSP(InAxxs,varargin)

    % This function  finds a spatial filter based on the the CSP of
        % InAxxs{1} and InAxxs{1}. 
    %   <options>:
        % freq_range: Frequnecies of signal considered. 
        %             Default: all, except 0 Hz
        % do_whitening: Apply whitening and dimensionality deflation prior to
        %               estimation of spatial filter. Default: true.
        % rank_ratio:   Threshold for dimensionality reduction based on 
        %               eigenvalue spectrum. Default: 10^-4
    % OUTPUT:
        % OutAxx: Data in component space in Axx format
        % W: Spatial filter
        % A: Activation pattern
        % D: Eigenvalues
    % 
    % Written by Sebastian Bosse, 10.8.2018

    
opt	= ParseArgs(varargin,...
    'do_whitening', true, ...
    'rank_ratio', 10^-4, ... 
    'freq_range', InAxxs{1}.dFHz*[1:(InAxxs{1}.nFr-1)] ...
    );

if length(InAxxs)~=2
    error('CSP: CSP not implemented for more classes than 2')
end
if InAxxs{1}.nFr~=InAxxs{2}.nFr
    error('CSP: Number of frequency bins not identical in the two given classes')
end
if InAxxs{1}.dFHz~=InAxxs{2}.dFHz
    error('CSP: Frequency resolution not identical in the two given classes')
end
if InAxxs{1}.nT ~= InAxxs{2}.nT
    error('number of time points differs between classes') ;
end
if InAxxs{1}.nCh ~= InAxxs{2}.nCh
    error('number of channels differs between classes') ;
end
if InAxxs{1}.nTrl ~= InAxxs{2}.nTrl
    warning('number of trials differs between classes') ;
end

cmplx_signal = cell(length(InAxxs),1) ;
    
freq_idxs = 1+opt.freq_range/InAxxs{1}.dFHz ; % shift as 0Hz has idx 1

for class_idx = 1:length(InAxxs)
    if freq_idxs(1)==1 % 0Hz is considered
        cmplx_signal{class_idx} = cat(1, InAxxs{class_idx}.Cos(freq_idxs(2:end),:,:),InAxxs{class_idx}.Cos(freq_idxs,:,:)) ... % even real part 
                    + cat(1, -InAxxs{class_idx}.Sin(freq_idxs(2:end),:,:),InAxxs{class_idx}.Sin(freq_idxs,:,:)); % odd imag part

    else
        cmplx_signal{class_idx} = cat(1, InAxxs{class_idx}.Cos(freq_idxs,:,:),InAxxs{class_idx}.Cos(freq_idxs,:,:)) ... % even real part 
                    + cat(1, -InAxxs{class_idx}.Sin(freq_idxs,:,:),InAxxs{class_idx}.Sin(freq_idxs,:,:)); % odd imag part    
    end
    
end    

C = zeros(size(cmplx_signal{1},2), size(cmplx_signal{1},2), length(InAxxs)) ;
% calculate covariance matrices in fourier domain
for class_idx = 1:length(InAxxs)
    C(:,:,class_idx) = conj(reshape(permute(cmplx_signal{class_idx},[2,1,3]),size(cmplx_signal{class_idx},2),[]))*(reshape(permute(cmplx_signal{class_idx},[2,1,3]),size(cmplx_signal{class_idx},2),[])');
    % normalize by number of observations
    C(:,:,class_idx) = C(:,:,class_idx)/  size(cmplx_signal{class_idx},3);
end

if opt.do_whitening % and deflate matrix dimensionality
    [V, D] = eig(sum(C,3));
    [ev, desc_idxs] = sort(diag(D), 'descend');
    V = V(:,desc_idxs);
    
    dims = sum((ev/sum(ev))>opt.rank_ratio) ;
    P = V(:,1:dims)*diag(1./sqrt(ev(1:dims)));
else
    dims = size(C,1) ;
    P = eye(dims);

end

C_w = zeros(size(P,2),size(P,2),length(InAxxs)) ;

for class_idx = 1:length(InAxxs)
    C_w(:,:,class_idx) = P' * C(:,:,class_idx) * P ;
    
    if max(max(abs(C(:,:,class_idx)-C(:,:,class_idx)')))>10^-10
        error('CSP: Whitened  covariance matrix is not symmetric')
    end
    C_w(:,:,class_idx) = (C_w(:,:,class_idx)+C_w(:,:,class_idx)')/2 ;
end






if opt.do_whitening % calculation in 'white space'
    [W,D] = eig(C_w(:,:,1)-C_w(:,:,2)) ;
else % calculation in channel space
    [W,D] = eig(C_w(:,:,1)-C_w(:,:,2),sum(C_w,3)) ;
end
[D,sorted_idx] = sort(diag(D),'descend') ;
W = W(:,sorted_idx);
W = P * W;

A = sum(C,3) * W *pinv(W'*sum(C,3)*W);

if max(imag(A(:)))>0
    error('should not be complex!!!')
end


OutAxxs = InAxxs ;
for class_idx = 1:length(InAxxs)
    % project to csp domain
    temp = W'*reshape(permute(InAxxs{class_idx}.Cos,[2,1,3]),size(InAxxs{class_idx}.Cos,2),[]);
    OutAxxs{class_idx}.Cos = permute(reshape(temp,dims,size(InAxxs{class_idx}.Cos,1),size(InAxxs{class_idx}.Cos,3)),[2,1,3]);
    temp = W'*reshape(permute(InAxxs{class_idx}.Sin,[2,1,3]),size(InAxxs{class_idx}.Sin,2),[]);
    OutAxxs{class_idx}.Sin = permute(reshape(temp,dims,size(InAxxs{class_idx}.Sin,1),size(InAxxs{class_idx}.Sin,3)),[2,1,3]);
    OutAxxs{class_idx}.Amp = abs(OutAxxs{class_idx}.Cos +1i *OutAxxs{class_idx}.Sin);
    temp = W'*reshape(permute(InAxxs{class_idx}.Wave,[2,1,3]),size(InAxxs{class_idx}.Wave,2),[]);
    OutAxxs{class_idx}.Wave = permute(reshape(temp,dims,size(InAxxs{class_idx}.Wave,1),size(InAxxs{class_idx}.Wave,3)),[2,1,3]);
    
end

