function [sampMu,sampCovMat,smaller_eigenvec,smaller_eigenval,larger_eigenvec,larger_eigenval,phi] = eigFourierCoefs(xyData)
% [sampMu,sampCovMat,smaller_eigenvec,smaller_eigenval,larger_eigenvec,larger_eigenval,phi] = eigFourierCoefs(xyData)
% 
% Perform the eigenvalue decomposition on the 2D data xyData and return
% some useful calculations (e.g. for fitting an ellipse to the
% distribution or for conducting a 2D statistical test). 
%
% This function assumes that the 2D data in xyData are Fourier coefficients
% organized with the real coefficients in the 1st column and the imaginary
% coefficients in the 2nd column. Rows = samples.

dims = size(xyData);
N = dims(1);
if dims(2) ~= 2
    error('input data sampMust be a matrix of 2D row samples');
end
if N < 2
    error('input data sampMust contain at least 2 samples');
end

srIx = 1;
siIx = 2;

sampMu = mean(xyData);
sampCovMat = cov([xyData(:,srIx),xyData(:,siIx)]);
[eigenvec, eigenval] = eig(sampCovMat);

% sort the eigenvectors by their eigenvalues
[orderedVals,eigAscendIx] = sort(diag(eigenval));
smaller_eigenvec = eigenvec(:,eigAscendIx(1));
larger_eigenvec = eigenvec(:,eigAscendIx(2));
smaller_eigenval = orderedVals(1);
larger_eigenval = orderedVals(2);

phi = atan2(larger_eigenvec(2), larger_eigenvec(1));
% This angle is between -pi and pi, shift to 0 and 2pi:
if(phi < 0)
    phi = phi + 2*pi;
end