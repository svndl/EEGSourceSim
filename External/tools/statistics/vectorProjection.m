function [project_amp,project_real, project_imag] = vectorProjection( real_in, imag_in, testFig  )
    if nargin < 3
        testFig = false;
    else
    end
    
    if any(size(real_in) ~= size(imag_in))
        msg = '\n real and imaginary data must be same size! \n';
        error(msg);
    else
    end
    
    inSize = size(real_in);
    numConds = prod( inSize(2:end) );
    numSubs = size(real_in,1);    
    project_amp = nan(size(real_in));
    project_real = nan(size(real_in));
    project_imag = nan(size(real_in));
    
    for c = 1:numConds
        xyData = cat(2,real_in(:,c),imag_in(:,c));
        not_nan = ~any(isnan(xyData),2);
        xyData = xyData(not_nan,:);
        xyMean = repmat(mean(xyData,1),size(xyData,1),1);
        lenC = dot(xyData, xyMean, 2) ./ sum(xyMean .* xyMean, 2);
        xyOut = bsxfun(@times, lenC, xyMean);
        project_amp(not_nan,c) = sqrt(xyOut(:,1).^2+xyOut(:,2).^2).*sign(lenC);
        project_real(not_nan,c) = xyOut(:,1);
        project_imag(not_nan,c) = xyOut(:,2);
        if testFig
            figure;
            outMean = mean(xyOut,1);
            plot([0,xyMean(1,1)],[0,xyMean(1,2)],'-kx');
            hold on
            plot(outMean(1,1),outMean(1,2),'go');
            arrayfun(@(x,y) plot(xyData(x,1),xyData(x,2),'ro'), 1:numSubs,'uni',false);
            arrayfun(@(x,y) plot(xyOut(x,1),xyOut(x,2),'bo'), 1:numSubs,'uni',false);
            arrayfun(@(x,y) plot([xyData(x,1),xyOut(x,1)],[xyData(x,2),xyOut(x,2)],'k-'), 1:numSubs,'uni',false);
            axis equal
            xlim([-4,4])
            ylim([-4,4])
        else
        end
    end
end

