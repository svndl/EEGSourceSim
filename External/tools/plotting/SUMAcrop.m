function SUMAcrop(imNames,preCrop,blackBg)
    if nargin < 3
        blackBg = false;
    else
    end
    if nargin < 2
        preCrop = 0;
    else
    end
    if numel(preCrop)==1
        preCrop = repmat(preCrop,1,4);
    elseif numel(preCrop)==2
        preCrop(3) = preCrop(1);
        preCrop(4) = preCrop(2);
    else
    end
    for a=1:length(imNames)
        inPut = imread(imNames{a});
        if ndims(inPut) == 2
            inPut = repmat(inPut,[1,1,3]);
        else
        end
        inPut = inPut(preCrop(2):(end-preCrop(4)),preCrop(1):(end-preCrop(3)),1:3);
        tempIm = sum(inPut,3);
        imSize = size(tempIm);
        tempX = sum(tempIm,1)';
        tempY = sum(tempIm,2);
        if blackBg
            xW = 0*3*imSize(1);
            yW = 0*3*imSize(2);
            xIdx = tempX>xW;
            yIdx = tempY>yW;
        else
            xW = 255*3*imSize(1);
            yW = 255*3*imSize(2);
            xIdx = tempX<xW;
            yIdx = tempY<yW;
        end
        xNew = length(find(xIdx));
        yNew = length(find(yIdx));
        allIdx=repmat(xIdx,1,imSize(1))'.*repmat(yIdx,1,imSize(2));
        allIdx = repmat(allIdx,[1,1,3]);
        cropIm = reshape(inPut(logical(allIdx)),yNew,xNew,3);
        if blackBg
            tLayer = double(sum(cropIm,3)>(0*3));
        else
            tLayer = double(sum(cropIm,3)<(255*3));
        end
        newName = [imNames{a}(1:end-4),'_crop.png'];
        imwrite(cropIm,newName,'png','Alpha',tLayer)
    end
end

