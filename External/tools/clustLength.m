function lenOut = clustLength(inData)
    spanLocs = bwlabel(inData);   %identify contiguous ones
    spanLength = regionprops(spanLocs, 'area');  %length of each span
    lenOut = [ spanLength.Area];
    if isempty(lenOut)
        lenOut = 0;
    else       
    end
end

