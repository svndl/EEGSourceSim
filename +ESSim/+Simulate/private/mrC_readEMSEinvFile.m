function Inv = mrC_readEMSEinvFile(filename)
% Inv = mrC_readEMSEinvFile(filename)
% returns nChannels x nVertices matrix
% 
% based on emseReadInverse, modified to fread nRows x nCols bytes
% beginning at (assumed) end of header rather than by fseeking back from EOF.
% Thus, this implementation should read .inv files with or without the xml-ish footer.

% $Log: mrC_readEMSEinvFile.m,v $
% Revision 1.3  2009/11/18 01:55:03  nicholas
% merged into main branch
%
% Revision 1.2.2.2  2009/11/12 20:45:18  nicholas
% changed 128 channel check from error to warning
%
% Revision 1.2.2.1  2009/07/21 18:17:31  nicholas
% *** empty log message ***
%
% Revision 1.2  2008/11/07 00:03:00  ales
% Further squashed the emse style inverse binary/text reading bug
%
	fid = fopen(filename,'rb','ieee-le');
	if fid == -1
		error('Error opening %s',filename)
	end
	% Get the magic number
	magicNum = upper(fscanf(fid,'%c',8));
	if strcmp(magicNum,'454D5345') % Magic number is OK, we're reading a real inverse file
		% Next read in the major and minor revs, and other header fields.
		% Based on the file format description in Appendix A of EMSE's help file,
		% we expect exactly ten elements in Header, with the dimensions of the inverse
		% matrix in the 9th and 10th position.  Here's what can go wrong:
		% 1) SSI might revise inverse file header field structure without warning us.
		% 2) There will be two extra fields if "cortical thinning was used", whatever that means.
		% For now, this implementation simply checks whether fscanf returns less than expected number of
		% Header elements, otherwise throwing an error.  It falls to you, dear reader,
		% to implement handling of the remaining possibilities listed above should they ever occur.
		[Header,nHeader] = fscanf(fid,'%d',10);
		% fscanf is not robust to bytes following the last header field that have degenerate ASCII values;
		% so we use fgetl, which seems to behave correctly;
                % These lines did not completely fix the bug. Changed to an explicit fseek,
		% see line below if block: fseek(fid,1,0)
		% nHeader = nHeader + 1;
		% Header(nHeader) = str2num(fgetl(fid));

		if nHeader ~= 10
			error('Expecting 10 header elements, found %d in %s',nHeader,filename)
        	end	
		
		%This line explicity sets the file read position to what we think is the begining of good data.
	        fseek(fid,1,0);
        
		nRows = Header(9);
		nCols = Header(10);
		if nCols ~= 128
% 			error('Expecting 128 columns in inverse, found %d in %s',nCols,filename);
			warning('Expecting 128 columns in inverse, found %d in %s',nCols,filename);
		end
		[Inv,nInv] = fread(fid,nCols*nRows,'float64',0,'ieee-le');
		fclose(fid);
		if nInv ~= nCols*nRows
			error( 'Size of inverse (%d) does not match dimensions in file header (%d*%d=%d)',nInv,nRows,nCols,nRows*nCols)
		end
		Inv = reshape( Inv, nCols, nRows );		% nChannels x nVertices
	else
		error('Magic# in %s = %s, expecting 454D5345',filename,magicNum)
	end
end

