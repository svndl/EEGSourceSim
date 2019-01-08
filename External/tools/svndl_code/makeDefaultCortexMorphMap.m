function [mapMtx totalTrans] = makeDefaultCortexMorphMap(fromSub,toSub)
%function [mapMtx] = makeDefaultCortexMorphMap(fromSub,toSub)
%
% This function makes a morphing matrix that maps values from 1 subject to another
% Freesurfer based MNE morphing matrices must have been computed already
% with command:
% mne_make_morph_maps --to skeri0055_fs4 --from skeri0001_fs4
%
%example usage:
%[mapMtx] = makeDefaultCortexMorphMap('skeri0001','skeri0055');

% $Log: makeDefaultCortexMorphMap.m,v $
% Revision 1.4  2009/06/26 21:12:30  ales
% mrSimScript does arbitrary waveforms
% skeriDefaultSimParameters has updated help
% makeDefaultCortexMorphMap returns an identity matrix for mapping an individual on to itself
% makeForwardMatrixFromMne has some added help
%
% Revision 1.3  2009/05/29 16:58:14  ales
% Changed something. I don't remember what now. sorry.
%
% Revision 1.2  2009/05/21 17:04:23  ales
% Added auto log message comments
%

% tell user what you're doing
fprintf('making morphing matrix %s to template %s \n',fromSub,toSub)

% Read stuff in

toSub = [toSub '_fs4'];
fromSub = [fromSub '_fs4'];

freesurfDir = getpref('freesurfer','SUBJECTS_DIR');


toSubSrcSpaceFile = fullfile(freesurfDir,toSub,'bem',[ toSub '-ico-5p-src.fif']);
fromSubSrcSpaceFile = fullfile(freesurfDir,fromSub,'bem',[ fromSub '-ico-5p-src.fif']);

if strcmp(fromSub,toSub)   
    toSubSrc = mne_read_source_spaces(toSubSrcSpaceFile);   
    mapMtx = speye(sum([toSubSrc.nuse])); %<- tricky use of [] to make a matrix from structure
    totalTrans = [];
    return;
end

toSubSrc = mne_read_source_spaces(toSubSrcSpaceFile);
fromSubSrc = mne_read_source_spaces(fromSubSrcSpaceFile);


[hi2hi{1},hi2hi{2}] = mne_read_morph_map(fromSub,toSub,freesurfDir);


%% Make lowrez transfer matrix

for iHemi = 1:2,

fromDecIdx = double(fromSubSrc(iHemi).vertno); %Decimation index for fromSub
toDecIdx   = double(toSubSrc(iHemi).vertno); %Decimation index for toSub

nHiFrom = double(fromSubSrc(iHemi).np); %Num  vertices hires in from subject
nHiTo = double(toSubSrc(iHemi).np); %Num  vertices hires in "to" subject

nLoTo = double(toSubSrc(iHemi).nuse);
nLoFrom = double(fromSubSrc(iHemi).nuse);

idxFrom = nearpoints(fromSubSrc(iHemi).rr',fromSubSrc(iHemi).rr(fromDecIdx,:)');
idxTo = nearpoints(toSubSrc(iHemi).rr',toSubSrc(iHemi).rr(toDecIdx,:)');

%Sparse mapping matrix from low res to hi res surface using nearest
%neighbour interpolation.
fromLo2Hi = sparse(1:nHiFrom,idxFrom,ones(size(idxFrom)),nHiFrom,nLoFrom);
toLo2Hi = sparse(1:nHiTo,idxTo,ones(size(idxTo)),nHiTo,nLoTo);


%Indexing matrix from Hi res "to" subject to low res "to" subject
toHi2Lo = sparse(1:nLoTo,toDecIdx,ones(nLoTo,1),nLoTo,nHiTo);



%Fix if missing some end columns because of a zero mapping
[i j s] = find(hi2hi{iHemi});
hi2hi{iHemi} = sparse(i,j,s,nHiTo,nHiFrom);
 
totalTrans{iHemi} = toHi2Lo*(hi2hi{iHemi}*fromLo2Hi);

end


%Make the hemi parts into 1 big matrix suitable for our default cortex.
mapMtx = [totalTrans{1}  sparse(size(totalTrans{1},1),size(totalTrans{2},2)); ...
    sparse(size(totalTrans{2},1),size(totalTrans{1},2)) totalTrans{2}];


