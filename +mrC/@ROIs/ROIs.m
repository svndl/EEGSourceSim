classdef ROIs
    % This class defines a data type for storing and retrieving ROIs in mrC. 
    
    %----------------------------------------------------------------------
    % Author: Elham Barzegaran, 06/14/2018
    %======================================================================
       
    properties
        ROIList % a cell array of roi structure that keeps all the information about ROIs
        % roi structure includes roi.Type, roi.Name, roi.Hemi, roi.meshIndices, roi.Date, roi.Comment 
        subID
    end
    
    properties (Dependent)
        Atlases   % How Many atlases are there
        ROINum  % how ROIs in general
    end
    
    methods
        %-----------------------Initializing-------------------------------
        function obj = ROIs(subID,anatDir)
            if nargin == 0
                % just  loading
                return
            end
            % gets subject IDs and Anatomy Path and generate ROIs class
            % with all ROIs available for that subject
            %%%%%%%%%%%%%%%%%%%%%%%%this function should be later updated
            %%%%%%%%%%%%%%%%%%%%%%%% to replace ROIFromSuma%%%%%%%%%%%%%%%%
            ROIList = [];  
            Atlases = [];
            if ~exist('subID','var') || isempty(subID)
                subID = 'NA';
            end
            obj.subID = subID;
            if exist('subID','var')
                if ~exist('anatDir','var') || isempty(anatDir) || ~exist(anatDir,'dir')
                    anatDir = getpref('mrCurrent','AnatomyFolder');
                    if ~exist(anatDir,'dir')
                        warning('No anatomy folder is found');
                        obj.ROIList = ROIList;
                        return
                    end 
                end
                if ~exist(fullfile(anatDir,subID),'dir')
                    if ~strcmp(subID,'NA') % if it is not an empty class
                        warning('Subject not found in anatomy folder');
                    end
                    obj.ROIList = ROIList;
                    return
                end

                roiDir = fullfile(anatDir,subID,'Standard','meshes'); % where matlab files are stored
                roiFolders = subfolders(roiDir);
                
                roifoldind = cellfun(@(x) ~isempty(x),strfind(roiFolders,'benson')) |...
                    cellfun(@(x) ~isempty(x),strfind(lower(roiFolders),'wang')) |...
                    cellfun(@(x) ~isempty(x),strfind(lower(roiFolders),'kgs')) | ...
                    cellfun(@(x) ~isempty(x),strfind(lower(roiFolders),'glass'));
                
                roiFolders = roiFolders(roifoldind);
                
                for f = 1:numel(roiFolders)
                    if strfind(lower(roiFolders{f}),'_rois')
                       nroiList =  subfiles(fullfile(fullfile(roiDir,roiFolders{f}),'/*.mat'),1);
                       
                       if strfind(lower(roiFolders{f}),'benson'),% Only loads the 0-80 eccentricity files for benson
                           ind = cellfun(@(x) ~isempty(x),strfind(nroiList,'0-80'));
                           nroiList = nroiList(ind);
                           if isempty(nroiList)
                               warning('No Benson ROI with [0 80] eccentricity was found. To use Benson ROIs in this class generate them using mrC.RoiFromSuma');
                               mrC.RoiFromSuma(subID,'mode','benson','plotting',false,'ecc_range',[0 80]);
                               nroiList =  subfiles(fullfile(fullfile(roiDir,roiFolders{f}),'/*.mat'),1);
                               ind = cellfun(@(x) ~isempty(x),strfind(nroiList,'0-80'));
                               nroiList = nroiList(ind);
                           end
                       end
                       currRoi = [];
                       if ~isempty(nroiList)
                           currRoi = cellfun(@(x) load(x),nroiList);
                           FNames = cellfun(@(x) x.name, struct2cell(currRoi),'UniformOutput', false);
                           Idx = strfind(FNames,'-');
                           iseccfield = cell2mat(cellfun(@(x) isfield(x,'eccData'), struct2cell(currRoi),'UniformOutput', false));
                           eccfield = cell(size(currRoi));
                           eccfield(iseccfield) = cellfun(@(x) x.eccData, struct2cell(currRoi(iseccfield)),'UniformOutput', false);
                           
                           % save wangatlas as wang
                           if strfind(lower(currRoi(1).ROI.type),'wang')
                               Type = 'wang';
                           else
                               Type = lower(currRoi(1).ROI.type);
                           end
                           
                           
                           % create the ROI strcture
                           currRoi = struct('Type', Type, ...%cellfun(@(x) x.ype, struct2cell(currRoi),'UniformOutput', false),...
                               'Name',cellfun(@(x,y) x(1:y(end)-1), FNames,Idx,'UniformOutput', false),...
                               'Hemi',cellfun(@(x,y) x(y(end)+1:end), FNames,Idx,'UniformOutput', false),...
                               'meshIndices',cellfun(@(x) x.meshIndices, struct2cell(currRoi),'UniformOutput', false),...
                               'eccData',eccfield',...
                               'Date',cellfun(@(x) x.date, struct2cell(currRoi),'UniformOutput', false),...
                               'Comment',cellfun(@(x) x.comment, struct2cell(currRoi),'UniformOutput', false)...
                               );

                       end
                       ROIList = [ROIList currRoi];
                    end
                end
            end
            obj.ROIList = ROIList;

        end
        
        %----------------------Dependent variables-------------------------
        function value = get.ROINum(obj)
            value = numel(obj.ROIList);
        end
        
        function value = get.Atlases(obj)
            if ~isempty(obj.ROIList)
                value = unique({obj.ROIList.Type});
            else
                value = [];
            end
            
        end
        
        %---------------------add-remove-retrieve atlases------------------
        function [obj] = getAtlasROIs(obj,atlas,ecc_range)
            %returns the list of ROIs in one or several atlases
            
            if obj.ROINum~=0,
            
                atlas = lower(atlas);
                if ~strcmpi(atlas,'benson') || ~exist('ecc_range','var')
                    % returns ROIs in the atlas
                    List = obj.ROIList;
                    Ind = strcmpi(atlas,{List.Type});
                    obj.ROIList = obj.ROIList(Ind);

                else
                    if sum(ecc_range<0) || sum(ecc_range>80) || (ecc_range(1)>ecc_range(2))
                        error('ecc_range should be defined as [ecc_min ecc_max] and should be within the range of [0 80]');
                    end
                    % get benson ROIs
                    List = obj.ROIList;
                    Ind = strcmpi(atlas,{List.Type});
                    obj.ROIList = obj.ROIList(Ind);


                    % Adjust ROIs with the specified eccentricity
                    List = obj.ROIList;
                    Ind = cellfun(@(x) ~isempty(x),strfind({List.Name},'0-80'));% read the benson atlas with 0-80 degree
                    List = List(Ind);

                    % if cannot find them generate them with RoiFromSuma
                    if isempty(List)
                        warning('Benson ROIs with eccentricity [0 80] is not found');
                        display('Generating Benson ROIs with  eccentricity [0 80]');
                        mrC.RoiFromSuma(obj.subID,'mode','benson','plotting',false,'ecc_range',[0 80]);
                        obj = mrC.ROIs(obj.subID);
                        obj = obj.getAtlasROIs('benson');

                    end

                     % Adjust ROIs with the specified eccentricity
                    List = obj.ROIList;
                    Ind = cellfun(@(x) ~isempty(x),strfind({List.Name},'0-80'));% read the benson atlas with 0-80 degree
                    List = List(Ind);
                    Ind = cellfun(@(x) x>ecc_range(1) & x<ecc_range(2),{List.eccData},'UniformOutput',false);
                    NewMInd = cellfun(@(x,y) x(y),{List.meshIndices},Ind,'UniformOutput',false);
                    Newecc = cellfun(@(x,y) x(y),{List.eccData},Ind,'UniformOutput',false);
                    [List(:).meshIndices] = deal(NewMInd{:});
                    [List(:).eccData] = deal(Newecc{:});

                    % Adjust the Names
                    Names = {List.Name};
                    Names = cellfun(@(x,y) x(1:y),Names,strfind(Names,'_'),'UniformOutput',false);
                    Names = cellfun(@(x) [x num2str(ecc_range(1)) '-' num2str(ecc_range(2))],Names,'UniformOutput',false);
                    [List(:).Name] = deal(Names{:});

                    obj.ROIList = List;
                end
            end
        end
        
        function obj = removeAtlas(obj,atlas)
            %remove atlases from from current objec
            List = obj.ROIList;
            Ind = strcmpi(atlas,{List.Type});
            List(Ind)=[];
            obj.ROIList= List;
        end
        
%         function obj = addAtlas(obj,atlaslist)
%             % adds or replace atlases to the current object
%         end
        
        function obj = mergROIs(obj,obj2)
            
            if obj.ROINum==0 || obj2.ROINum==0
                if obj.ROINum==0,
                    obj = obj2;
                end
                return;
            end
            
            if ~strcmp(obj.subID,obj2.subID)
                error('Merg is not possible, not the same subject');
            end
            list1 = cellfun(@(x,y,z) [x '_' y '_' z],{obj.ROIList.Type},{obj.ROIList.Name},{obj.ROIList.Hemi},'UniformOutput',false);
            list2 = cellfun(@(x,y,z) [x '_' y '_' z],{obj2.ROIList.Type},{obj2.ROIList.Name},{obj2.ROIList.Hemi},'UniformOutput',false);
            
            % merg and remove the repeated ones with newer ones
            rep = cellfun(@(x) strcmpi(x,list2),list1,'UniformOutput', false);%(1) find the repeated ROIs (with the same atlas)
            [ind1,ind2] = find(cat(1,rep{:}));
           
            %(2) select the newer ones
            if ~isempty(ind1)
                date1 = {obj.ROIList.Date};date1 = date1(ind1);
                date2 = {obj2.ROIList.Date};date2 = date2(ind2);
                repl = cellfun(@(x,y) datenum(x)<datenum(y),date1,date2);
                obj.ROIList(ind1(repl))= obj2.ROIList(ind2(repl));
            end
            
            % add the newer ones
            aind = 1:numel(list2);aind(ind2)=[];
            obj.ROIList = [obj.ROIList obj2.ROIList(aind)];
        end
        
        
        %--------------------------Search ROIs-----------------------------
        function [obj, ROIInd, ROIinfo] = searchROIs(obj,ROIname,Atlas,Hemi)
            % Syntax:  [obj, ROIInd, ROIinfo] = obj.searchROIs(ROIname,Atlas,Hemi)
            %
            % ROIname: is a string or an array of strings.
            % atlas: is an optional input indicate atlas name to look in
            % Hemi: is an optional input indicate the hemisphere to search,
            %       can be an string or an array of strings, L(for left),
            %       R (for right) and B (for both)
            %-----------------------------------
            List = obj.ROIList;
            if isempty(List),% if the ROI class is empty
                ROIInd = [];
                ROIinfo = [];
                return;
            end
            
            if ~iscell(ROIname),
                if strcmp(ROIname,'all'),
                    Ind = true(1,obj.ROINum);
                else
                    Ind = cellfun(@(x) ~isempty(x),strfind(lower({List.Name}),lower(ROIname)));
                end
            else
                arrInd = cellfun(@(x) strfind(lower({List.Name}),x),lower(ROIname),'UniformOutput',false);
                Ind = sum(cellfun(@(x) ~isempty(x),cat(1,arrInd{:})));
            end
            
            
            if exist('Atlas','var') && ~isempty(Atlas)
                Ind2 = cellfun(@(x) ~isempty(x),strfind(lower({List.Type}),lower(Atlas)));
                Ind = (Ind>0) & Ind2;
            end
            
            if exist('Hemi','var') && ~isempty(Hemi)
                if ~iscell(Hemi),
                    if strcmpi(Hemi,'B'),
                        Ind2 = cellfun(@(x) ~isempty(x),strfind(lower({List.Hemi}),lower('L'))) | cellfun(@(x) ~isempty(x),strfind(lower({List.Hemi}),lower('R')));
                    else
                        Ind2 = cellfun(@(x) ~isempty(x),strfind(lower({List.Hemi}),lower(Hemi)));
                    end
                    Ind = (Ind>0) & Ind2;
                else
                    arrInd2 = cellfun(@(x) strfind(lower({List.Hemi}),x),lower(Hemi),'UniformOutput',false);
                    Ind2 = sum(cellfun(@(x) ~isempty(x),cat(1,arrInd2{:})));
                    Ind = (Ind>0) & (Ind2>0);
                end
            end
            
            ROIInd = find(Ind)';
            ROIinfo = List(ROIInd);%[{List(Ind).Name}' {List(Ind).Hemi}' {List(Ind).Type}' ];
            ROIIndc = num2cell(ROIInd); [ROIinfo(:).ROIind] = deal(ROIIndc{:});
            obj.ROIList = obj.ROIList(ROIInd);
            
            % Check the order of ROIs
            
        end
        
        
        function [obj] = selectROIs(obj,ROIInd)
            % this function might not be necessary
            obj.ROIList = obj.ROIList(ROIInd);
        end
        
        function [NameList] = getFullNames(obj,mode)
            
            if (obj.ROINum==0), % if the ROI class is empty
                NameList = [];
                return;
            end
            
            if ~exist('mode','var'),
                mode ='atlas';
            end
            if strcmp(mode,'atlas')
                NameList = cellfun(@(x,y,z) [x '_' y '_' z],{obj.ROIList.Type},{obj.ROIList.Name},{obj.ROIList.Hemi},'UniformOutput',false);
            elseif strcmp(mode,'noatlas')
                NameList = cellfun(@(y,z) [y '_' z],{obj.ROIList.Name},{obj.ROIList.Hemi},'UniformOutput',false);
            elseif strcmp(mode,'noatlashemi')
                NameList = cellfun(@(y,z) [y],{obj.ROIList.Name},'UniformOutput',false);
            else
                warning('The defined mode is wrong, use atlas or noatlas');
                NameList=[];
            end
        end
        
        %-----------------------convert to matrix format-------------------
        function [roiChunk, NameList] = ROI2mat(obj,nsource)
            roiChunk = zeros(nsource,obj.ROINum);
            for r = 1:obj.ROINum % Read each Roi and check if they exist
                roiChunk(obj.ROIList(r).meshIndices,r)=1;
            end
            [NameList] = getFullNames(obj);
        end
        
        %-------------------save and load ROIs class object----------------
        function saveROIs(obj,anatDir)
            if ~exist('anatDir','var') || isempty(anatDir) || ~exist(anatDir,'dir')
                anatDir = getpref('mrCurrent','AnatomyFolder');
                if ~exist(anatDir,'dir')
                    warning('No anatomy folder is found');
                    return
                end 
            end
            if ~exist(fullfile(anatDir,obj.subID),'dir')
                warning('Subject not found in anatomy folder');
                return
            end
            roiDir = fullfile(anatDir,obj.subID,'Standard','meshes');
            save(fullfile(roiDir,'ROIsClass'),'obj');            
        end
        
        function obj = loadROIs(obj,subID,anatDir)
            % make an empty object and load ROIs into it?
            if ~exist('anatDir','var') || isempty(anatDir) || ~exist(anatDir,'dir')
                anatDir = getpref('mrCurrent','AnatomyFolder');
                if ~exist(anatDir,'dir')
                    warning('No anatomy folder is found');
                    return
                end 
            end
            if ~exist(fullfile(anatDir,subID),'dir')
                warning('Subject not found in anatomy folder');
                return
            end
            roiDir = fullfile(anatDir,subID,'Standard','meshes');
            if exist(fullfile(roiDir,'ROIsClass.mat'),'file')
                load(fullfile(roiDir,'ROIsClass.mat'),'obj');
            else
                obj = mrC.ROIs(subID,anatDir);
                obj.saveROIs(anatDir);
            end
        end
        function s = saveobj(obj)
            s.ROIList = obj.ROIList ;
            s.subID = obj.subID ;
        end
    end
         
    methods(Static)
      function obj = loadobj(s)
            if isstruct(s)
                newObj = mrC.ROIs() ;
                newObj.ROIList = s.ROIList ;
                newObj.subID = s.subID ;
                obj = newObj ;
            else
                obj = s ;
            end
      end
  end
    
end