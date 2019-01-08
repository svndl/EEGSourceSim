classdef axx
    % This class defines standard data type that is exported by xDiva. Axx
    % contains spectral and time domain representation of the multichannel
    % steady-state brain response signal.
    
    %----------------------------------------------------------------------
    % Author: Peter Kohler,
    % Last modification: Elham Barzegaran, 03/26/2018
    %======================================================================

    properties
        cndNmb
        nTrl % how many trials
    end
    
    properties (Dependent)
        nT   % how many time points: for a period of a single stimulus cycle
        nCh  % how many channels
    end
    
    properties
        dTms % time resolution in miliseconds
        dFHz % frequency resolution
        nFr  % how many frequencies
        i1F1 % fundamental frequency 1
        i1F2 % fundamental frequency 1
        DataUnitStr %like 'microVolts' 
        Amp % Amplitude spectrum of EEG: a nFr x nCh x nTrl matrix
        Cos % Cosine part of spectrum (real): a nFr x nCh x nTrl matrix
        Sin % Sine part of spectrum (imaginary): a nFr x nCh x nTrl matrix
        SpecPValue
        SpecStdErr
        Cov
        Wave % Averaged EEG in time domain: a nT x nCh x nTrl matrix: 
        %%% Note that nT might be different from the time points used for
        %%% spectrum estimation, since nT is the length of single stimulus 
        %%% cycle, while time window for spectrum is selected so that dFHz 
        %%% is about 0.5 Hz and time window lengths is an integer of stimulus cycle
    end
    
    methods
        function obj = axx(axxStrct,AVR)
            if nargin<2
                AVR=1;% average over trials?
            end
             % class constructor
            if nargin > 0
                % unchanged values, use first strct (will also work if oldStrct is just one 
                obj.dTms = axxStrct(1).dTms;
                obj.dFHz = axxStrct(1).dFHz; 
                obj.i1F1 = axxStrct(1).i1F1;
                obj.i1F2 = axxStrct(1).i1F2;
                obj.DataUnitStr = axxStrct(1).DataUnitStr;
                obj.cndNmb = axxStrct(1).cndNmb;
                obj.nTrl = axxStrct(1).nTrl;
                obj.nFr = axxStrct(1).nFr;

                % average together the following values Wave, Amp, Cos, Sin, Cov, SpecPValue, SpecStdErr
                if AVR ==1,
                    obj.Wave = mean(cat(3,axxStrct.Wave),3);
                    ampVal = mean(cat(3,axxStrct.Amp),3);
                    obj.Amp = ampVal(1:obj.nFr,:);
                    cosVal = mean(cat(3,axxStrct.Cos),3);
                    obj.Cos = cosVal(1:obj.nFr,:);
                    sinVal = mean(cat(3,axxStrct.Sin),3);
                    obj.Sin = sinVal(1:obj.nFr,:);
                    if isfield(axxStrct, 'SpecPValue')
                        specpVal = mean(cat(3,axxStrct.SpecPValue),3);
                        obj.SpecPValue = specpVal(1:obj.nFr,:);
                    else
                        obj.SpecPValue = [];
                    end
                    if isfield(axxStrct, 'SpecStdErr')
                        specstdVal = mean(cat(3,axxStrct.SpecStdErr),3);
                        obj.SpecStdErr = specstdVal(1:obj.nFr,:);
                    else
                        obj.SpecStdErr = [];
                    end
                    if isfield(axxStrct, 'Cov')
                        obj.Cov = mean(cat(3,axxStrct.Cov),3);
                    else
                        obj.SpecStdErr = [];
                    end
                elseif AVR==0,
                    obj.Wave = axxStrct.Wave;
                    obj.Amp = axxStrct.Amp(1:obj.nFr,:,:);
                    obj.Cos = axxStrct.Cos(1:obj.nFr,:,:);
                    obj.Sin = axxStrct.Sin(1:obj.nFr,:,:);
                    
                    if isfield(axxStrct, 'SpecPValue')
                       obj.SpecPValue = xxStrct.SpecPValue(1:obj.nFr,:,:);
                    else
                        obj.SpecPValue = [];
                    end
                    
                    if isfield(axxStrct, 'SpecStdErr')
                        obj.SpecStdErr = axxStrct.SpecStdErr(1:obj.nFr,:,:);
                    else
                        obj.SpecStdErr = [];
                    end
                    
                    if isfield(axxStrct, 'Cov')
                        obj.Cov = axxStrct.Cov;
                    else
                        obj.SpecStdErr = [];
                    end
                end
            end
        end 
        function value = get.nT(obj)
            value = size(obj.Wave,1);
        end
        function value = get.nCh(obj)
            value = size(obj.Wave,2);
        end
        
        function identify(thisAxx)
            if isempty(thisAxx.Wave)
                disp('I am an axx file. I am empty.');
            else
                disp('I am an axx file. I contain data.');
            end
        end
        
        function writetofile(thisAxx,FilePath)
            Fields = fieldnames(thisAxx);
            for f = 1:numel(Fields)
                eval(['Sstruct.' Fields{f} '=' 'thisAxx.' Fields{f} ';']);
            end
            save(FilePath,'-struct','Sstruct');
        end
        function s = saveobj(obj)
            s.cndNmb = obj.cndNmb ;
            s.nTrl = obj.nTrl ;
            s.nT = obj.nT ;
            s.nCh = obj.nCh ;
            s.dTms = obj.dTms ;
            s.dFHz = obj.dFHz ;
            s.nFr   = obj.nFr ;
            s.i1F1  = obj.i1F1 ;
            s.i1F2  = obj.i1F2 ;
            s.DataUnitStr = obj.DataUnitStr;
            s.Amp  = obj.Amp ;
            s.Cos  = obj.Cos ;
            s.Sin  = obj.Sin ;
            s.SpecPValue = obj.SpecPValue ;
            s.SpecStdErr = obj.SpecStdErr ;
            s.Cov = obj.Cov ;
            s.Wave  = obj.Wave;
        end
        
      function outAxx = SelectTrials(obj,TIdx)
        % Select the trials indicated uin TIdx vector.
        if max(TIdx)>obj.nTrl || min(TIdx)<=0
            error('Wrong trial indexes');
        end

        obj.nTrl = numel(TIdx);
        obj.Amp = obj.Amp(:,:,TIdx);
        obj.Cos = obj.Cos(:,:,TIdx);
        obj.Sin = obj.Sin(:,:,TIdx);
        obj.Wave = obj.Wave(:,:,TIdx);
        outAxx = obj;
      end
      
      function outAxx = MergeTrials(obj1,obj2)
        % Select the trials indicated uin TIdx vector.
        if obj1.nT~=obj2.nT || obj1.nCh~=obj2.nCh || obj1.dTms~=obj2.dTms || obj1.dFHz~=obj2.dFHz
            error('Axx classes do not match: Time and frequency features should be the same...');
        end

        outAxx = obj1;
        outAxx.nTrl = obj1.nTrl+obj2.nTrl;
        outAxx.Amp = cat(3,obj1.Amp,obj2.Amp);
        outAxx.Cos = cat(3,obj1.Cos,obj2.Cos);
        outAxx.Sin = cat(3,obj1.Cos,obj2.Sin);
        outAxx.Wave = cat(3,obj1.Wave,obj2.Wave);

    end
    end
    
      
  methods(Static)
      function obj = loadobj(s)
            if isstruct(s)
                newObj = mrC.axx() ;
                newObj.cndNmb = s.cndNmb ;
                newObj.nTrl = s.nTrl ;
%                 newObj.set.nT(s.nT) ;
%                 newObj.set.nCh(s.nCh) ;
                newObj.dTms = s.dTms ;
                newObj.dFHz = s.dFHz ;
                newObj.nFr   = s.nFr ;
                newObj.i1F1  = s.i1F1 ;
                newObj.i1F2  = s.i1F2 ;
                newObj.DataUnitStr = s.DataUnitStr;
                newObj.Amp  = s.Amp ;
                newObj.Cos  = s.Cos ;
                newObj.Sin  = s.Sin ;
                newObj.SpecPValue = s.SpecPValue ;
                newObj.SpecStdErr = s.SpecStdErr ;
                newObj.Cov = s.Cov ;
                newObj.Wave  = s.Wave;
                obj = newObj ;
            else
                obj = s ;
            end
      end
  end
    
end