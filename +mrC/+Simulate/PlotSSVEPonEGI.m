function PlotSSVEPonEGI(EEGAxx, SignalType, SavePath, SaveName,signalFF,RoiList,subIDs, subjects,Mode)


% INPUT:
    % EEGAXX: is the data that mrC.Simulate.SimulateProject returns
    % SignalType: Determines if the plots are in phase or amplitude
                % [Amplitude]/Phase
                
                
    % subjects: [Individuals]/Average
    % Mode:     [Interactive]/Simple
%% set default valuse

if ~exist('Mode','var') || isempty(Mode)
    Mode = 'Simple';
end

if ~exist('SignalType','var') || isempty(SignalType)
    SignalType = 'Amplitude';
end

if ~exist('signalFF','var')|| isempty(signalFF)
      signalFF = 1;
end
 
if ~exist('subjects','var') || isempty(subjects)
    subjects = 'Individuals';
end
%-------------------Calculate EEG spectrum---------------------------------
    sub1 = find(~cellfun(@isempty,EEGAxx),1);
    freq = 0:EEGAxx{sub1}.dFHz:EEGAxx{sub1}.dFHz*(EEGAxx{sub1}.nFr-1); % frequncy labels, based on fft

    for s = 1:length(subIDs)
        if ~isempty(EEGAxx{s})
            % --------------PLOT: interactive head and spectrum plots-----------
            SDEEG{s} = EEGAxx{s}.Cos+(EEGAxx{s}.Sin*1i);%EEGAxx{s}.Amp;% it is important which n is considered for fft
            if strcmp(subjects,'Individuals'),
                if strcmp(SignalType,'Amplitude')
                    H = mrC.Simulate.PlotEEG(abs(SDEEG{s}),freq,SavePath,subIDs{s},RoiList,signalFF,SignalType,[],Mode);% Plot individuals
                elseif strcmp(SignalType,'Phase')
                    H = mrC.Simulate.PlotEEG(wrapTo2Pi(angle(SDEEG{s})),freq,SavePath,subIDs{s},RoiList,signalFF,SignalType,[],Mode);% Plot individuals
                end
                if ~isempty(H)
                    set(H,'PaperPositionMode','manual');
                    set(H,'units','centimeters');
                    set(H, 'PaperPosition',[1 1 12 5]);   
                    print(fullfile(SavePath,[ 'SimEEG_Subject_' subIDs{s} '_Electrode_Freq' num2str(signalFF) 'Hz_' SignalType '_' SaveName '.tif']),'-dtiff','-r300');
                end
            end 
        end 
    end
    
    % Plot average over individuals
    MSDEEG = mean(cat(4,SDEEG{:}),4);
    
    if strcmp(SignalType,'Amplitude')
        h = mrC.Simulate.PlotEEG(abs(MSDEEG),freq,SavePath,'average over all  ',RoiList,signalFF,SignalType,[],Mode);
    
    elseif strcmp(SignalType,'Phase')    
        h = mrC.Simulate.PlotEEG(wrapTo2Pi(angle(MSDEEG)),freq,SavePath,'average over all  ',RoiList,signalFF,SignalType,[],Mode);
        
    end
    if ~isempty(h)
        
        set(h,'PaperPositionMode','manual')
        set(h,'units','centimeters')
        set(h, 'PaperPosition',[1 1 12 5]);    
        print(fullfile(SavePath,[ 'SimEEG_Subject_Average_Electrode_Freq' num2str(signalFF(1)) 'Hz_' SignalType '_' SaveName  '.tif']),'-dtiff','-r300');
    end
end