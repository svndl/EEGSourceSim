function steadyStateWaveMovie(varargin)    
    % Description:	generate steady-state waveform movies
    % 
    % Syntax:	steadyStateWaveMovie(<options>)
    % <options>
    %   freq    - 1 x n array of doubles indicating the frequencies to plot
    %    
    %   colors    -  n x 3 array of doubles indicating the colors to use
    % 
    %   slowdown    - scalar indicating how much to slow down the movie
    %                   t means t x slower [10]
    %
    %   sampleRate    - scalar indicating temporal sampling rate of movie [10 Hz]
    %
    %   outFormat    - string indicating whether to output gif (['gif']) or
    %                  avi movie ('avi')
    %
    %   filename    - string indicating the output path and filename of the
    %                saved file
    
    setenv('DYLD_LIBRARY_PATH','')

    %% PARSE ARGS
    opt	= ParseArgs(varargin,...
            'freq', [5,6,1], ...
            'colors', [], ...
            'slowdown', 10, ...
            'sampleRate', 20, ...
            'movieDur',[], ...
            'makeMovie',true, ...
            'outFormat', 'gif', ...
            'figAspectRatio',  [4,2], ...
            'phaseDeg',  0, ...
            'printSingleWaves', false, ...
            'filename',  '~/Desktop/test' ...
            );
    
    if isempty(opt.colors)
        cBrewer = load('colorBrewer');
        opt.colors = cBrewer.rgb10; % repeat disparity color
    else
    end
    
    if isempty(opt.movieDur)
        % movie duration in seconds, 
        % x cycles of smallest frequency multiplied by the slowdown factor
        opt.movieDur = (1:10)./min(opt.freq);
        % at least one second long
        opt.movieDur = min(opt.movieDur(opt.movieDur>=1));
    else
    end
    
%     if isempty(opt.slowdown)
%         if opt.makeMovie
%             opt.slowdown = 10;
%         else
%             opt.slowdown = 1;
%         end
%     else
%     end
    
    % total number of samples
    numSamples = opt.movieDur * opt.slowdown * opt.sampleRate;  
    
    if ~ismember(opt.outFormat,{'gif','avi'})
        error('unknown format %s',opt.outFormat);
    else
    end
    
    % get rid of suffix, if user added it
    if strcmp(opt.filename(end-3:end),['.',opt.outFormat])
        opt.filename = opt.filename(1:end-4);
    else
    end
    
    % compute the sine waves
    numFreq = length(opt.freq);
    amp = 1;
    T = linspace(0,opt.movieDur,numSamples+1);
    T = T(1:numSamples);
    sinFxn = zeros(numSamples,numFreq);
    if numel(opt.phaseDeg) < numFreq
        opt.phaseDeg = repmat(opt.phaseDeg(1),1,numFreq);
    else
    end
    for f = 1:numFreq
        circFreq = 2*pi*opt.freq(f);
        sinFxn(:,f) = amp.*sin( circFreq.*T + opt.phaseDeg(f) * pi/180);
    end
    
    % draw it
    lWidth = 10;
    gcaOpts = {'tickdir','out','ticklength',[0.0,0.0],'box','off','linewidth',lWidth};

    for f = 1:numFreq
        sinFig = figure;
        figPos = get(sinFig,'pos');
        figPos(3) = figPos(4)*opt.figAspectRatio(1);
        figPos(4) = figPos(4)*opt.figAspectRatio(2);
        set(sinFig,'pos',figPos);
        hold on
        for z = 1:numSamples
            plotH = plot(T(1:z),sinFxn(1:z,f),'color','k','linewidth',lWidth);
            set(gca,gcaOpts{:},'xtick',[],'ytick',[], 'visible', 'off','Color',[1 1 1]);
            warning('off','all')
            xlim([-.05,opt.movieDur+.05]);
            ylim([-1.05,1.05]);
            drawnow
            pause(.1);
            frame = getframe(sinFig);
            im = frame2im(frame); 
            A = rgb2ind(im,256);
            A = A(5:(end-5),5:(end-5)); % crop outer edges
            bgIdx = mode(A); 
            lineIdx = mode(A(A ~= bgIdx)); % assume that second most often value is line
            curA = uint8(zeros( size(A) ));
            curA(  A == lineIdx ) = f;
            movieFrames(:,:,f,z) = curA;
            delete(plotH);
        end
        close(sinFig);
    end
        
    movieMap(1,:) = [1 1 1]; % make background white (will be transparent in gif) 
    movieMap(2:numFreq+1,:) = opt.colors(1:numFreq,:);
    movieFrames(:,:,end+1,:) = squeeze(max(movieFrames,[],3));
   
    % crop images
    lastFrame = movieFrames(:,:,end,end);
    colIdx = find(sum(lastFrame == 0,1) ~= size(lastFrame,1));
    rowIdx = find(sum(lastFrame == 0,2) ~= size(lastFrame,2));    
    movieFrames = movieFrames(rowIdx,colIdx,:,:);
    
    for n = 1:(numFreq+1)
        if n == (numFreq+1)
            curFilename = opt.filename;
            % always print full set of waveforms
            printWave = true;
        else
            curFilename = sprintf('%s_%0.3fHz',opt.filename,opt.freq(n));
            if opt.printSingleWaves
                printWave = true;
            else
                printWave = false;
            end
        end
        if printWave
            % write last frame as single image
            if exist([curFilename,'_lastframe.gif'],'file')
                delete([curFilename,'gif']);
            end
            imwrite(movieFrames(:,:,n,end),movieMap,[curFilename,'_lastframe.gif'],'transparentColor',0);
            if opt.makeMovie
                if strcmp(opt.outFormat,'gif')
                    if exist([curFilename,'gif'],'file')
                        delete([curFilename,'gif']);
                    end
                    for k = 1:numSamples
                        if k==1
                            gifOpts = {'LoopCount',Inf};
                        else
                            gifOpts = {'WriteMode','append'};
                        end
                        imwrite(movieFrames(:,:,n,k),movieMap,[curFilename,'.gif'],'gif','DelayTime',1/opt.sampleRate,gifOpts{:},'transparentColor',0);
                    end
                else
                    if exist([opt.filename,'avi'],'file')
                        delete([curFilename,'avi']);
                    end
                    vidObj = VideoWriter([curFilename,'.avi'],'Indexed AVI');
                    vidObj.FrameRate = opt.sampleRate; % frames per second
                    vidObj.Colormap = movieMap;
                    open(vidObj);
                    for k = 1:size(movieFrames, 4)
                       % Write each frame to the file.
                       writeVideo(vidObj,movieFrames(:,:,n,k));
                    end
                    close(vidObj);
                end
            else
            end
        else
        end
    end
    warning('on','all')
end
