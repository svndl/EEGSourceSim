% ADD required toolboxes
clear; clc;

SimFolder = fileparts(pwd);
addpath(genpath(SimFolder));
%% Prepare the results folders
FigPath = 'Figures';
ResultPath = 'ResultData';
if ~exist(fullfile(pwd,FigPath),'dir'),mkdir(FigPath);end
if ~exist(fullfile(pwd,ResultPath),'dir'),mkdir(ResultPath);end

%% Prepare Project path and ROIs
DestPath = fullfile(SimFolder,'Examples','ExampleData');
AnatomyPath = fullfile(DestPath,'anatomy');
ProjectPath = fullfile(DestPath,'FwdProject');
[Inverse,subIDs_Inverse] = mrC.Simulate.ReadInverses(ProjectPath,'mneInv_bem_gcv_regu_TWindow_0_1334_wangROIsCorr.inv');
subIDs_Inverse = subIDs_Inverse(cellfun(@(x) ~isempty(x),Inverse));
clear Inverse;

% Pre-select ROIs
[RoiList,subIDs] = mrC.Simulate.GetRoiClass(ProjectPath,AnatomyPath,subIDs_Inverse);% 13 subjects with Wang atlab 
Wang_RoiList = cellfun(@(x) {x.getAtlasROIs('wang')},RoiList);
Inds = 1:50;Inds([1:14 17:24 27:28])= [];
Wang_RoiList = cellfun(@(x) {x.selectROIs(Inds)},Wang_RoiList);
Ind2 = [5:18 1:4 19:22 25:26];
Wang_RoiList = cellfun(@(x) {x.selectROIs(Ind2)},Wang_RoiList);

%% Generate Resolution matrices
FilePath = fullfile(ResultPath,'LocalizationExampleData_Paper.mat');
do_new_data_generation = true;

if ~exist(FilePath,'file') || do_new_data_generation
    [CrossTalk1,Error1,ROISource1,~,~,~] = mrC.Simulate.ResolutionMatrices(ProjectPath,'subSelect',subIDs,...
        'rois',Wang_RoiList,'roiType','wang','anatomyPath',AnatomyPath,'doAUC',true,'inverse','mneInv_bem_snr_100.inv');
    
    [CrossTalk2,Error2,ROISource2,ScalpData,LIST,subIDs] = mrC.Simulate.ResolutionMatrices(ProjectPath,'subSelect',subIDs,...
        'rois',Wang_RoiList,'roiType','wang','anatomyPath',AnatomyPath,'doAUC',true,'inverse','mneInv_bem_gcv_regu_TWindow_0_1334_wangROIsCorr.inv');
    save(FilePath,'CrossTalk1','Error1','ROISource1','CrossTalk2','Error2','ROISource2','ScalpData','LIST','subIDs');
else
    load(FilePath);
end

LIST2 = cellfun(@(x) x(6:end-2),LIST,'uni',false);
LIST2 = strrep(LIST2,'_','-');
LIST2 = LIST2(1:2:end);
% Preparing the color
if ~exist(fullfile(ResultPath,'ROI_colors_Paper.mat'),'file')
    Colors = distinguishable_colors(numel(LIST)/2);
    CC = zeros(size(Colors,1)*2,3);% for two hemisphere or two inverses
    CC(1:2:end)=Colors;
    CC(2:2:end)=Colors+repmat(-[.1 .1 .1],[size(Colors,1) 1]);CC(CC<0)=0;
    save(fullfile(ResultPath,'ROI_colors_Paper.mat'),'Colors','CC');
else
    load(fullfile(ResultPath,'ROI_colors_Paper.mat'));
end

% prepare xticklabel data
for i = 1:numel(LIST2)
    XTLabel{i} = ['\color[rgb]{' num2str(Colors(i,1)) ',' num2str(Colors(i,2)) ',' num2str(Colors(i,3)) '}' LIST2{i}];
end


%% Plot ROIs, Cross Talk and 
BrainFromSuma = true;
FIG1 = figure;
set(FIG1,'units','normalized','color',[.8 .8 .8]);
FS = 11;
INDs = Inds(Ind2);

if ~BrainFromSuma
    S1 = subplot(2,3,1); mrC.Simulate.VisualizeSourceRoi2('nl-0048',[AnatomyPath '/'],'wang',INDs,'anterior','B',CC,[]);
    set(S1,'position',get(S1,'position')+[-.025 -.05 .05 .05]);
    subplot(2,3,4), mrC.Simulate.VisualizeSourceRoi2('nl-0048',[AnatomyPath '/'],'wang',INDs,'ventral','B',CC,1:2:numel(INDs));
else
    im1 = imread(fullfile('private','SumaBrain','eb_dorsolateral_pial_crop.png'));
    im2 = imread(fullfile('private','SumaBrain','eb_ventromedial_pial_crop.png'));
    S1 = subplot(3,2,1);imagesc(im1(:,:,1:3),'AlphaData',squeeze(sum(im1,3)~=0));axis off %equal   
    S2 = subplot(3,2,2);imagesc(im2(:,:,1:3),'AlphaData',squeeze(sum(im2,3)~=0));axis off %equal
    set(S1,'position',get(S1,'position')+[-.04 -.02 .05 .04]);
    set(S2,'position',get(S2,'position')+[-.04 -.02 .05 .04]);
    
    im1 = imread(fullfile('private','SumaBrain','eb_dorsolateral_inf_crop.png'));
    im2 = imread(fullfile('private','SumaBrain','eb_ventromedial_inf_crop.png'));
    S1 = subplot(3,2,3);imagesc(im1(:,:,1:3),'AlphaData',squeeze(sum(im1,3)~=0));axis off %equal   
    S2 = subplot(3,2,4);imagesc(im2(:,:,1:3),'AlphaData',squeeze(sum(im2,3)~=0));axis off %equal
    set(S1,'position',get(S1,'position')+[-.04 -.02 .05 .04]);
    set(S2,'position',get(S2,'position')+[-.04 -.02 .05 .04]);
end

% Plot Cross Talk Matrices
CrossTalk1 = cellfun(@(x) x./repmat(max(x),[size(x,1) 1]),CrossTalk1,'uni',false);% normalize
CrossTalk2 = cellfun(@(x) x./repmat(max(x),[size(x,1) 1]),CrossTalk2,'uni',false);% normalize
CT1 = (cat(3,CrossTalk1{:}));CT1 = (CT1(1:2:end,1:2:end,:)+CT1(2:2:end,2:2:end,:))./2;
CT2 = (cat(3,CrossTalk2{:}));CT2 = (CT2(1:2:end,1:2:end,:)+CT2(2:2:end,2:2:end,:))./2;

for i = 1:2
    eval(['CTM1 = mean(cat(3,CrossTalk' num2str(i) '{:}),3);']);
    S = subplot(3,2,i+4);
    eval(['imagesc(abs(mean(CT' num2str(i) ',3)));']);
    colormap(jmaColors('coolhot'));
    colormap('gray')
    %caxis([-max(abs(CTMM(:))) max(abs(CTMM(:)))]);
    caxis([0 1]);
    
    set(gca,'ytick',1:numel(LIST2),'yticklabel',XTLabel,'xtick',1:numel(LIST),'xticklabel',XTLabel,'fontsize',FS-2);
    if exist('xtickangle'), xtickangle(90); end
    xlabel('receiving ROI','fontsize',FS);
    for r1 = 1:12
        for r2 = 1:12
            [~,p_CT(r1,r2)]=ttest(CT1(r1,r2,:),CT2(r1,r2,:));
            if p_CT(r1,r2)<(0.01)%/(12*11))
                if i==1 
                    if abs(mean(CT1(r1,r2,:)))<abs(mean(CT2(r1,r2,:)))
                        text(r2-.2,r1+.2,'*','fontsize',15,'color','r')
                    end
                else
                    if abs(mean(CT1(r1,r2,:)))>abs(mean(CT2(r1,r2,:)))
                        text(r2-.2,r1+.2,'*','fontsize',15,'color','r')
                    end
                end
            end
        end
    end
    clear p_CT;
    if i==1
        ylabel('Seed ROI','fontsize',FS);
        T = title('Minimum Norm','fontsize',FS,'fontweight','bold');
        set(S,'position',get(S,'position')+[0.02 0 0 0]);
    else
        T = title('Minimum Norm + FACE','fontsize',FS,'fontweight','bold');
        CB = colorbar;
        set(S,'position',get(S,'position')+[0 0 0 0]);
        
        set(CB,'position',get(CB,'position')+[-.012 0 -.015 0],'Ticks',0:.5:1);
        
    end
    set(T,'position',get(T,'position')+[0 -.5 0])
    
    
end

%
set(FIG1,'PaperPosition',[1 1 5 7]);
print(fullfile(FigPath,'SourceEstimation1'),'-r300','-dtiff')
set(gcf, 'Color', 'w');
set(FIG1,'Units','Inch')
FIG1.Position = [1 1 5 8];
export_fig(FIG1,fullfile(FigPath,'SourceEstimation1'),'-pdf')

%% AUC relative Enegry, Focalization Error

% prepare data
ROINum = Wang_RoiList{1,1}.ROINum;
for i = 1:2
    eval(['Error = Error' num2str(i) ';']);
    for sub = 1:numel(subIDs)
        TP = Error{sub}.TP;
        FP = Error{sub}.FP;
        FN = Error{sub}.FN;
        Recall = TP./(TP+FN);
        Precision = TP./(TP+FP);
        for r = 1:size(TP,1)
            AUC(r,sub) = trapz(Recall(r,:),Precision(r,:));
        end
        Relative(:,sub) = Error{sub}.Relative;
        Focal(:,sub) = Error{sub}.Focalization;
    end
    
    AUC2(:,:,i) = (AUC(1:2:end,:)+AUC(2:2:end,:))/2; % Average over left and right
    Focal2(:,:,i) = (Focal(1:2:end,:)+Focal(2:2:end,:))/2; % Average over left and right
    Relative2(:,:,i) = (Relative(1:2:end,:)+Relative(2:2:end,:))/2; % Average over left and right
end
% Stats

for r = 1:size(Relative2,1)
    [~,p_AUC(r)] = ttest(AUC2(r,:,1),AUC2(r,:,2));
    [~,p_Relative(r)] = ttest(Relative2(r,:,1),Relative2(r,:,2));
    [~,p_Focal(r)] = ttest(Focal2(r,:,1),Focal2(r,:,2));
end

FS = 11;
FIG3 = figure;
C = [.6 .6 .6];%Colors;%

% AREA UNDER CURVE
AUCM = squeeze(mean(AUC2,2));
AUCS = squeeze(std(AUC2,[],2))./sqrt(size(AUC2,2));
S1 = subplot(3,1,1); 
B = bar(AUCM);
set(B(1),'FaceColor',[0.3,.3,.3])
set(B(2),'FaceColor',[.7,.7,.7])
set(gca,'xtick',1:numel(LIST2)*2,'xticklabel',XTLabel,'fontsize',FS-1,'ytick',0.2:.2:1,'yticklabel',arrayfun(@(x) num2str(x),0.2:.2:1,'uni',false));
ylim([0.0 0.75]);
hold on;
XT = [(1:size(AUCM,1))-.15; (1:size(AUCM,1))+.15];
errorbar(reshape(XT,1,numel(XT)),reshape(AUCM',1,numel(AUCM)),reshape(AUCS',1,numel(AUCS)),'.k');
YL = ylabel('AUCPR','fontsize',FS,'fontweight','bold'); YLP = get(YL,'position');
set(YL,'position',get(YL,'position')+[-0.2 0 0]); 
set(S1,'position',get(S1,'position')+[0 0 0 .04]);
% plot significant stars
for r = 1:size(XT,2)
    if p_AUC(r)<0.05
        line(XT(:,r),[max(AUCM(r,:)+AUCS(r,:)) max(AUCM(r,:)+AUCS(r,:))]+0.04,'color','k','linewidth',1)
        text(r-.08,max(AUCM(r,:)+AUCS(r,:))+.06,'*')
    end
end
legend({'MN','MN + FACE'},'location','northwest')

% RELATIVE ENERGY
RelativeM = squeeze(mean(Relative2,2));
RelativeS = squeeze(std(Relative2,[],2))./sqrt(size(Relative2,2));
S2 = subplot(3,1,2); 
B = bar(RelativeM);
set(B(1),'FaceColor',[0.3,.3,.3])
set(B(2),'FaceColor',[.7,.7,.7])
set(gca,'xtick',1:numel(LIST2)*2,'xticklabel',XTLabel,'fontsize',FS-1,'ytick',0:.1:.3,'yticklabel',arrayfun(@(x) num2str(x),0:.1:.3,'uni',false));
ylim([0 .25]);
hold on;
XT = [(1:size(RelativeM,1))-.15; (1:size(RelativeM,1))+.15];
errorbar(reshape(XT,1,numel(XT)),reshape(RelativeM',1,numel(RelativeM)),reshape(RelativeS',1,numel(RelativeS)),'.k');
YL = ylabel('Relative Energy','fontsize',FS,'fontweight','bold');
YLPT = get(YL,'position');
set(YL,'position',[YLP(1) YLPT(2:end)]); 
set(YL,'position',get(YL,'position')+[-0.15 0 0]);
set(S2,'position',get(S2,'position')+[0 -.02 0 .04]);
% plot significant stars
for r = 1:size(XT,2)
    if p_Relative(r)<0.05
        line(XT(:,r),[max(RelativeM(r,:)+RelativeS(r,:)) max(RelativeM(r,:)+RelativeS(r,:))]+0.01,'color','k','linewidth',1)
        text(r-.08,max(RelativeM(r,:)+RelativeS(r,:))+.015,'*')
    end
end


% FOCALIZATION ERROR
FocalM = squeeze(mean(Focal2,2));
FocalS = squeeze(std(Focal2,[],2))./sqrt(size(Focal2,2));
S3 = subplot(3,1,3); 
B = bar(FocalM);
set(B(1),'FaceColor',[0.3,.3,.3])
set(B(2),'FaceColor',[.7,.7,.7])
set(gca,'xtick',1:numel(LIST2)*2,'xticklabel',XTLabel,'fontsize',FS-1,'ytick',0.2:.2:1,'yticklabel',arrayfun(@(x) num2str(x),0.2:.2:1,'uni',false));
ylim([0.2 .9]);
hold on;
XT = [(1:size(FocalM,1))-.15; (1:size(FocalM,1))+.15];
errorbar(reshape(XT,1,numel(XT)),reshape(FocalM',1,numel(FocalM)),reshape(FocalS',1,numel(FocalS)),'.k');
YL = ylabel('Focalization Error','fontsize',FS,'fontweight','bold');
YLPT = get(YL,'position');
set(YL,'position',[YLP(1) YLPT(2:end)]); 
set(S3,'position',get(S3,'position')+[0 -.04 0 .04]);
% plot significant stars
for r = 1:size(XT,2)
    if p_Focal(r)<0.05
        line(XT(:,r),[max(FocalM(r,:)+FocalS(r,:)) max(FocalM(r,:)+FocalS(r,:))]+0.04,'color','k','linewidth',1)
        text(r-.08,max(FocalM(r,:)+FocalS(r,:))+.06,'*')
    end
end

set(FIG3,'PaperPosition',[1 1 5 5]);
print(fullfile(FigPath,'SourceEstimation2'),'-r300','-dtiff')
set(gcf, 'Color', 'w');
set(FIG3,'Units','Inch')
set(FIG3,'Position',[1 1 6 6]);
export_fig(FIG3,fullfile(FigPath,'SourceEstimation2'),'-pdf')
