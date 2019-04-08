% estimate and write out functional model and  parameters for
% distance-dependence of coherence (according to Fig. 4 from 'Multi-scale
% analysis of neural activity in humans: Implications for 
% micro-scale electrocorticography'
clear all

model_fun.expo = @(p,x)( (1) *exp(-p(1)*x) );
% model_fun.lorentzian = @(p,x)( (1)*(p(1)^2./(p(1)^2+x.^2))    ) ;
% model_fun.gauss = @(p,x) ( (1) * exp(-(x.^2) )/p(1)  ) ;
% model_fun.power_law = @(p,x) ( (p(2)*x+1).^-p(1)) ;

% % data taken from Fig. 4
% x = [0,0.38,0.5,0.54,1:0.5:3] ;% in mm
% y.delta = [7., 6.3, 6.1, 6.05, 5.5,  5.3,  5,    4.9,  4.8]/7; % strange number where read from print in cm
% y.theta = [7., 6.1, 5.8, 5.7,  5.1,  4.8,  4.5,  4.4,  4.3]/7;
% y.beta =  [7., 6.3, 6.1, 6.05, 5.75, 5.5,  5.25, 5.1,  5]/7;
% y.alpha = [7., 5.1, 4.9, 4.8,  4.45, 4.2,  3.9,  3.8,  3.6]/7;
% y.gamma = [7., 3.2, 2.9, 2.75, 2.3,  2.05, 1.9,  1.75, 1.6]/7;

% from microEcog3, fig. 3 in supplementary material
% x =       [0., 0.38, 0.5,   0.54, 1.,   1.5,  2.,   2.5,  3.,   3.5,  4.,   4.5,  5.] ;% in mm
% y.delta = [1., 0.85, 0.825, 0.81, 0.72, 0.67, 0.62, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61] ;
% y.theta = [1., 0.9,  0.87,  0.85, 0.8,  0.73, 0.7,  0.67, 0.67, 0.67, 0.67, 0.67, 0.67] ;
% y.beta =  [1., 0.73, 0.7,   0.68, 0.63, 0.6,  0.55, 0.63, 0.5,  0.48, 0.46, 0.44, 0.4 ];
% y.alpha = [1., 0.9,  0.87,  0.81, 0.82, 0.79, 0.75, 0.725, 0.71, 0.68, 0.66, 0.62, 0.6 ] ;
% y.gamma = [1., 0.45, 0.42,  0.39, 0.33, 0.3 , 0.27, 0.25,  0.22, 0.21,  0.2,  0.2, 0.19];

%% data taken from micro ECoG3 in suplimentary material
x       = [ 0., 0.5,  1.,   1.5,  2.,   2.5,  3.,  3.5,  4.] ;
y.delta = [ 1., 0.96, 0.94, 0.9,  0.88, 0.84, 0.8, 0.78, 0.76 ];
y.theta = [ 1., 0.96, 0.94, 0.9,  0.88, 0.84, 0.82, 0.8, 0.78 ];
y.alpha = [ 1., 0.95, 0.91, 0.86, 0.84, 0.78, 0.72, 0.67, 0.61];
y.beta  = [ 1., 0.92, 0.84, 0.76, 0.68, 0.61, 0.59, 0.56, 0.48] ;
y.gamma = [ 1., 0.84, 0.7,  0.68, 0.48, 0.4, 0.36, 0.32, 0.28];

band_freqs.delta = [0,4];
band_freqs.theta = [4,8];
band_freqs.alpha = [8,12];
band_freqs.beta = [12,30];
band_freqs.gamma = [30,80];
%%
FigPath = '../../../Examples/Figures';
fig_all_fits=figure ;
freq_band_names = fieldnames(y) ;
model_names = fieldnames(model_fun) ;
FS = 12; % fontsize
x_temp = [0:0.5:70];

colors = {'r','g','b','k','c','m','y'};
linestyles = {'-','--',':','-.'};

handles = [];
for i = 1:length(freq_band_names)
    this_bands_mse = 100000;
    this_band_name = freq_band_names{i};
    for this_model_idx = 1:length(model_names)
        subplot(1,length(model_names),this_model_idx );
        this_model_name = model_names{this_model_idx} ;
        scatter(x,y.(freq_band_names{i}),colors{i})
        hold on 
        if strcmp(this_model_name,'power_law')
             if strfind(version,'R2011')
             [this_p,R,J,CovB,MSE] =nlinfit(x(1:end),y.(this_band_name)(1:end),model_fun.(this_model_name),[1.,1.,1.   ,0.]);
             else
             [this_p,R,J,CovB,MSE,ErrorModelInfo] =nlinfit(x(1:end),y.(this_band_name)(1:end),model_fun.(this_model_name),[1.,1.,1.   ,0.]);
             end
        else
            if strfind(version,'R2011')
            [this_p,R,J,CovB,MSE] =nlinfit(x(1:end),y.(this_band_name)(1:end),model_fun.(this_model_name),[1.,1.,0  ,0.]);
            else
            [this_p,R,J,CovB,MSE,ErrorModelInfo] =nlinfit(x(1:end),y.(this_band_name)(1:end),model_fun.(this_model_name),[1.,1.,0  ,0.]);
            end
        end
        sprintf('band: %s, model: %s, MSE=%f, mean(R)= %f',this_band_name,this_model_name, MSE,mean(R.^2))
        model_params.(this_band_name).(this_model_name) = this_p;
        h=plot(x_temp,model_fun.(this_model_name)( model_params.(this_band_name).(this_model_name) ,x_temp),colors{i},'linewidth',2);
        xlim([0,50])
        if this_model_idx==1
            handles = [handles,h];
        end
        
        if MSE<this_bands_mse
            this_bands_mse = MSE ;
            best_model.(this_band_name).model_params = this_p;
            best_model.(this_band_name).fun = model_fun.(this_model_name);
            best_model.(this_band_name).model_name = this_model_name ;
        end
        
        
        %title(this_model_name)
    end

end
set(gca,'fontsize',FS,'ytick',0:.2:1);
xlabel('Spatial Distance [mm]','fontsize',FS);
ylabel('Coherence','fontsize',FS);
legend(handles,freq_band_names,'fontsize',FS);

set(fig_all_fits,'paperposition',[10 10 4 3]);
set(fig_all_fits,'Unit','Inch','position',[10 10 4 3],'color','w');

print(fullfile(FigPath,'CoherenceFitFunction.tif'),'-r300','-dtiff');
export_fig(fig_all_fits,fullfile(FigPath,'CoherenceFitFunction'),'-pdf');

save('spatial_decay_models_coherence','model_fun','model_params','best_model','band_freqs')
%%
% plot best model
if false
    addpath('../../../External/tools/BrewerMap/')
    set(0,'defaultfigurecolor',[1 1 1])

    plot_dir = fullfile('/Users/bosse/dev/matlab/mrC/Examples/plots') ;
    if ~exist(plot_dir)
        mkdir(plot_dir)
    end

    sorted_freq_band_names = sort(freq_band_names);

    alw = 0.75;    % AxesLineWidth
    fsz = 11;      % Fontsize
    lw = 1.5;      % LineWidth
    msz = 8;       % MarkerSize

    x_temp = 0:0.01:70 ;
    colors = brewermap(length(freq_band_names),'Set2');
    fig = figure;
    legend_entries = {} ;

    for i = 1:length(sorted_freq_band_names )
        legend_entries{end+1} = sprintf('\\%s-band',sorted_freq_band_names {i});
        plot(x_temp,best_model.(sorted_freq_band_names {i}).fun(...
            best_model.(sorted_freq_band_names {i}).model_params ,x_temp),'color',colors(i,:),'linewidth',lw);
        hold on
    end
    xlabel('Spatial Distance [mm]')
    ylabel('Coherence')
    ylim([0,1])
    xlim([0,70])
    legend(legend_entries,'FontSize',fsz,'location','northeastoutside')
    set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    fig_width = 8.5/2 ;
    fig_height = 5/8*fig_width ;
    fig.PaperPositionMode = 'manual';
    fig.PaperUnits = 'inches';
    fig.Units = 'inches';
    fig.PaperPosition = [0, 0, fig_width,fig_height];
    fig.PaperSize = [fig_width,fig_height];
    fig.Position = [0.0, 0.0, fig_width, fig_height];
    fig.Resize = 'off';
    fig.InvertHardcopy = 'off';
    filename = fullfile(plot_dir,'spatial_decay_coherence');
    print(filename,'-dpdf');
end
