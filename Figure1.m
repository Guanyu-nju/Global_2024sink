%% Mask create
clc
clear
% Global mask
pixel_mask=importdata("E:\phd_file\Tropical_2024\Tropical_mask\RECCAP_mask.tif");


% 创建环境
data=geotiffread('E:\phd_file\Boreal_North_America\region_lu.tif');
info=geotiffinfo('E:\phd_file\Boreal_North_America\region_lu.tif');
[m,n] = size(data);
k=1;
for i = 1:m
    for j = 1:1
        [lat,lon]= pix2latlon(info.RefMatrix, i, j);   %读取栅格数据第1列所有行的纬度；
        lat_(k,:)=lat; %将纬度数据存储为1列；
        k=k+1;
    end
end

k=1;
for ii = 1:1
    for jj = 1:n
        [lat,lon]= pix2latlon(info.RefMatrix, ii, jj);   %读取栅格数据第1行所有行的经度；
        lon_(k,:)=lon;  %将经度数据存储为1列；
        k=k+1;
    end
end
[lon1,lat1]=meshgrid(lon_,lat_);

Boundry = shaperead("E:\phd_file\Tropical_2024\Tropical_mask\Land\arcgis_globe.shp");
bou_canX = [Boundry(:).X];
bou_canY= [Boundry(:).Y];
clearvars -except lon1 lat1 pixel_mask bou_canX bou_canY m n

Global_background=pixel_mask;
Global_background(Global_background<1 | Global_background>10)=nan;
Global_background(Global_background==8)=5;
Global_background(Global_background==9)=8;
Global_background(Global_background==10)=9;

LUCC_mask=importdata("E:\phd_file\Tropical_2024\land_cover\MODIS\MCD12C1_LUCC_2023.tif");
LUCC_mask=double(LUCC_mask);

forest_mask=LUCC_mask;
forest_mask(forest_mask==1 | forest_mask==2 |forest_mask==3 |forest_mask==4 |forest_mask==5)=666;
forest_mask(forest_mask~=666)=nan;
forest_mask(forest_mask==666)=1;

Non_forest_mask=LUCC_mask;
Non_forest_mask(Non_forest_mask==10 | Non_forest_mask==12 | Non_forest_mask==14 | ...
    Non_forest_mask==6 | Non_forest_mask==7 | Non_forest_mask==8 | Non_forest_mask==9)=666;
Non_forest_mask(Non_forest_mask~=666)=nan;
Non_forest_mask(Non_forest_mask==666)=1;

LUCC_mask(forest_mask==1 | Non_forest_mask==1)=1;
LUCC_mask(LUCC_mask~=1)=nan;

pixel_mask=LUCC_mask.*pixel_mask;

North_America_mask=pixel_mask;
North_America_mask(North_America_mask~=1)=nan;

South_America_mask=pixel_mask;
South_America_mask(South_America_mask~=2)=nan;
South_America_mask(~isnan(South_America_mask))=1;

Russia_mask=pixel_mask;
Russia_mask(Russia_mask~=3)=nan;
Russia_mask(~isnan(Russia_mask))=1;

Europe_mask=pixel_mask;
Europe_mask(Europe_mask~=4)=nan;
Europe_mask(~isnan(Europe_mask))=1;

WS_Asia_mask=pixel_mask;
WS_Asia_mask(WS_Asia_mask~=5 & WS_Asia_mask~=8)=nan;
WS_Asia_mask(~isnan(WS_Asia_mask))=1;

Africa_mask=pixel_mask;
Africa_mask(Africa_mask~=6)=nan;
Africa_mask(~isnan(Africa_mask))=1;

East_Asia_mask=pixel_mask;
East_Asia_mask(East_Asia_mask~=7)=nan;
East_Asia_mask(~isnan(East_Asia_mask))=1;

SE_Asia_mask=pixel_mask;
SE_Asia_mask(SE_Asia_mask~=9)=nan;
SE_Asia_mask(~isnan(SE_Asia_mask))=1;

Oceania_mask=pixel_mask;
Oceania_mask(Oceania_mask~=10)=nan;
Oceania_mask(~isnan(Oceania_mask))=1;


pixel_mask(pixel_mask<1 | pixel_mask>10)=nan;
pixel_mask(~isnan(pixel_mask))=1;

Tropical_mask=pixel_mask;
Tropical_mask(1:67,:)=nan;
Tropical_mask(114:180,:)=nan;

Boreal_mask=pixel_mask;
Boreal_mask(40:180,:)=nan;

North_temperate_mask=pixel_mask;
North_temperate_mask(1:39,:)=nan;
North_temperate_mask(68:180,:)=nan;

South_temperate_mask=pixel_mask;
South_temperate_mask(1:113,:)=nan;
South_temperate_mask(151:180,:)=nan;

Global_mask=pixel_mask;

area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;
%%
for year=2015:2024

    % NBE
    BEPS_GFAS_landflux=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\BEPS_GFAS\annual\NBE_BEPS_GFAS_' num2str(year) '.tif']);
    BEPS_GFED_landflux=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\BEPS_GFED\annual\NBE_BEPS_GFED_' num2str(year) '.tif']);
    CASA_GFAS_landflux=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\CASA_GFAS\annual\NBE_CASA_GFAS_' num2str(year) '.tif']);
    CASA_GFED_landflux=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\CASA_GFED\annual\NBE_CASA_GFED_' num2str(year) '.tif']);
    Mean_NBE=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\mean_value\annual\NBE_' num2str(year) '.tif']);
 
    % All fire
    Mean_total_fire=importdata(['E:\phd_file\Tropical_2024\fire emission\mean_all_carbon_value\annual\Fire_' num2str(year) '.tif']);

    % GPP
    GOSIF_GPP=importdata(['E:\phd_file\Tropical_2024\GPP\GOSIF\yearly\1degree\GOSIF_GPP_' num2str(year) '.tif']);
    FluxSat_GPP=importdata(['E:\phd_file\Tropical_2024\GPP\FluxSat\yearly\FluxSat_GPP_' num2str(year) '.tif']);
    Mean_GPP=importdata(['E:\phd_file\Tropical_2024\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    Mean_GPP_uncertainty=importdata(['E:\phd_file\Tropical_2024\GPP\mean_value\uncertainty\yearly\GPP_uncertainty_' num2str(year) '.tif']);
    % TER
    Mean_TER=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\TER\annual\TER_' num2str(year) '.tif']);   

    % landflux list
    Mean_landflux_list(year-2014)=nansum(nansum(Mean_NBE.*area_grid/(10^15)));
    GCAS_landflux_list(1,year-2014)=nansum(nansum(BEPS_GFAS_landflux.*area_grid/(10^15)));
    GCAS_landflux_list(2,year-2014)=nansum(nansum(BEPS_GFED_landflux.*area_grid/(10^15)));
    GCAS_landflux_list(3,year-2014)=nansum(nansum(CASA_GFAS_landflux.*area_grid/(10^15)));
    GCAS_landflux_list(4,year-2014)=nansum(nansum(CASA_GFED_landflux.*area_grid/(10^15)));
    
   % Fire list
    Mean_fire_list(year-2014)=nansum(nansum(Mean_total_fire.*area_grid/(10^15)));


    % GPP list
    GOSIF_GPP_list(year-2014)=nansum(nansum(GOSIF_GPP.*area_grid/(10^15)));
    FluxSat_GPP_list(year-2014)=nansum(nansum(FluxSat_GPP.*area_grid/(10^15)));
    Mean_GPP_list(year-2014)=nansum(nansum(Mean_GPP.*area_grid/(10^15)));   

    num_ensembles = 100;                % 集合大小（50或100）
    for i = 1:num_ensembles
        % 生成标准正态分布随机场（均值为0，标准差为1）
        random_field = randn(180, 360);
        % 扰动公式：GPP_ensemble = GPP + GPP_uncertainty * random_field
        GPP_ensemble(:, :, i) = Mean_GPP + Mean_GPP_uncertainty .* random_field;
    end

    for ii=1:100
        GPP_temp=GPP_ensemble(:, :, ii);
        GPP_list_temp(ii) = nansum(nansum(GPP_temp.*area_grid/(10^15)));
    end
    GPP_std_altra(year-2014)=std(GPP_list_temp); clear GPP_ensemble GPP_list_temp

    % TER list
    Mean_TER_list(year-2014)=nansum(nansum(Mean_TER.*area_grid/(10^15)));

end


landflux_std=importdata("E:\phd_file\Tropical_2024\GCAS_2015-2024\Regional_uncertainty\annual\Global_land.txt");
landflux_std=landflux_std.data;
landflux_std=landflux_std(2:11,2:5);
landflux_std=mean(landflux_std,2);
landflux_std=landflux_std';

Fire_std=Mean_fire_list*0.2;
for i=1:length(GOSIF_GPP_list)
    GPP_std(i)=std(bootstrp(1000,@mean,[GOSIF_GPP_list(i),FluxSat_GPP_list(i)]));
end
GPP_std=sqrt(power(GPP_std,2)+power(GPP_std_altra,2));

TER_std=sqrt(power(Fire_std,2)+power(landflux_std,2)+power(GPP_std,2));

Global_data=importdata("E:\phd_file\Tropical_2024\CGR\CO2增长率_4Models.xlsx");
Global_data=Global_data.data;
MBL_data=Global_data(1:10,2);
MBL_data=MBL_data';

%%
f=figure
t = tiledlayout(2,2);
t.TileSpacing = 'tight';
t.Padding = 'compact';
set(gcf,'unit','pixels','position',[1000,774,807,464]);

nexttile
cn=bar(MBL_data,'Linestyle','none', 'FaceColor', 'flat')
for i=1:length(cn.CData)-1
    cn.CData(i,:) = [163,207,242]/255;
end
cn.CData(end,:) = [248 189 71]/255;
ylim([1,4])
set(gca,'YTick', [1:1:4]);
xlim([0.2,10.8])
set(gca,'XTick',[1:1:10],'FontName','Arial','fontsize',12)
set(gca,'XTickLabel',[],'FontName','Arial','fontsize',12)
ylabel('CGR (ppm yr^{-1})','FontName','Arial','FontSize',12);
text('string','a','Units','normalized','position',[-0.193372603257099 1.03434643966688 0],'FontName','Arial','FontSize',18,'fontweight','bold')

x=1:length(Mean_landflux_list);
nexttile
p1=plot(Mean_landflux_list,'-','LineWidth',3,'MarkerSize',15,'color','k'); hold on
fill([x, fliplr(x)], [Mean_landflux_list, fliplr(Mean_landflux_list+landflux_std)],'k','linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [Mean_landflux_list, fliplr(Mean_landflux_list-landflux_std)],'k','linestyle', 'none', 'FaceAlpha',0.2);
ylim([-3.3,0.4])
set(gca,'YTick', [-3:1:0]);
xlim([0.2,10.8])
set(gca,'XTick',[1:1:10],'FontName','Arial','fontsize',12)
set(gca,'XTickLabel',[],'FontName','Arial','fontsize',12)
ylabel('NBE (PgC yr^{-1})','FontName','Arial','FontSize',12);
text('string','b','Units','normalized','position',[-0.124459036984632 1.03434643966688 0],'FontName','Arial','FontSize',18,'fontweight','bold')

nexttile
p2=plot(Mean_GPP_list,'-','LineWidth',3,'MarkerSize',15,'color',[17,119,51]/255); hold on
fill([x, fliplr(x)], [Mean_GPP_list, fliplr(Mean_GPP_list+GPP_std)],[17,119,51]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [Mean_GPP_list, fliplr(Mean_GPP_list-GPP_std)],[17,119,51]/255,'linestyle', 'none', 'FaceAlpha',0.2);

p3=plot(Mean_TER_list,'-','LineWidth',3,'MarkerSize',15,'color',[255,165,0]/255); hold on
fill([x, fliplr(x)], [Mean_TER_list, fliplr(Mean_TER_list+TER_std)],[255,165,0]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [Mean_TER_list, fliplr(Mean_TER_list-TER_std)],[255,165,0]/255,'linestyle', 'none', 'FaceAlpha',0.2);

ylim([119,135])
set(gca,'YTick', [120:5:135]);
xlim([0.2,10.8])
set(gca,'XTick',[1:1:10],'FontName','Arial','fontsize',12)
set(gca,'XTickLabel',[2015:1:2024],'FontName','Arial','fontsize',12)
ylabel('GPP and TER (PgC yr^{-1})','FontName','Arial','FontSize',12);
legend([p2 p3],{'GPP','TER'},'NumColumns',1,'FontName','Arial','FontSize',12,'Box','off','Location','northwest')
text('string','c','Units','normalized','position',[-0.193372603257099 1.03434643966688 0],'FontName','Arial','FontSize',18,'fontweight','bold')

nexttile
p4=plot(Mean_fire_list,'-','LineWidth',3,'MarkerSize',15,'color','r'); hold on
fill([x, fliplr(x)], [Mean_fire_list, fliplr(Mean_fire_list+Fire_std)],'r','linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [Mean_fire_list, fliplr(Mean_fire_list-Fire_std)],'r','linestyle', 'none', 'FaceAlpha',0.2);
ylim([1,3])
% set(gca,'YTick', [-0.1:0.1:0.2]);
xlim([0.2,10.8])
set(gca,'XTick',[1:1:10],'FontName','Arial','fontsize',12)
set(gca,'XTickLabel',[2015:1:2024],'FontName','Arial','fontsize',12)
ylabel('Fire (PgC yr^{-1})','FontName','Arial','FontSize',12);
text('string','d','Units','normalized','position',[-0.124459036984632 1.03434643966688 0],'FontName','Arial','FontSize',18,'fontweight','bold')

result=['E:\phd_file\Tropical_2024\Result\V4\Carbon_budget_line.png']
% print(result,f,'-r600','-dpng');
