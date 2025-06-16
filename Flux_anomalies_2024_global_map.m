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
forest_mask(forest_mask==1 | forest_mask==2 |forest_mask==3 |forest_mask==4 |forest_mask==5| forest_mask==8)=666;
forest_mask(forest_mask~=666)=nan;
forest_mask(forest_mask==666)=1;

Non_forest_mask=LUCC_mask;
Non_forest_mask(Non_forest_mask==10 | Non_forest_mask==6 | Non_forest_mask==7 | Non_forest_mask==9)=666;
Non_forest_mask(Non_forest_mask~=666)=nan;
Non_forest_mask(Non_forest_mask==666)=1;

Cropland_mask=LUCC_mask;
Cropland_mask(Cropland_mask==12 | Cropland_mask==14)=666;
Cropland_mask(Cropland_mask~=666)=nan;
Cropland_mask(Cropland_mask==666)=1;


LUCC_mask(forest_mask==1 | Non_forest_mask==1| Cropland_mask==1)=1;
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
regions = {'forest', 'Non_forest', 'Cropland'};

for year=2015:2024

    % NBE
    % BEPS_GFAS_landflux=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\BEPS_GFAS\annual\NBE_BEPS_GFAS_' num2str(year) '.tif']);
    % BEPS_GFED_landflux=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\BEPS_GFED\annual\NBE_BEPS_GFED_' num2str(year) '.tif']);
    % CASA_GFAS_landflux=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\CASA_GFAS\annual\NBE_CASA_GFAS_' num2str(year) '.tif']);
    % CASA_GFED_landflux=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\CASA_GFED\annual\NBE_CASA_GFED_' num2str(year) '.tif']);
    Mean_NBE=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\mean_value\annual\NBE_' num2str(year) '.tif']);
    % All fire
    Mean_total_fire=importdata(['E:\phd_file\Tropical_2024\fire emission\mean_all_carbon_value\annual\Fire_' num2str(year) '.tif']);
    % NEE
    Mean_NEE=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\NEE\mean_value\annual\NEE_' num2str(year) '.tif']);
    % GPP
    GOSIF_GPP=importdata(['E:\phd_file\Tropical_2024\GPP\GOSIF\yearly\1degree\GOSIF_GPP_' num2str(year) '.tif']);
    FluxSat_GPP=importdata(['E:\phd_file\Tropical_2024\GPP\FluxSat\yearly\FluxSat_GPP_' num2str(year) '.tif']);
    Mean_GPP=importdata(['E:\phd_file\Tropical_2024\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    % TER
    Mean_TER=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\TER\annual\TER_' num2str(year) '.tif']);

    for i=1:length(regions)

        % Fire list
        region_name = regions{i};
        mask_var = [region_name '_mask'];

        Fire_var = [region_name '_Mean_fire_list'];
        % landflux_var = [region_name '_landflux_list'];
        Mean_landflux_var = [region_name '_Mean_landflux_list'];
        GOSIF_GPP_var = [region_name '_GOSIF_GPP_list'];
        FluxSat_GPP_var = [region_name '_FluxSat_GPP_list'];
        Mean_GPP_var = [region_name '_Mean_GPP_list'];
        Mean_TER_var = [region_name '_Mean_TER_list'];

        eval([Fire_var '(year-2014) = nansum(nansum(Mean_total_fire .* area_grid .* ' mask_var '/(10^15)));']);
        eval([Mean_landflux_var '(year-2014) = nansum(nansum(Mean_NBE .* area_grid .* ' mask_var '/(10^15)));']);
        % eval([landflux_var '(1,year-2014) = nansum(nansum(BEPS_GFAS_landflux .* area_grid .* ' mask_var '/(10^15)));']);
        % eval([landflux_var '(2,year-2014) = nansum(nansum(BEPS_GFED_landflux .* area_grid .* ' mask_var '/(10^15)));']);
        % eval([landflux_var '(3,year-2014) = nansum(nansum(CASA_GFAS_landflux .* area_grid .* ' mask_var '/(10^15)));']);
        % eval([landflux_var '(4,year-2014) = nansum(nansum(CASA_GFED_landflux .* area_grid .* ' mask_var '/(10^15)));']);
        eval([GOSIF_GPP_var '(year-2014) = nansum(nansum(GOSIF_GPP .* area_grid .* ' mask_var '/(10^15)));']);
        eval([FluxSat_GPP_var '(year-2014) = nansum(nansum(FluxSat_GPP .* area_grid .*' mask_var '/(10^15)));']);
        eval([Mean_GPP_var '(year-2014) = nansum(nansum(Mean_GPP .* area_grid .* ' mask_var '/(10^15)));']);
        eval([Mean_TER_var '(year-2014) = nansum(nansum(Mean_TER .* area_grid .* ' mask_var '/(10^15)));']);


    end
end

for i=1:length(regions)

    region_name = regions{i};

        eval(['landflux_temp = importdata(''E:\phd_file\Tropical_2024\GCAS_2015-2024\Regional_uncertainty\annual\' region_name '.txt'');']);
        eval('landflux_temp = landflux_temp.data;');
        eval('landflux_temp = landflux_temp(2:11,2:5);');
        eval('landflux_temp = mean(landflux_temp,2);');
        eval('landflux_temp = landflux_temp'';');
        
        landflux_var = [region_name '_landflux_2024anomaly_std'];
        eval([landflux_var '= sqrt(power(landflux_temp(end),2)+power(landflux_temp(end-2),2));']);
        clear landflux_temp

end

for i=1:length(regions)

    region_name = regions{i};
    mask_var = [region_name '_mask'];
    Mean_GPP_var = [region_name '_GPP_2024anomaly_std_altra'];

    GPP2024=importdata("E:\phd_file\Tropical_2024\GPP\mean_value\year\GPP_2024.tif");
    GPP2024_uncertainty=importdata("E:\phd_file\Tropical_2024\GPP\mean_value\uncertainty\yearly\GPP_uncertainty_2024.tif");

    GPP2022=importdata("E:\phd_file\Tropical_2024\GPP\mean_value\year\GPP_2022.tif");
    GPP2022_uncertainty=importdata("E:\phd_file\Tropical_2024\GPP\mean_value\uncertainty\yearly\GPP_uncertainty_2022.tif");

    num_ensembles = 100;                % 集合大小（50或100）

    for cc = 1:num_ensembles
        % 生成标准正态分布随机场（均值为0，标准差为1）
        random_field1 = randn(180, 360);
        random_field2 = randn(180, 360);
        % 扰动公式：GPP_ensemble = GPP + GPP_uncertainty * random_field
        GPP_2024ensemble(:, :, cc) = GPP2024 + GPP2024_uncertainty .* random_field1;
        GPP_2022ensemble(:, :, cc) = GPP2022 + GPP2022_uncertainty .* random_field2;

    end


    for ii=1:100
        GPP_2024temp=GPP_2024ensemble(:, :, ii);
        GPP_2022temp=GPP_2022ensemble(:, :, ii);
        eval(['GPP_2024list_temp(ii) = nansum(nansum(GPP_2024temp.*area_grid.*' mask_var '/(10^15)));']);
        eval(['GPP_2022list_temp(ii) = nansum(nansum(GPP_2022temp.*area_grid.*' mask_var '/(10^15)));']);
    end
    eval([Mean_GPP_var '= sqrt(power(std(GPP_2024list_temp),2)+power(std(GPP_2022list_temp),2));' ]);

    clear GPP_2024temp  GPP_2022temp GPP_2024list_temp GPP_2022list_temp

end

% 2024 Anomaly relative to 2022
for i=1:length(regions)

    % Fire list
    region_name = regions{i};

    Fire_var = [region_name '_Fire_2024anomaly'];
    landflux_var = [region_name '_landflux_2024anomaly'];
    Mean_GPP_var = [region_name '_GPP_2024anomaly'];
    Mean_TER_var = [region_name '_TER_2024anomaly'];

    eval([Fire_var '=' region_name '_Mean_fire_list(:,end)-' region_name '_Mean_fire_list(:,end-2);']);
    eval([landflux_var '=' region_name '_Mean_landflux_list(:,end)-' region_name '_Mean_landflux_list(:,end-2);']);
    eval([Mean_GPP_var '=' region_name '_Mean_GPP_list(:,end)-' region_name '_Mean_GPP_list(:,end-2);']);
    eval([Mean_TER_var '=' region_name '_Mean_TER_list(:,end)-' region_name '_Mean_TER_list(:,end-2);']);


end
% STD
for i=1:length(regions)

    % Fire list
    region_name = regions{i};

    Fire_var = [region_name '_Fire_2024anomaly_std'];
    % landflux_var = [region_name '_landflux_2024anomaly_std'];
    Mean_GPP_var = [region_name '_GPP_2024anomaly_std'];
    Mean_GPP_var_altra = [region_name '_GPP_2024anomaly_std_altra'];
    Mean_TER_var = [region_name '_TER_2024anomaly_std'];

    % eval([landflux_var '=std([' region_name '_landflux_list(1,end)-' region_name '_landflux_list(1,end-2),' ...
    %     region_name '_landflux_list(2,end)-' region_name '_landflux_list(2,end-2),' region_name '_landflux_list(3,end)-' region_name '_landflux_list(3,end-2),' ...
    %     region_name '_landflux_list(4,end)-' region_name '_landflux_list(4,end-2)]);']);

    eval([Fire_var '=' region_name '_Fire_2024anomaly*0.2;']);
    eval([Mean_GPP_var '=std(bootstrp(1000,@mean,[' region_name '_GOSIF_GPP_list(end)-' region_name ...
        '_GOSIF_GPP_list(end-2),' region_name '_FluxSat_GPP_list(end)-' region_name '_FluxSat_GPP_list(end-2)]));']);
    eval([Mean_GPP_var '=sqrt(power(' Mean_GPP_var ',2)+power(' Mean_GPP_var_altra ',2));']);
    eval([Mean_TER_var '=sqrt(power(' region_name '_Fire_2024anomaly_std,2)+power(' region_name '_landflux_2024anomaly_std,2)+power(' region_name '_GPP_2024anomaly_std,2));']);

end

NBE_result=[forest_landflux_2024anomaly,Non_forest_landflux_2024anomaly,Cropland_landflux_2024anomaly];
NBE_result_std=[forest_landflux_2024anomaly_std,Non_forest_landflux_2024anomaly_std,Cropland_landflux_2024anomaly_std];

GPP_result=[forest_GPP_2024anomaly,Non_forest_GPP_2024anomaly,Cropland_GPP_2024anomaly];
GPP_result_std=[forest_GPP_2024anomaly_std,Non_forest_GPP_2024anomaly_std,Cropland_GPP_2024anomaly_std];

TER_result=[forest_TER_2024anomaly,Non_forest_TER_2024anomaly,Cropland_TER_2024anomaly];
TER_result_std=[forest_TER_2024anomaly_std,Non_forest_TER_2024anomaly_std,Cropland_TER_2024anomaly_std];

Fire_result=[forest_Fire_2024anomaly,Non_forest_Fire_2024anomaly,Cropland_Fire_2024anomaly];
Fire_result_std=[forest_Fire_2024anomaly_std,Non_forest_Fire_2024anomaly_std,Cropland_Fire_2024anomaly_std];
%% 2024 anomaly

NBE_year2024=importdata("E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\mean_value\annual\NBE_2024.tif").*pixel_mask;
NEE_year2024=importdata("E:\phd_file\Tropical_2024\GCAS_2015-2024\NEE\mean_value\annual\NEE_2024.tif").*pixel_mask;
GPP_year2024=importdata("E:\phd_file\Tropical_2024\GPP\mean_value\year\GPP_2024.tif").*pixel_mask;
TER_year2024=importdata("E:\phd_file\Tropical_2024\GCAS_2015-2024\TER\annual\TER_2024.tif").*pixel_mask;
Fire_year2024=importdata("E:\phd_file\Tropical_2024\fire emission\mean_all_carbon_value\annual\Fire_2024.tif").*pixel_mask;

NBE_year2022=importdata("E:\phd_file\Tropical_2024\GCAS_2015-2024\NBE\mean_value\annual\NBE_2022.tif").*pixel_mask;
NEE_year2022=importdata("E:\phd_file\Tropical_2024\GCAS_2015-2024\NEE\mean_value\annual\NEE_2022.tif").*pixel_mask;
GPP_year2022=importdata("E:\phd_file\Tropical_2024\GPP\mean_value\year\GPP_2022.tif").*pixel_mask;
TER_year2022=importdata("E:\phd_file\Tropical_2024\GCAS_2015-2024\TER\annual\TER_2022.tif").*pixel_mask;
Fire_year2022=importdata("E:\phd_file\Tropical_2024\fire emission\mean_all_carbon_value\annual\Fire_2022.tif").*pixel_mask;

% 2024 calculate
NBE_Anomaly_2024=NBE_year2024-NBE_year2022;
NEE_Anomaly_2024=NEE_year2024-NEE_year2022;
GPP_Anomaly_2024=GPP_year2024-GPP_year2022;
TER_Anomaly_2024=TER_year2024-TER_year2022;
Fire_Anomaly_2024=Fire_year2024-Fire_year2022;

%%
% 面积占比
count_sum=sum(sum(~isnan(NBE_Anomaly_2024)));
count1=sum(sum(NBE_Anomaly_2024>0));
result=count1/(count_sum)

%
forest_landflux_2024anomaly/(forest_landflux_2024anomaly+Non_forest_landflux_2024anomaly+Cropland_landflux_2024anomaly)
Non_forest_landflux_2024anomaly/(forest_landflux_2024anomaly+Non_forest_landflux_2024anomaly+Cropland_landflux_2024anomaly)
Cropland_landflux_2024anomaly/(forest_landflux_2024anomaly+Non_forest_landflux_2024anomaly+Cropland_landflux_2024anomaly)



%%

ff=figure
t = tiledlayout(2,2);
set(gcf,'unit','pixels','position',[691,69,1051,662]);

t.TileSpacing = 'tight';
t.Padding = 'tight';
t.OuterPosition = [0 0.025 1 0.975]; % 为图例留出底部空间

a=nexttile
% part 1 *********************************************************************************
m_proj('robinson','long',[-180 180],'lat',[-90 90]); hold on
m_pcolor(lon1,lat1,NBE_Anomaly_2024,'linestyle','none');
m_coast('linewidth',1,'color','k'); %填充矢量边界
% m_line(bou_canX,bou_canY,'linewidth',1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[],'ytick',[],...
    'linewidth',1,'xaxisloc','bottom','yaxisloc','left','FontSize',12,'xticklabels',[],'yticklabels',[],'linestyle','none'); %对齐格网
colormap(a,flipud(nclCM(150,12)));
caxis([-100,100]);
t1=caxis;
t1=linspace(t1(1),t1(2),5);
h=colorbar('location','southoutside','FontName','Arial','FontSize',12,'ytick',t1);
set(get(h,'title'),'string','NBE anomalies (gC m^{-2} yr^{-1})','fontsize',12);
% title('GOSIF GPP','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[0.00502655592092706 0.953054389856126 0],'FontName','Arial','FontSize',18,'fontweight','bold')
% text('string','NBP','Units','normalized','position',[0.905202275694663 0.840959082692752 0],'FontName','Arial','FontSize',14,'fontweight','bold')

c=nexttile
% part 1 *********************************************************************************
m_proj('robinson','long',[-180 180],'lat',[-90 90]); hold on
m_pcolor(lon1,lat1,GPP_Anomaly_2024,'linestyle','none');
m_coast('linewidth',1,'color','k'); %填充矢量边界
% m_line(bou_canX,bou_canY,'linewidth',1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[],'ytick',[],...
    'linewidth',1,'xaxisloc','bottom','yaxisloc','left','FontSize',12,'xticklabels',[],'yticklabels',[],'linestyle','none'); %对齐格网
colormap(c,nclCM(399,12));
caxis([-200,200]);
t1=caxis;
t1=linspace(t1(1),t1(2),5);
h=colorbar('location','southoutside','FontName','Arial','FontSize',12,'ytick',t1);
set(get(h,'title'),'string','GPP anomalies (gC m^{-2} yr^{-1})','fontsize',12);
% title('GOSIF GPP','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','b','Units','normalized','position',[0.00502655592092706 0.953054389856126 0],'FontName','Arial','FontSize',18,'fontweight','bold')
% text('string','GPP','Units','normalized','position',[0.905202275694663 0.840959082692752 0],'FontName','Arial','FontSize',14,'fontweight','bold')

d=nexttile
% part 1 *********************************************************************************
m_proj('robinson','long',[-180 180],'lat',[-90 90]); hold on
m_pcolor(lon1,lat1,TER_Anomaly_2024,'linestyle','none');
m_coast('linewidth',1,'color','k'); %填充矢量边界
% m_line(bou_canX,bou_canY,'linewidth',1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[],'ytick',[],...
    'linewidth',1,'xaxisloc','bottom','yaxisloc','left','FontSize',12,'xticklabels',[],'yticklabels',[],'linestyle','none'); %对齐格网
colormap(d,flipud(nclCM(399,12)));
caxis([-200,200]);
t1=caxis;
t1=linspace(t1(1),t1(2),5);
h=colorbar('location','southoutside','FontName','Arial','FontSize',12,'ytick',t1);
set(get(h,'title'),'string','TER anomalies (gC m^{-2} yr^{-1})','fontsize',12);
% title('GOSIF GPP','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[0.00502655592092706 0.953054389856126 0],'FontName','Arial','FontSize',18,'fontweight','bold')
% text('string','TER','Units','normalized','position',[0.905202275694663 0.840959082692752 0],'FontName','Arial','FontSize',14,'fontweight','bold')

e=nexttile
% part 1 *********************************************************************************
m_proj('robinson','long',[-180 180],'lat',[-90 90]); hold on
m_pcolor(lon1,lat1,Fire_Anomaly_2024,'linestyle','none');
m_coast('linewidth',1,'color','k'); %填充矢量边界
% m_line(bou_canX,bou_canY,'linewidth',1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[],'ytick',[],...
    'linewidth',1,'xaxisloc','bottom','yaxisloc','left','FontSize',12,'xticklabels',[],'yticklabels',[],'linestyle','none'); %对齐格网
colormap(e,nclCM(326,12));
caxis([-40,40]);
t1=caxis;
t1=linspace(t1(1),t1(2),5);
h=colorbar('location','southoutside','FontName','Arial','FontSize',12,'ytick',t1);
set(get(h,'title'),'string','Fire anomalies (gC m^{-2} yr^{-1})','fontsize',12);
% title('GOSIF GPP','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','d','Units','normalized','position',[0.00502655592092706 0.953054389856126 0],'FontName','Arial','FontSize',18,'fontweight','bold')
% text('string','Fire','Units','normalized','position',[0.905202275694663 0.840959082692752 0],'FontName','Arial','FontSize',14,'fontweight','bold')


a1=axes('Position',[0.0944880315256473 0.692954682779458,0.063831234304131,0.175204485021809]);
c1 = bar(1,NBE_result(1),'Facecolor',[176,187,208]/255,'LineStyle','none'); hold on
c2 = bar(2,NBE_result(2),'Facecolor',[126,223,223]/255,'LineStyle','none'); hold on
c3 = bar(3,NBE_result(3),'Facecolor',[247,223,165]/255,'LineStyle','none'); hold on
errorbar(1:3,NBE_result,NBE_result_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',5,'Marker','none')
ylim([0,1.7])
set(gca,'yTick',[0:0.6:1.2],'FontName','Arial','fontsize',8)
set(gca,'xTick',[],'FontName','Arial','fontsize',8)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',8)
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',8)
box off
set(gca, 'Color', 'none')
ax = gca;
ax.XAxis.Visible = 'off'; % 隐藏Y轴线
lgd=legend([c1,c2,c3],{'Forest','Shrub/grassland','Cropland'},'NumColumns',3,'FontName','Arial','FontSize',12,'Box','off','Position',[0.341606273323407,0.000294309518182,0.351152067470111,0.045999998867512]);
% text('string','PgC yr^{-1}','Units','normalized','position',[-0.218758035014238,1.114628110253375,0],'FontName','Arial','FontSize',8)
lgd.ItemTokenSize = [30 8];

a2=axes('Position',[0.590950474579665 0.692827794561934,0.063831234304131,0.175204485021809]);
c1 = bar(1,GPP_result(1),'Facecolor',[176,187,208]/255,'LineStyle','none'); hold on
c2 = bar(2,GPP_result(2),'Facecolor',[126,223,223]/255,'LineStyle','none'); hold on
c3 = bar(3,GPP_result(3),'Facecolor',[247,223,165]/255,'LineStyle','none'); hold on
errorbar(1:3,GPP_result,GPP_result_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',5,'Marker','none')
ylim([-0.4,1.4])
set(gca,'yTick',[0:0.5:1],'FontName','Arial','fontsize',8)
set(gca,'xTick',[],'FontName','Arial','fontsize',8)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',8)
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',8)
box off
set(gca, 'Color', 'none')
ax = gca;
ax.XAxis.Visible = 'off'; % 隐藏Y轴线
% text('string','PgC yr^{-1}','Units','normalized','position',[-0.218758035014238,1.114628110253375,0],'FontName','Arial','FontSize',8)

a3=axes('Position',[0.0955587695105384 0.203444108761329,0.063831234304131,0.175204485021809]);
c1 = bar(1,TER_result(1),'Facecolor',[176,187,208]/255,'LineStyle','none'); hold on
c2 = bar(2,TER_result(2),'Facecolor',[126,223,223]/255,'LineStyle','none'); hold on
c3 = bar(3,TER_result(3),'Facecolor',[247,223,165]/255,'LineStyle','none'); hold on
errorbar(1:3,TER_result,TER_result_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',5,'Marker','none')
% ylim([0,1.5])
set(gca,'yTick',[0:1:2],'FontName','Arial','fontsize',8)
set(gca,'xTick',[],'FontName','Arial','fontsize',8)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',8)
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',8)
box off
set(gca, 'Color', 'none')
ax = gca;
ax.XAxis.Visible = 'off'; % 隐藏Y轴线
% text('string','PgC yr^{-1}','Units','normalized','position',[-0.218758035014238,1.114628110253375,0],'FontName','Arial','FontSize',8)

a4=axes('Position',[0.594517847325391 0.206465256797583 0.063831234304131,0.175204485021809]);
c1 = bar(1,Fire_result(1),'Facecolor',[176,187,208]/255,'LineStyle','none'); hold on
c2 = bar(2,Fire_result(2),'Facecolor',[126,223,223]/255,'LineStyle','none'); hold on
c3 = bar(3,Fire_result(3),'Facecolor',[247,223,165]/255,'LineStyle','none'); hold on
errorbar(1:3,Fire_result,Fire_result_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',5,'Marker','none')
ylim([0,0.55])
% set(gca,'yTick',[0:0.6:1.2],'FontName','Arial','fontsize',8)
set(gca,'xTick',[],'FontName','Arial','fontsize',8)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',8)
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',8)
box off
set(gca, 'Color', 'none')
ax = gca;
ax.XAxis.Visible = 'off'; % 隐藏Y轴线
% text('string','PgC yr^{-1}','Units','normalized','position',[-0.218758035014238 1.114628110253375 0],'FontName','Arial','FontSize',8)

result=['E:\phd_file\Tropical_2024\Result\V4\Flux_anomalies_2024_global_map.png']
% print(result,ff,'-r600','-dpng');

