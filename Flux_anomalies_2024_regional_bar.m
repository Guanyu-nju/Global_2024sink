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
%% 2024 anomaly

regions = {'North_America', 'South_America', 'Russia', 'Europe', 'WS_Asia', 'Africa', 'East_Asia', 'SE_Asia', 'Oceania', 'Tropical', 'Boreal', 'North_temperate', 'South_temperate'};

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

        Fire_var = [region_name '_Mean_fire_list'];
        % landflux_var = [region_name '_landflux_list'];
        Mean_landflux_var = [region_name '_Mean_landflux_list'];
        GOSIF_GPP_var = [region_name '_GOSIF_GPP_list'];
        FluxSat_GPP_var = [region_name '_FluxSat_GPP_list'];
        Mean_GPP_var = [region_name '_Mean_GPP_list'];
        Mean_TER_var = [region_name '_Mean_TER_list'];

        mask_var = [region_name '_mask'];

        eval([Fire_var '(year-2014) = nansum(nansum(Mean_total_fire .* area_grid .* ' mask_var '/(10^15)));']);
        eval([Mean_landflux_var '(year-2014) = nansum(nansum(Mean_NBE .* area_grid .* ' mask_var '/(10^15)));']);
        % eval([landflux_var '(1,year-2014) = nansum(nansum(BEPS_GFAS_landflux .* area_grid .* ' mask_var '/(10^15)));']);
        % eval([landflux_var '(2,year-2014) = nansum(nansum(BEPS_GFED_landflux .* area_grid .* ' mask_var '/(10^15)));']);
        % eval([landflux_var '(3,year-2014) = nansum(nansum(CASA_GFAS_landflux .* area_grid .* ' mask_var '/(10^15)));']);
        % eval([landflux_var '(4,year-2014) = nansum(nansum(CASA_GFED_landflux .* area_grid .* ' mask_var '/(10^15)));']);
        eval([GOSIF_GPP_var '(year-2014) = nansum(nansum(GOSIF_GPP .* area_grid .* ' mask_var '/(10^15)));']);
        eval([FluxSat_GPP_var '(year-2014) = nansum(nansum(FluxSat_GPP .* area_grid .* ' mask_var '/(10^15)));']);
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
% Carbon budget result
for i=1:length(regions)

    % Fire list
    region_name = regions{i};

    Carbon_budget_var = [region_name '_Carbon_budegt_2024anomaly'];
    Carbon_budget_std_var = [region_name '_Carbon_budegt_2024anomaly_std'];

    eval([Carbon_budget_var '=[' region_name '_landflux_2024anomaly,' region_name '_GPP_2024anomaly,' region_name '_TER_2024anomaly,' region_name '_Fire_2024anomaly];']);
    eval([Carbon_budget_std_var '=[' region_name '_landflux_2024anomaly_std,' region_name '_GPP_2024anomaly_std,' region_name '_TER_2024anomaly_std,' region_name '_Fire_2024anomaly_std];']);
end

Latitudinal_Carbon_budget_2024anomaly=[Boreal_Carbon_budegt_2024anomaly;North_temperate_Carbon_budegt_2024anomaly;...
    Tropical_Carbon_budegt_2024anomaly;South_temperate_Carbon_budegt_2024anomaly];

Latitudinal_Carbon_budget_2024anomaly_std=[Boreal_Carbon_budegt_2024anomaly_std;North_temperate_Carbon_budegt_2024anomaly_std;...
    Tropical_Carbon_budegt_2024anomaly_std;South_temperate_Carbon_budegt_2024anomaly_std];
%%
ff=figure
set(gcf,'unit','pixels','position',[634,176,1446,959]);


a=axes('Position',[0.300506911907889,0.304483837330553,0.422176352269151,0.420229405630865]);
% a=axes('Position',[0.266251728907331 0.274244004171011 0.470262793914246 0.464025026068822]);

% part 1 *********************************************************************************
m_proj('robinson','long',[-180 180],'lat',[-90 90]); hold on
m_pcolor(lon1,lat1,Global_background,'linestyle','none');
m_coast('linewidth',1,'color','k'); %填充矢量边界
% m_line(bou_canX,bou_canY,'linewidth',1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[],'ytick',[],...
    'linewidth',1,'xaxisloc','bottom','yaxisloc','left','FontSize',12,'xticklabels',[],'yticklabels',[],'linestyle','none'); %对齐格网
mycolor=[151,202,61; 214,230,167; 174,198,232; 209,235,250; 246,140,64; 204,152,102; 251,186,120; 244,211,97; 116,168,218]/255;
colormap(mycolor);


customColors = [
    0.17 0.42 0.7;    % 绿色
    0.0667 0.4667 0.2;    % 蓝色
    1 0.6471 0;    % 黄色
    0.698 0.1333 0.1333 % 紫色
    ];

b=axes('Position',[0.044790727736028,0.511991657977059,0.214545371849035,0.22450915816199]);
cn = bar(North_America_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,North_America_Carbon_budegt_2024anomaly,North_America_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-1,1])
set(gca,'yTick',[-1:0.5:1],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','North America','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','d','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

c=axes('Position',[0.044790727736028,0.264590385746401,0.214545371849035,0.22450915816199]);
cn = bar(South_America_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,South_America_Carbon_budegt_2024anomaly,South_America_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-1,1])
set(gca,'yTick',[-1:0.5:1],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','South America','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','f','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

d=axes('Position',[0.259866799658575,0.761209593326381,0.214545371849035,0.22450915816199]);
cn = bar(Europe_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,Europe_Carbon_budegt_2024anomaly,Europe_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-0.63,0.63])
set(gca,'yTick',[-0.6:0.3:0.6],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','Europe','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','a','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

e=axes('Position',[0.520724337694536,0.761209593326381,0.214545371849035,0.22450915816199]);
cn = bar(WS_Asia_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,WS_Asia_Carbon_budegt_2024anomaly,WS_Asia_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-1,1])
set(gca,'yTick',[-1:0.5:1],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','West/Central/South Asia','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','b','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

f=axes('Position',[0.781581875730498,0.761209593326381,0.214545371849035,0.22450915816199]);
cn = bar(Russia_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,Russia_Carbon_budegt_2024anomaly,Russia_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-0.6,0.6])
set(gca,'yTick',[-0.6:0.3:0.6],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','Russia','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','c','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

g=axes('Position',[0.781581875730498,0.511991657977059,0.214545371849035,0.22450915816199]);
cn = bar(East_Asia_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,East_Asia_Carbon_budegt_2024anomaly,East_Asia_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-0.68,0.68])
set(gca,'yTick',[-0.6:0.3:0.6],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','East Asia','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','e','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

h=axes('Position',[0.781581875730498,0.264590385746401,0.214545371849035,0.22450915816199]);
cn = bar(SE_Asia_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,SE_Asia_Carbon_budegt_2024anomaly,SE_Asia_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-0.45,0.45])
set(gca,'yTick',[-0.4:0.2:0.4],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','Southeast Asia','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','g','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

i=axes('Position',[0.781581875730498,0.017666506705733,0.214545371849035,0.22450915816199]);
cn = bar(Oceania_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,Oceania_Carbon_budegt_2024anomaly,Oceania_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-0.26,0.26])
set(gca,'yTick',[-0.2:0.1:0.2],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','Oceania','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','j','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

j=axes('Position',[0.520724337694536,0.017666506705733,0.214545371849035,0.22450915816199]);
cn = bar(Africa_Carbon_budegt_2024anomaly,'LineStyle','none'); hold on
set(cn, 'FaceColor', 'flat'); % 设置为单独配色模式
set(cn, 'CData', customColors); % 应用颜色矩阵
errorbar(1:4,Africa_Carbon_budegt_2024anomaly,Africa_Carbon_budegt_2024anomaly_std,'LineStyle','none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
ylim([-1.2,1.2])
set(gca,'yTick',[-1.2:0.6:1.2],'FontName','Arial','fontsize',12)
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','Africa','Units','normalized','position',[0.0376597117821248 0.0752726737914311 0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','i','Units','normalized','position',[0.0344339052183526 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

k=axes('Position',[0.044790727736028,0.017666506705733,0.419801250133955,0.22450915816199]);
bb = bar(Latitudinal_Carbon_budget_2024anomaly,'FaceColor','flat','LineStyle','none'); hold on
bb(1).CData = customColors(1,:);
bb(2).CData = customColors(2,:);
bb(3).CData = customColors(3,:);
bb(4).CData = customColors(4,:);
%分组误差棒
[M,N]=size(Latitudinal_Carbon_budget_2024anomaly);
for i=1:N
    xx(:,i)=bb(i).XEndPoints';
end
h2=errorbar(xx(:,:),Latitudinal_Carbon_budget_2024anomaly,Latitudinal_Carbon_budget_2024anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
pp=plot([xx(1,end)+(xx(2,1)-xx(1,end))/2,xx(1,end)+(xx(2,1)-xx(1,end))/2],[-1,5],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
pp=plot([xx(2,end)+(xx(2,1)-xx(1,end))/2,xx(2,end)+(xx(2,1)-xx(1,end))/2],[-1,5],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
pp=plot([xx(3,end)+(xx(2,1)-xx(1,end))/2,xx(3,end)+(xx(2,1)-xx(1,end))/2],[-1,5],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
xlim([0.5,4.5])
ylim([-0.7,2.2])
set(gca,'YTick', [-0.5:0.5:2],'FontName','Arial','FontSize',14);
set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',[],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
lgd=legend(bb,{'\DeltaNBE','\DeltaGPP','\DeltaTER','\DeltaFire'},'NumColumns',2,'FontName','Arial','FontSize',12,'Box','on','Position',[0.0450410544822098 0.840458811261731 0.140989374286808 0.0583941602887245])
text('string','50–90°N','Units','normalized','position',[0.064018855109967,0.052016858984925,0],'FontName','Arial','FontSize',14,'fontweight','bold')
text('string','23–50°N','Units','normalized','position',[0.307840930892504,0.052016858984925,0],'FontName','Arial','FontSize',14,'fontweight','bold')
text('string','23°N–23°S','Units','normalized','position',[0.551663006675041,0.052016858984925,0],'FontName','Arial','FontSize',14,'fontweight','bold')
text('string','23–60°S','Units','normalized','position',[0.813606993495469,0.052016858984925,0],'FontName','Arial','FontSize',14,'fontweight','bold')
% text('string','h','Units','normalized','position',[0.019606887096442 0.921784283350763 0],'FontName','Arial','FontSize',18,'fontweight','bold')

annotation('arrow',[0.398340248962656 0.260719225449516],...
    [0.598582898852972 0.635036496350365]);
annotation('arrow',[0.441908713692946 0.260719225449516],...
    [0.465110531803962 0.363920750782065]);
annotation('arrow',[0.529045643153527 0.36237897648686],...
    [0.609010427528676 0.758081334723671]);
annotation('arrow',[0.580221300138313 0.622406639004149],...
    [0.587112617309698 0.759124087591241]);
annotation('arrow',[0.632780082987552 0.741355463347165],...
    [0.577727841501564 0.614181438998957]);
annotation('arrow',[0.645228215767635 0.738589211618257],...
    [0.5130771637122 0.454640250260688]);
annotation('arrow',[0.671507607192255 0.760719225449516],...
    [0.451554744525547 0.227320125130344]);
annotation('arrow',[0.535961272475795 0.606500691562932],...
    [0.494307612095933 0.2429614181439]);
annotation('arrow',[0.612724757952974 0.778008298755187],...
    [0.639250260688217 0.793534932221064]);
result=['E:\phd_file\Tropical_2024\Result\V4\Flux_anomalies_2024_regional_bar.png']
% print(result,ff,'-r600','-dpng');
