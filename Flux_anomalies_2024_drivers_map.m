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

pixel_mask(pixel_mask<1 | pixel_mask>10)=nan;
pixel_mask(~isnan(pixel_mask))=1;

dominant_factor=importdata("E:\phd_file\Tropical_2024\domnant_factor.tif").*pixel_mask;

Positive_Fire_mask=dominant_factor;
Positive_Fire_mask(Positive_Fire_mask~=1)=nan;

Positive_TER_mask=dominant_factor;
Positive_TER_mask(Positive_TER_mask~=2)=nan;
Positive_TER_mask(~isnan(Positive_TER_mask))=1;

Negative_GPP_mask=dominant_factor;
Negative_GPP_mask(Negative_GPP_mask~=3)=nan;
Negative_GPP_mask(~isnan(Negative_GPP_mask))=1;

Positive_GPP_mask=dominant_factor;
Positive_GPP_mask(Positive_GPP_mask~=4)=nan;
Positive_GPP_mask(~isnan(Positive_GPP_mask))=1;

Negative_TER_mask=dominant_factor;
Negative_TER_mask(Negative_TER_mask~=5)=nan;
Negative_TER_mask(~isnan(Negative_TER_mask))=1;

Negative_Fire_mask=dominant_factor;
Negative_Fire_mask(Negative_Fire_mask~=6)=nan;
Negative_Fire_mask(~isnan(Negative_Fire_mask))=1;


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

area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;

%% 2024 anomaly

regions = {'Positive_Fire', 'Positive_TER', 'Negative_GPP', 'Positive_GPP', 'Negative_TER'};

for year=2000:2024

    TEM_temp=importdata(['E:\phd_file\Tropical_2024\Air_temperature\annual\TEM_' num2str(year) '.tif']);
    VPD_temp=importdata(['E:\phd_file\Tropical_2024\VPD\annual\VPD_' num2str(year) '.tif']);
    SR_temp=importdata(['E:\phd_file\Tropical_2024\Solar_Radiation\annual\SR_' num2str(year) '.tif']);

    for i=1:length(regions)

        region_name = regions{i};

        TEM_var = [region_name '_TEM_list'];
        VPD_var = [region_name '_VPD_list'];
        SR_var = [region_name '_SR_list'];

        mask_var = [region_name '_mask'];

        eval([TEM_var '(year-1999) = nansum(nansum(TEM_temp .* area_grid .* ' mask_var '))/nansum(nansum(area_grid .*' mask_var '));']);
        eval([VPD_var '(year-1999) = nansum(nansum(VPD_temp .* area_grid .* ' mask_var '))/nansum(nansum(area_grid .*' mask_var '));']);
        eval([SR_var '(year-1999) = nansum(nansum(SR_temp .* area_grid .* ' mask_var '))/nansum(nansum(area_grid .*' mask_var '));']);

    end
end

for year=2003:2024

    TWS_temp=importdata(['E:\phd_file\Tropical_2024\TWS\year\long_term\1degree\TWS_' num2str(year) '.tif']);

    for i=1:length(regions)

        region_name = regions{i};

        TWS_var = [region_name '_TWS_list'];
        mask_var = [region_name '_mask'];
        eval([TWS_var '(year-2002) = nansum(nansum(TWS_temp .* area_grid .* ' mask_var '))/nansum(nansum(area_grid .*' mask_var '));']);

    end

end
for i=1:length(regions)

    region_name = regions{i};
    TWS_var = [region_name '_TWS_list'];
    eval([TWS_var '(2017-2003+1:2018-2003+1)=[];']);

end

% 2024 Anomaly relative to 2022
for i=1:length(regions)

    % Fire list
    region_name = regions{i};

    TEM_var = [region_name '_TEM_anomalies'];
    VPD_var = [region_name '_VPD_anomalies'];
    SR_var = [region_name '_SR_anomalies'];
    TWS_var = [region_name '_TWS_anomalies'];

    eval([TEM_var '=(' region_name '_TEM_list(end)-' region_name '_TEM_list(end-2))/nanstd(' region_name '_TEM_list(1:end-3));']);
    eval([VPD_var '=(' region_name '_VPD_list(end)-' region_name '_VPD_list(end-2))/nanstd(' region_name '_VPD_list(1:end-3));']);
    eval([SR_var '=(' region_name '_SR_list(end)-' region_name '_SR_list(end-2))/nanstd(' region_name '_SR_list(1:end-3));']);
    eval([TWS_var '=(' region_name '_TWS_list(end)-' region_name '_TWS_list(end-2))/nanstd(' region_name '_TWS_list(1:end-3));']);

end
% Climate result
for i=1:length(regions)


    region_name = regions{i};

    Climate_var = [region_name '_Climate_2024anomaly'];
    eval([Climate_var '=[' region_name '_TEM_anomalies,' region_name '_TWS_anomalies,' region_name '_VPD_anomalies,' region_name '_SR_anomalies];']);
end

Driver_Climate_2024anomaly=[Positive_Fire_Climate_2024anomaly;Positive_TER_Climate_2024anomaly;...
    Negative_GPP_Climate_2024anomaly;Positive_GPP_Climate_2024anomaly;Negative_TER_Climate_2024anomaly];

Ratio_result=[sum(sum(~isnan(Positive_Fire_mask)))/sum(sum(~isnan(dominant_factor))),sum(sum(~isnan(Positive_TER_mask)))/sum(sum(~isnan(dominant_factor))),...
    sum(sum(~isnan(Negative_GPP_mask)))/sum(sum(~isnan(dominant_factor))),...
    sum(sum(~isnan(Positive_GPP_mask)))/sum(sum(~isnan(dominant_factor))),sum(sum(~isnan(Negative_TER_mask)))/sum(sum(~isnan(dominant_factor))),sum(sum(~isnan(Negative_Fire_mask)))/sum(sum(~isnan(dominant_factor)))]*100;
%%
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

        Tropical_landflux_var = ['Tropical_' region_name '_Mean_landflux_list'];
        Boreal_landflux_var = ['Boreal_' region_name '_Mean_landflux_list'];
        North_temperate_landflux_var = ['North_temperate_' region_name '_Mean_landflux_list'];
        South_temperate_landflux_var = ['South_temperate_' region_name '_Mean_landflux_list'];

        eval([Tropical_landflux_var '(year-2014) = nansum(nansum(Mean_NBE .* area_grid .* Tropical_mask .* ' mask_var '/(10^15)));']);
        eval([Boreal_landflux_var '(year-2014) = nansum(nansum(Mean_NBE .* area_grid .* Boreal_mask .* ' mask_var '/(10^15)));']);
        eval([North_temperate_landflux_var '(year-2014) = nansum(nansum(Mean_NBE .* area_grid .* North_temperate_mask .* ' mask_var '/(10^15)));']);
        eval([South_temperate_landflux_var '(year-2014) = nansum(nansum(Mean_NBE .* area_grid .* South_temperate_mask .* ' mask_var '/(10^15)));']);


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

    Tropical_landflux_var = ['Tropical_' region_name '_landflux_2024anomaly'];
    Boreal_landflux_var = ['Boreal_' region_name '_landflux_2024anomaly'];
    North_temperate_landflux_var = ['North_temperate_' region_name '_landflux_2024anomaly'];
    South_temperate_landflux_var = ['South_temperate_' region_name '_landflux_2024anomaly'];

    eval([Tropical_landflux_var '=Tropical_' region_name '_Mean_landflux_list(:,end)-Tropical_' region_name '_Mean_landflux_list(:,end-2);']);
    eval([Boreal_landflux_var '=Boreal_' region_name '_Mean_landflux_list(:,end)-Boreal_' region_name '_Mean_landflux_list(:,end-2);']);
    eval([North_temperate_landflux_var '=North_temperate_' region_name '_Mean_landflux_list(:,end)-North_temperate_' region_name '_Mean_landflux_list(:,end-2);']);
    eval([South_temperate_landflux_var '=South_temperate_' region_name '_Mean_landflux_list(:,end)-South_temperate_' region_name '_Mean_landflux_list(:,end-2);']);

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

Driver_Carbon_budget_2024anomaly=[Positive_Fire_Carbon_budegt_2024anomaly;Positive_TER_Carbon_budegt_2024anomaly;...
    Negative_GPP_Carbon_budegt_2024anomaly;Positive_GPP_Carbon_budegt_2024anomaly;Negative_TER_Carbon_budegt_2024anomaly];

Driver_Carbon_budget_2024anomaly_std=[Positive_Fire_Carbon_budegt_2024anomaly_std;Positive_TER_Carbon_budegt_2024anomaly_std;...
    Negative_GPP_Carbon_budegt_2024anomaly_std;Positive_GPP_Carbon_budegt_2024anomaly_std;Negative_TER_Carbon_budegt_2024anomaly_std];
regions = {'Positive_Fire', 'Positive_TER', 'Negative_GPP', 'Positive_GPP', 'Negative_TER'};

Driver_NBE_contribution=[ Boreal_Positive_Fire_landflux_2024anomaly, Boreal_Positive_TER_landflux_2024anomaly,Boreal_Negative_GPP_landflux_2024anomaly,Boreal_Positive_GPP_landflux_2024anomaly,Boreal_Negative_TER_landflux_2024anomaly;...
    North_temperate_Positive_Fire_landflux_2024anomaly, North_temperate_Positive_TER_landflux_2024anomaly,North_temperate_Negative_GPP_landflux_2024anomaly,North_temperate_Positive_GPP_landflux_2024anomaly,North_temperate_Negative_TER_landflux_2024anomaly;...
    Tropical_Positive_Fire_landflux_2024anomaly, Tropical_Positive_TER_landflux_2024anomaly,Tropical_Negative_GPP_landflux_2024anomaly,Tropical_Positive_GPP_landflux_2024anomaly,Tropical_Negative_TER_landflux_2024anomaly;...
    South_temperate_Positive_Fire_landflux_2024anomaly, South_temperate_Positive_TER_landflux_2024anomaly,South_temperate_Negative_GPP_landflux_2024anomaly,South_temperate_Positive_GPP_landflux_2024anomaly,South_temperate_Negative_TER_landflux_2024anomaly];
%%
ff=figure
set(gcf,'unit','pixels','position',[1000,583,1003,655]);

t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'tight';
a=nexttile
% part 1 *********************************************************************************
m_proj('robinson','long',[-180 180],'lat',[-90 90]); hold on
m_pcolor(lon1,lat1,dominant_factor,'linestyle','none');
m_coast('linewidth',1,'color','k'); %填充矢量边界
% m_line(bou_canX,bou_canY,'linewidth',1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[],'ytick',[],...
    'linewidth',1,'xaxisloc','bottom','yaxisloc','left','FontSize',12,'xticklabels',[],'yticklabels',[],'linestyle','none'); %对齐格网
caxis([1,6]);
mycolor=[153,32,24;
    255,225,12;
    195,122,27;
    96,185,5;
    23,80,35;
    68,179,209]/255;

colormap(mycolor);
cb=colorbar('Ticks',[1+5/12:5/6:6],...
    'TickLabels',{'+Fire','+TER','–GPP','+GPP','–TER','–Fire'},'FontName','Arial','FontSize',12,'location','southoutside')
% set(get(cb,'ylabel'),'string','Dominant contributor','fontsize',12);
% title('Dominant contributor','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.0744876165355192 1.02314626882292 0],'FontName','Arial','FontSize',18,'fontweight','bold')

nexttile
b=bar(Driver_NBE_contribution,'stacked','LineStyle','none'); hold on
b(1).FaceColor = [153,32,24]/255;
b(2).FaceColor = [255,225,12]/255;
b(3).FaceColor = [195,122,27]/255;
b(4).FaceColor = [96,185,5]/255;
b(5).FaceColor = [23,80,35]/255;
ylim([-1,2])
set(gca,'yTick',[-1:1:2],'FontName','Arial','fontsize',12)
set(gca,'XTickLabel',{'50–90°N','23–50°N','23°N–23°S','23–60°S'},'FontName','Arial','fontsize',12,'XTickLabelRotation',30)
ylabel('NBE anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
text('string','b','Units','normalized','position',[-0.103897538855785 0.982407326802595 0],'FontName','Arial','FontSize',18,'fontweight','bold')

nexttile
customColors = [
    0.17 0.42 0.7;    % 绿色
    0.0667 0.4667 0.2;    % 蓝色
    1 0.6471 0;    % 黄色
    0.698 0.1333 0.1333 % 紫色
    ];
bb = bar(Driver_Carbon_budget_2024anomaly,'FaceColor','flat','LineStyle','none'); hold on
bb(1).CData = customColors(1,:);
bb(2).CData = customColors(2,:);
bb(3).CData = customColors(3,:);
bb(4).CData = customColors(4,:);
%分组误差棒
[M,N]=size(Driver_Carbon_budget_2024anomaly);
for i=1:N
    xx(:,i)=bb(i).XEndPoints';
end
h2=errorbar(xx(:,:),Driver_Carbon_budget_2024anomaly,Driver_Carbon_budget_2024anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',5,'Marker','none')
pp=plot([xx(1,end)+(xx(2,1)-xx(1,end))/2,xx(1,end)+(xx(2,1)-xx(1,end))/2],[-4,6],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
pp=plot([xx(2,end)+(xx(2,1)-xx(1,end))/2,xx(2,end)+(xx(2,1)-xx(1,end))/2],[-4,6],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
pp=plot([xx(3,end)+(xx(2,1)-xx(1,end))/2,xx(3,end)+(xx(2,1)-xx(1,end))/2],[-4,6],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
pp=plot([xx(4,end)+(xx(2,1)-xx(1,end))/2,xx(4,end)+(xx(2,1)-xx(1,end))/2],[-4,6],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
ylim([-4,6])
set(gca,'xTicklabel',{'+Fire','+TER','–GPP','+GPP','–TER'},'FontName','Arial','fontsize',12)
set(gca,'yTick',[-4:2:6],'FontName','Arial','fontsize',12)
ylabel('Flux anomalies (PgC yr^{-1})','FontName','Arial','FontSize',12)
lgd=legend(bb,{'\DeltaNBE','\DeltaGPP','\DeltaTER','\DeltaFire'},'NumColumns',4,'FontName','Arial','FontSize',12,'Box','off','Location','southoutside')
text('string','c','Units','normalized','position',[-0.103897538855785 0.982407326802595 0],'FontName','Arial','FontSize',18,'fontweight','bold')




nexttile
customColors = [
    186 123 52;    % 红色
    58 93 126;    % 绿色
    228 79 47;    % 黄色
    183 155 36 % 紫色
    ]/255;
bb = bar(Driver_Climate_2024anomaly,'FaceColor','flat','LineStyle','none'); hold on
bb(1).CData = customColors(1,:);
bb(2).CData = customColors(2,:);
bb(3).CData = customColors(3,:);
bb(4).CData = customColors(4,:);
%分组误差棒
[M,N]=size(Driver_Climate_2024anomaly);
for i=1:N
    xx(:,i)=bb(i).XEndPoints';
end

pp=plot([xx(1,end)+(xx(2,1)-xx(1,end))/2,xx(1,end)+(xx(2,1)-xx(1,end))/2],[-5,5],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
pp=plot([xx(2,end)+(xx(2,1)-xx(1,end))/2,xx(2,end)+(xx(2,1)-xx(1,end))/2],[-5,5],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
pp=plot([xx(3,end)+(xx(2,1)-xx(1,end))/2,xx(3,end)+(xx(2,1)-xx(1,end))/2],[-5,5],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on
pp=plot([xx(4,end)+(xx(2,1)-xx(1,end))/2,xx(4,end)+(xx(2,1)-xx(1,end))/2],[-5,5],'--','LineWidth',1,'MarkerSize',3,'color','k');hold on

% xlim([0.5,4.5])
ylim([-4.5,4.5])
set(gca,'yTick',[-4:2:4],'FontName','Arial','fontsize',12)
% set(gca,'xTick',[],'FontName','Arial','fontsize',12)
set(gca,'xTicklabel',{'+Fire','+TER','–GPP','+GPP','–TER'},'FontName','Arial','fontsize',12)

ylabel('\sigma','FontName','Arial','FontSize',12)
lgd=legend(bb,{'Temperature','TWS','VPD','SR'},'NumColumns',4,'FontName','Arial','FontSize',12,'Box','off','Location','southoutside')
text('string','d','Units','normalized','position',[-0.103897538855785 0.982407326802595 0],'FontName','Arial','FontSize',18,'fontweight','bold')


z=axes('Position',[0.105040950983325 0.680916030534351 0.0863847718481807 0.170992366412214]);
bbb=bar(Ratio_result,'FaceColor', 'flat')
set(bbb, 'CData', [153,32,24;
    255,225,12;
    195,122,27;
    96,185,5;
    23,80,35;
    68,179,209]/255); % 应用颜色矩阵
box off

set(gca,'yTick',[0:10:40],'FontName','Arial','fontsize',8)
set(gca,'xTicklabel',{''},'FontName','Arial','fontsize',8)
ylabel('Area ratio (%)','FontName','Arial','FontSize',8)
set(gca, 'Color', 'none')


result=['E:\phd_file\Tropical_2024\Result\V4\Flux_anomalies_2024_drivers_map.png']
% print(result,ff,'-r600','-dpng');
