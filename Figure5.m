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
dominant_factor=importdata("E:\phd_file\Tropical_2024\domnant_factor.tif").*pixel_mask;
Positive_TER_mask=dominant_factor;
Positive_TER_mask(Positive_TER_mask~=2)=nan;
Positive_TER_mask(~isnan(Positive_TER_mask))=1;

Negative_GPP_mask=dominant_factor;
Negative_GPP_mask(Negative_GPP_mask~=3)=nan;
Negative_GPP_mask(~isnan(Negative_GPP_mask))=1;


regions = {'Global'};

for year=2015:2024

    for month=1:12

        GPP_temp=importdata(['E:\phd_file\Tropical_2024\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);
        TER_temp=importdata(['E:\phd_file\Tropical_2024\GCAS_2015-2024\TER\monthly\TER_' num2str(year) '_' num2str(month) '.tif']);
        TEM_temp=importdata(['E:\phd_file\Tropical_2024\Air_temperature\monthly\ERA5_land_air_temperature_' num2str(year) '_' num2str(month) '.tif']);
        TWS_temp=importdata(['E:\phd_file\Tropical_2024\TWS\GLDAS\month\GLDAS_TWS\GLDAS_TWS_' num2str(year) '_' num2str(month) '.tif']);


        GPP_monthly_temp(:,:,month)=GPP_temp;
        TER_monthly_temp(:,:,month)=TER_temp;
        TEM_monthly_temp(:,:,month)=TEM_temp;
        TWS_monthly_temp(:,:,month)=TWS_temp;

    end

    GPP_year(:,:,year-2014)=sum(GPP_monthly_temp(:,:,1:12),3);   
    TER_year(:,:,year-2014)=sum(TER_monthly_temp(:,:,1:12),3);
    TEM_year(:,:,year-2014)=nanmean(TEM_monthly_temp(:,:,1:12),3);
    TWS_year(:,:,year-2014)=nanmean(TWS_monthly_temp(:,:,1:12),3);

    clear GPP_monthly_temp TER_monthly_temp TEM_monthly_temp TWS_monthly_temp

end

for year=1:10

    for i=1:length(regions)

        region_name = regions{i};
        mask_var = [region_name '_mask'];

        GPP_var = [region_name '_Forest_GPP_year_Positive_TER_list'];
        TER_var = [region_name '_Forest_TER_year_Positive_TER_list'];
        TEM_var = [region_name '_Forest_TEM_year_Positive_TER_list'];
        TWS_var = [region_name '_Forest_TWS_year_Positive_TER_list'];
        
        eval([GPP_var '(year) = nansum(nansum(GPP_year(:,:,year).*area_grid .* Positive_TER_mask .* forest_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*forest_mask.*' mask_var '));']);
        eval([TER_var '(year) = nansum(nansum(TER_year(:,:,year).*area_grid .* Positive_TER_mask .*forest_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*forest_mask.*' mask_var '));']);
        eval([TEM_var '(year) = nansum(nansum(TEM_year(:,:,year).*area_grid .* Positive_TER_mask .*forest_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*forest_mask.*' mask_var '));']);
        eval([TWS_var '(year) = nansum(nansum(TWS_year(:,:,year).*area_grid .* Positive_TER_mask .*forest_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*forest_mask.*' mask_var '));']);
          
        GPP_var = [region_name '_Non_forest_GPP_year_Positive_TER_list'];
        TER_var = [region_name '_Non_forest_TER_year_Positive_TER_list'];
        TEM_var = [region_name '_Non_forest_TEM_year_Positive_TER_list'];
        TWS_var = [region_name '_Non_forest_TWS_year_Positive_TER_list'];
        
        eval([GPP_var '(year) = nansum(nansum(GPP_year(:,:,year).*area_grid .* Positive_TER_mask .* Non_forest_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*Non_forest_mask.*' mask_var '));']);
        eval([TER_var '(year) = nansum(nansum(TER_year(:,:,year).*area_grid .* Positive_TER_mask .*Non_forest_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*Non_forest_mask.*' mask_var '));']);
        eval([TEM_var '(year) = nansum(nansum(TEM_year(:,:,year).*area_grid .* Positive_TER_mask .*Non_forest_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*Non_forest_mask.*' mask_var '));']);
        eval([TWS_var '(year) = nansum(nansum(TWS_year(:,:,year).*area_grid .* Positive_TER_mask .*Non_forest_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*Non_forest_mask.*' mask_var '));']);


        GPP_var = [region_name '_Cropland_GPP_year_Positive_TER_list'];
        TER_var = [region_name '_Cropland_TER_year_Positive_TER_list'];
        TEM_var = [region_name '_Cropland_TEM_year_Positive_TER_list'];
        TWS_var = [region_name '_Cropland_TWS_year_Positive_TER_list'];
        
        eval([GPP_var '(year) = nansum(nansum(GPP_year(:,:,year).*area_grid .* Positive_TER_mask .* Cropland_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*Cropland_mask.*' mask_var '));']);
        eval([TER_var '(year) = nansum(nansum(TER_year(:,:,year).*area_grid .* Positive_TER_mask .*Cropland_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*Cropland_mask.*' mask_var '));']);
        eval([TEM_var '(year) = nansum(nansum(TEM_year(:,:,year).*area_grid .* Positive_TER_mask .*Cropland_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*Cropland_mask.*' mask_var '));']);
        eval([TWS_var '(year) = nansum(nansum(TWS_year(:,:,year).*area_grid .* Positive_TER_mask .*Cropland_mask.*' mask_var '))/nansum(nansum(area_grid .* Positive_TER_mask .*Cropland_mask.*' mask_var '));']);

    end

end

%%
f=figure
set(gcf,'unit','pixels','position',[725,391,675,601]);

t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
scatter(Global_Non_forest_TEM_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[14.5:0.2:15.2],'FontSize',12,'FontName','Arial','XLim',[14.5,15.2]);
set(gca, 'YTick',[800:50:950],'FontSize',12,'FontName','Arial','yLim',[800,950]);
xlabel('Temperature (°C)','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Non_forest_TEM_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Non_forest_TEM_year_Positive_TER_list),max(Global_Non_forest_TEM_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.63127519228908 0.166989437289201 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Non_forest_TEM_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,2);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Non_forest_TEM_year_Positive_TER_list(end),Global_Non_forest_TER_year_Positive_TER_list(end),'2024\rightarrow ','HorizontalAlignment','right','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','a','Units','normalized','position',[-0.234795589036422 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')

nexttile
scatter(Global_Non_forest_TWS_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[1085:10:1115],'FontSize',12,'FontName','Arial','xLim',[1085,1118]);
set(gca, 'YTick',[800:50:950],'FontSize',12,'FontName','Arial','yLim',[800,950]);
xlabel('TWS (cm)','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Non_forest_TWS_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Non_forest_TWS_year_Positive_TER_list),max(Global_Non_forest_TWS_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.63127519228908 0.166989437289201 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Non_forest_TWS_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,2);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Non_forest_TWS_year_Positive_TER_list(end),Global_Non_forest_TER_year_Positive_TER_list(end),'2024\rightarrow ','HorizontalAlignment','right','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','b','Units','normalized','position',[-0.234795589036422 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')


nexttile
scatter(Global_Non_forest_GPP_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[800:50:950],'FontSize',12,'FontName','Arial','xLim',[800,950]);
set(gca, 'YTick',[800:50:950],'FontSize',12,'FontName','Arial','yLim',[800,950]);
xlabel('GPP (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Non_forest_GPP_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,1);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Non_forest_GPP_year_Positive_TER_list),max(Global_Non_forest_GPP_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.63127519228908 0.166989437289201 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Non_forest_GPP_year_Positive_TER_list,Global_Non_forest_TER_year_Positive_TER_list,1);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Non_forest_GPP_year_Positive_TER_list(end),Global_Non_forest_TER_year_Positive_TER_list(end),'2024\rightarrow ','HorizontalAlignment','right','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','c','Units','normalized','position',[-0.234795589036422 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')
scatter([0:3:2000],[0:3:2000],'filled','MarkerFaceColor','k','SizeData',2); hold on

nexttile
% 标准化（Z-score）
Global_Non_forest_TWS_norm = normalize(Global_Non_forest_TWS_year_Positive_TER_list','zscore');
Global_Non_forest_TEM_norm = normalize(Global_Non_forest_TEM_year_Positive_TER_list','zscore');
Global_Non_forest_GPP_norm = normalize(Global_Non_forest_GPP_year_Positive_TER_list','zscore');

Global_Non_forest_X = [Global_Non_forest_TWS_norm.^2, Global_Non_forest_TEM_norm.^2, Global_Non_forest_TWS_norm, Global_Non_forest_TEM_norm, Global_Non_forest_GPP_norm, ones(size(Global_Non_forest_TER_year_Positive_TER_list'))];
% 拟合参数
Global_Non_forest_beta = Global_Non_forest_X \ Global_Non_forest_TER_year_Positive_TER_list';
% 计算预测值
Global_Non_forest_y_pred = Global_Non_forest_X * Global_Non_forest_beta;

Global_Non_forest_Climate_beta=Global_Non_forest_beta;Global_Non_forest_Climate_beta(5)=0;
Global_Non_forest_TER_climate_predict_result= Global_Non_forest_X * Global_Non_forest_Climate_beta;
Global_Non_forest_TER_predict_2024anomalies=Global_Non_forest_y_pred(end)-Global_Non_forest_y_pred(end-2);
Global_Non_forest_TER_climate_predict_2024anomalies=Global_Non_forest_TER_climate_predict_result(end)-Global_Non_forest_TER_climate_predict_result(end-2);
Global_Non_forest_result=[Global_Non_forest_TER_climate_predict_2024anomalies,Global_Non_forest_TER_predict_2024anomalies-Global_Non_forest_TER_climate_predict_2024anomalies];
b=bar(Global_Non_forest_result,'FaceColor','flat')
b.CData(1,:) = [0 204 255]/255;
b.CData(2,:) = [102 204 51]/255;
set(gca, 'YTick',[0:20:80],'FontSize',12,'FontName','Arial','yLim',[0,80]);

set(gca,'XTickLabel',{'Climate','GPP'},'FontName','Arial','fontsize',12)
ylabel('TER anomalies (gC m^{-2})','FontName','Arial','FontSize',12);
text('string','d','Units','normalized','position',[-0.234795589036422 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')

result=['E:\phd_file\Tropical_2024\Result\V4\Global_Non_forest_response_scatter_bar.png']
% print(result,f,'-r600','-dpng');

Global_Non_forest_result(1)/sum(Global_Non_forest_result)
Global_Non_forest_result(2)/sum(Global_Non_forest_result)

%%
f=figure
set(gcf,'unit','pixels','position',[725,391,675,601]);

t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
scatter(Global_Forest_TEM_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[13:0.4:14.2],'FontSize',12,'FontName','Arial','XLim',[13,14.2]);
set(gca, 'YTick',[1680:40:1820],'FontSize',12,'FontName','Arial','yLim',[1680,1820]);
xlabel('Temperature (°C)','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Forest_TEM_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Forest_TEM_year_Positive_TER_list),max(Global_Forest_TEM_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.63127519228908 0.166989437289201 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Forest_TEM_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,2);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Forest_TEM_year_Positive_TER_list(end),Global_Forest_TER_year_Positive_TER_list(end),'2024\rightarrow ','HorizontalAlignment','right','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','a','Units','normalized','position',[-0.274321280147048 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')

nexttile
scatter(Global_Forest_TWS_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[1555:15:1605],'FontSize',12,'FontName','Arial','xLim',[1555 1605]);
set(gca, 'YTick',[1680:40:1820],'FontSize',12,'FontName','Arial','yLim',[1680,1820]);
xlabel('TWS (cm)','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Forest_TWS_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Forest_TWS_year_Positive_TER_list),max(Global_Forest_TWS_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.0734239526196586 0.166989437289201 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Forest_TWS_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,2);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.0734239526196586 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.0734239526196586 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Forest_TWS_year_Positive_TER_list(end),Global_Forest_TER_year_Positive_TER_list(end),' \leftarrow2024','HorizontalAlignment','left','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','b','Units','normalized','position',[-0.274321280147048 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')


nexttile
scatter(Global_Forest_GPP_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[1680:40:1840],'FontSize',12,'FontName','Arial','xLim',[1680,1840]);
set(gca, 'YTick',[1680:40:1840],'FontSize',12,'FontName','Arial','yLim',[1680,1840]);
xlabel('GPP (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Forest_GPP_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,1);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Forest_GPP_year_Positive_TER_list),max(Global_Forest_GPP_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.63127519228908 0.218937489237253 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Forest_GPP_year_Positive_TER_list,Global_Forest_TER_year_Positive_TER_list,1);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Forest_GPP_year_Positive_TER_list(end),Global_Forest_TER_year_Positive_TER_list(end),'2024\rightarrow ','HorizontalAlignment','right','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','c','Units','normalized','position',[-0.274321280147048 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')
scatter([0:3:2000],[0:3:2000],'filled','MarkerFaceColor','k','SizeData',2); hold on

nexttile
% 标准化（Z-score）
Global_Forest_TWS_norm = normalize(Global_Forest_TWS_year_Positive_TER_list','zscore');
Global_Forest_TEM_norm = normalize(Global_Forest_TEM_year_Positive_TER_list','zscore');
Global_Forest_GPP_norm = normalize(Global_Forest_GPP_year_Positive_TER_list','zscore');

Global_Forest_X = [Global_Forest_TWS_norm.^2, Global_Forest_TEM_norm.^2, Global_Forest_TWS_norm, Global_Forest_TEM_norm, Global_Forest_GPP_norm, ones(size(Global_Forest_TER_year_Positive_TER_list'))];
% 拟合参数
Global_Forest_beta = Global_Forest_X \ Global_Forest_TER_year_Positive_TER_list';
% 计算预测值
Global_Forest_y_pred = Global_Forest_X * Global_Forest_beta;

Global_Forest_Climate_beta=Global_Forest_beta;Global_Forest_Climate_beta(5)=0;
Global_Forest_TER_climate_predict_result= Global_Forest_X * Global_Forest_Climate_beta;
Global_Forest_TER_predict_2024anomalies=Global_Forest_y_pred(end)-Global_Forest_y_pred(end-2);
Global_Forest_TER_climate_predict_2024anomalies=Global_Forest_TER_climate_predict_result(end)-Global_Forest_TER_climate_predict_result(end-2);
Global_Forest_result=[Global_Forest_TER_climate_predict_2024anomalies,Global_Forest_TER_predict_2024anomalies-Global_Forest_TER_climate_predict_2024anomalies];
b=bar(Global_Forest_result,'FaceColor','flat')
b.CData(1,:) = [0 204 255]/255;
b.CData(2,:) = [102 204 51]/255;
set(gca, 'YTick',[0:20:60],'FontSize',12,'FontName','Arial','yLim',[0,60]);

set(gca,'XTickLabel',{'Climate','GPP'},'FontName','Arial','fontsize',12)
ylabel('TER anomalies (gC m^{-2})','FontName','Arial','FontSize',12);
text('string','d','Units','normalized','position',[-0.274321280147048 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')

result=['E:\phd_file\Tropical_2024\Result\V4\Global_Forest_response_scatter_bar.png']
% print(result,f,'-r600','-dpng');

Global_Forest_result(1)/sum(Global_Forest_result)

%%
f=figure
set(gcf,'unit','pixels','position',[725,391,675,601]);

t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
scatter(Global_Cropland_TEM_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[16.7:0.3:17.6],'FontSize',12,'FontName','Arial','XLim',[16.7,17.6]);
set(gca, 'YTick',[1150:60:1390],'FontSize',12,'FontName','Arial','yLim',[1150,1390]);
xlabel('Temperature (°C)','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Cropland_TEM_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Cropland_TEM_year_Positive_TER_list),max(Global_Cropland_TEM_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.63127519228908 0.166989437289201 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Cropland_TEM_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,2);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Cropland_TEM_year_Positive_TER_list(end),Global_Cropland_TER_year_Positive_TER_list(end),'2024\rightarrow ','HorizontalAlignment','right','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','a','Units','normalized','position',[-0.274321280147048 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')

nexttile
scatter(Global_Cropland_TWS_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[1160:10:1200],'FontSize',12,'FontName','Arial','xLim',[1160,1200]);
set(gca, 'YTick',[1150:60:1390],'FontSize',12,'FontName','Arial','yLim',[1150,1390]);
xlabel('TWS (cm)','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Cropland_TWS_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Cropland_TWS_year_Positive_TER_list),max(Global_Cropland_TWS_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.684994200553543 0.166989437289201 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Cropland_TWS_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,2);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.684994200553543 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.684994200553543 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Cropland_TWS_year_Positive_TER_list(end),Global_Cropland_TER_year_Positive_TER_list(end),'2024\rightarrow ','HorizontalAlignment','right','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','b','Units','normalized','position',[-0.274321280147048 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')


nexttile
scatter(Global_Cropland_GPP_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
box on
set(gca, 'XTick',[1160:60:1400],'FontSize',12,'FontName','Arial','xLim',[1160,1400]);
set(gca, 'YTick',[1160:60:1400],'FontSize',12,'FontName','Arial','yLim',[1160,1400]);
xlabel('GPP (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
ylabel('TER (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',12);
% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(Global_Cropland_GPP_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,1);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(Global_Cropland_GPP_year_Positive_TER_list),max(Global_Cropland_GPP_year_Positive_TER_list));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.63127519228908 0.166989437289201 0],'FontName','Arial','FontSize',12)
P_value=F_test(Global_Cropland_GPP_year_Positive_TER_list,Global_Cropland_TER_year_Positive_TER_list,1);
P_value=fix(P_value * 100) / 100;
if P_value>0.01
    P_text=['\itP\rm = ' num2str(P_value)];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
else
    P_text=['\itP\rm < 0.01'];
    text('string',P_text,'Units','normalized','position',[0.63127519228908 0.0687438232541132 0],'FontName','Arial','FontSize',12)
end
text(Global_Cropland_GPP_year_Positive_TER_list(end),Global_Cropland_TER_year_Positive_TER_list(end),'2024\rightarrow ','HorizontalAlignment','right','FontSize',12,'FontName','Arial','fontweight','bold','Color','k')
text('string','c','Units','normalized','position',[-0.274321280147048 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')
scatter([0:3:2000],[0:3:2000],'filled','MarkerFaceColor','k','SizeData',2); hold on

nexttile
% 标准化（Z-score）
Global_Cropland_TWS_norm = normalize(Global_Cropland_TWS_year_Positive_TER_list','zscore');
Global_Cropland_TEM_norm = normalize(Global_Cropland_TEM_year_Positive_TER_list','zscore');
Global_Cropland_GPP_norm = normalize(Global_Cropland_GPP_year_Positive_TER_list','zscore');

Global_Cropland_X = [Global_Cropland_TWS_norm.^2, Global_Cropland_TEM_norm.^2, Global_Cropland_TWS_norm, Global_Cropland_TEM_norm, Global_Cropland_GPP_norm, ones(size(Global_Cropland_TER_year_Positive_TER_list'))];
% 拟合参数
Global_Cropland_beta = Global_Cropland_X \ Global_Cropland_TER_year_Positive_TER_list';
% 计算预测值
Global_Cropland_y_pred = Global_Cropland_X * Global_Cropland_beta;

Global_Cropland_Climate_beta=Global_Cropland_beta;Global_Cropland_Climate_beta(5)=0;
Global_Cropland_TER_climate_predict_result= Global_Cropland_X * Global_Cropland_Climate_beta;
Global_Cropland_TER_predict_2024anomalies=Global_Cropland_y_pred(end)-Global_Cropland_y_pred(end-2);
Global_Cropland_TER_climate_predict_2024anomalies=Global_Cropland_TER_climate_predict_result(end)-Global_Cropland_TER_climate_predict_result(end-2);
Global_Cropland_result=[Global_Cropland_TER_climate_predict_2024anomalies,Global_Cropland_TER_predict_2024anomalies-Global_Cropland_TER_climate_predict_2024anomalies];
b=bar(Global_Cropland_result,'FaceColor','flat')
b.CData(1,:) = [0 204 255]/255;
b.CData(2,:) = [102 204 51]/255;
set(gca, 'YTick',[0:20:80],'FontSize',12,'FontName','Arial','yLim',[0,80]);

set(gca,'XTickLabel',{'Climate','GPP'},'FontName','Arial','fontsize',12)
ylabel('TER anomalies (gC m^{-2})','FontName','Arial','FontSize',12);
text('string','d','Units','normalized','position',[-0.274321280147048 1.05272973700269 0],'FontName','Arial','FontSize',18,'fontweight','bold')

result=['E:\phd_file\Tropical_2024\Result\V4\Global_Cropland_response_scatter_bar.png']
% print(result,f,'-r600','-dpng');

Global_Cropland_result(1)/sum(Global_Cropland_result)
