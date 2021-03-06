clear all;
close all;
clc;

%FileNames  ---------------------------------------------------------------------------------------------------------------------------
logFileName = 'log.xlsx';
xrdFileName = 'xrd.xls';
griFileName = 'gri.xlsx';
sraFileName = 'sra.xls';

%mineral densities  -------------------------------------------------------------------------------------------------------------------
qtz=2.65; kfeld=2.52; pfeld=2.66;              %silicates
cal=2.71; dol=2.87;                              %carbonates
pyr=4.99; marc=4.87;                             %heavies
ill_smec=2.60; ill_mic=2.75; kaol=2.60; chl=2.94;  %clays
densities=[qtz, kfeld, pfeld, cal, dol, pyr, marc, ill_smec, ill_mic, kaol,chl];
%Constants-----------------------------------------------------------------------------------------------------------------------------
densityKerogen = 1.35;
densitySand = 2.65;
densityShale = 2.71;
densityFluid = 1;
%Indexes in files----------------------------------------------------------------------------------------------------------------------
logDepthIndex = 1;
logDreshIndex = 2;
logDtcIndex = 3;
logEcgrIndex = 4;
logNphiIndex = 5;
logRhobIndex = 6;
logUranIndex = 7;
xrdTocIndex = 3;
xrdDepthIndex = 1;
griDepthIndex = 2;
griBulkDensityIndex = 3;
griGrainDensityIndex = 7;
griPorosityIndex = 8;
griSaturationIndex=10;
sraDepthIndex=1;
sraTmaxIndex=6;
sraHIIndex=8;
sraOIIndex=9;
%NonClay Indexes----------------------------------------------------------------------------------------------------------------------
silicateIndex = [1 2 3];
carbonateIndex = [4 5];
heaviesIndex = [6 7];
clayIndex = [8 9 10 11];
%Plot ranges--------------------------------------------------------------------------------------------------------------------------
logRange = [15998 17603];
logEcgrUrXaxisRange = [0 400];
logDepthRange = [2800 3060];
logNphiXaxisRange = [-0.15 0.45];
logRhobXaxisRange = [1.95 2.95];
logDtcXaxisRange = [0 200];
logDreshXaxisRange = [0.1 1000];
logVshaleXaxisRange = [-1 2];
xrdNonClayRange = [4 10];
xrdClayRange = [11 14];
xrdTocPasseyXaxisRange = [-1 5];
griXrdPlotRange = [2 3];
porosityWithoutKeroXaxisRange = [-0.1 0.4];
saturationRange = [0 1];
picketplotXaxisRange=[0.1 1000];
picketplotYaxisRange=[0.01 1];
%Parameter Constants------------------------------------------------------------------------------------------------------------------
logSandLineHist = 80;
logShaleLineHist = 190;
clayFactor = 0.45;
dReshBaseHist = 35;
dtcBaseHist = 60;
levelOfMaturity = 12;
scalingFactor1 = 70;
scalingFactor2 = 0.75;
scalingFactor3 = 10;
scalingFactor4 = 0.75;
saturationExponentN=2;
TortuosityFactorA=0.65;
%-------------------------------------------------------------------------------------------------------------------------------------
LOG = xlsread(logFileName);
XRD = xlsread(xrdFileName);
GRI = xlsread(griFileName);
SRA = xlsread(sraFileName);

logRhob = LOG(logRange(1,1):logRange(1,2),logRhobIndex);
xrdToc = XRD(:,xrdTocIndex);
weightPercentKerogen = 1.1*xrdToc;
weightPercentKerogenNormFactor = (1-weightPercentKerogen/100);

%-------------------------------------------------------------------------------------------------------------------------------------
xrdDepth = XRD(:,xrdDepthIndex);
xrdNonClayWeightPercent = XRD(:,xrdNonClayRange(1,1):xrdNonClayRange(1,2));
xrdClayWeightPercent = XRD(:,xrdClayRange(1,1):xrdClayRange(1,2));
weightPercentCombine = cat(2,xrdNonClayWeightPercent, xrdClayWeightPercent);
numberOfMinerals = size(weightPercentCombine,2);

weightPercentsNormalized = zeros(size(weightPercentCombine));
weightByDensityNormalized = zeros(size(weightPercentCombine));
XRDGrainDensityWithoutKerogen = zeros(size(xrdClayWeightPercent,1),1);
XRDGrainDensityWithKerogen = zeros(size(xrdClayWeightPercent,1),1);

for i=1:size(weightPercentCombine,1) 
    temp = 0;
    for j=1:size(weightPercentCombine,2)    
   		weightPercentsNormalized(i,j) = weightPercentKerogenNormFactor(i,1).*weightPercentCombine(i,j);
   		weightByDensityNormalized(i,j) = weightPercentsNormalized(i,j)./densities(1,j);
   		temp = temp + weightByDensityNormalized(i,j);
    end
    XRDGrainDensityWithoutKerogen(i,1) = 100/temp;
    XRDGrainDensityWithKerogen(i,1) = 100/(temp + weightPercentKerogen(i,1)/densityKerogen);
end
xrdNonClayWeightNormSum = sum(weightPercentsNormalized(:,1:size(xrdNonClayWeightPercent,2)),2);
xrdClayWeightNormSum = sum(weightPercentsNormalized(:,(size(xrdNonClayWeightPercent,2)+1):size(weightPercentCombine,2)),2);

%-------------------------------------------------------------------------------------------------------------------------------------
griDepth = GRI(:,griDepthIndex);
griGrainDensity = GRI(:,griGrainDensityIndex);
griBulkDensity = GRI(:,griBulkDensityIndex);
griPorosity = GRI(:,griPorosityIndex);
griSaturation=GRI(:,griSaturationIndex);

c=[];
griXrdCommonWeightPercentage = [];
cIndex = 0;
for j= 1:length(GRI)       
    for  k= 1:length(XRD)
        if griDepth(j,1)==xrdDepth(k,1) %grain density xrd vs gri calibration at same depth data
        	cIndex = cIndex+1;
            c(cIndex,1) = griDepth(j,1);
            c(cIndex,2) = griGrainDensity(j,1);
            c(cIndex,3) = XRDGrainDensityWithKerogen(k,1);
            c(cIndex,4) = xrdNonClayWeightNormSum(k,1);
            c(cIndex,5) = xrdClayWeightNormSum(k,1);
            c(cIndex,6) = weightPercentKerogen(k,1);
            c(cIndex,8) = griBulkDensity(j,1);
            c(cIndex,9:8+numberOfMinerals) = weightByDensityNormalized(k,:);
            c(cIndex,9+numberOfMinerals) = XRDGrainDensityWithoutKerogen(k,1);  
            c(cIndex,9+numberOfMinerals+1) = griPorosity(j,1);
            griXrdCommonWeightPercentage(cIndex,:) = weightPercentsNormalized(k,:);
            c(cIndex,9+numberOfMinerals+2) = griSaturation(j,1);
        end
    end 
end
c(:,7) = c(:,6)/1.1;

%-------------------------------------------------------------------------------------------------------------------------------------
% start different figure
figure('units','normalized','outerposition',[0 0 1 1])

logdepth = LOG(logRange(1,1):logRange(1,2),logDepthIndex);
logECGR = LOG(logRange(1,1):logRange(1,2),logEcgrIndex);
UR = LOG(logRange(1,1):logRange(1,2),logUranIndex);

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,1) 	%GR
plot(logECGR,logdepth) %ECGR
axis ij
xlim([logEcgrUrXaxisRange(1,1) logEcgrUrXaxisRange(1,2)]);
ylim([logDepthRange(1,1) logDepthRange(1,2)]);
xlabel('Gamma Ray (gAPI)')
ylabel('Depth (meters)')
legend('Gamma Ray')
ax1 = gca;
ax1_pos = get(ax1,'Position'); % position of first axes
ax2 = axes('Position',[ax1_pos(1,1) ax1_pos(1,2) 0.94*ax1_pos(1,3) ax1_pos(1,4)],...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
hold on 
plot(UR,logdepth,'Parent',ax2,'Color','r')%UR
set(ax2,'XColor','r');
set(ax2,'YColor','r');
set(ax2,'YTick',[]);
axis ij
xlabel('Uranium')
legend('Uranium')
xlim([logEcgrUrXaxisRange(1,1) logEcgrUrXaxisRange(1,2)]);
ylim([logDepthRange(1,1) logDepthRange(1,2)]);
format long
hold on 	

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,2)%grain density

plot(c(:,2),c(:,1),'ok')
hold on 
plot(c(:,3),c(:,1),'*r')  %XRDGrainDensityWithKerogen
axis ij
xlim([griXrdPlotRange(1,1) griXrdPlotRange(1,2)]);
ylim([logDepthRange(1,1) logDepthRange(1,2)]);
set(gca,'YTick',[]);
xlabel('grain den')
format long 
hold on
legend('GRI','XRD')

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,3)%neutron poristy and bulk density 

logNphi = LOG(logRange(1,1):logRange(1,2),logNphiIndex);
plot(logNphi,logdepth,'g')%neutron
xlim([logNphiXaxisRange(1,1) logNphiXaxisRange(1,2)])
ylim([logDepthRange(1,1) logDepthRange(1,2)])
axis ij
set(gca,'XDir','reverse')
set(gca,'YTick',[]);
xlabel('density & neutron')
format long

axRhob = gca;
ax1_pos = get(axRhob,'Position'); 
ax2 = axes('Position',[ax1_pos(1,1) ax1_pos(1,2) 0.94*ax1_pos(1,3) ax1_pos(1,4)],...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
hold on
plot(logRhob,logdepth,'parent',ax2,'color','r')%densitclosy
axis ij
xlim([logRhobXaxisRange(1,1) logRhobXaxisRange(1,2)]);
ylim([logDepthRange(1,1) logDepthRange(1,2)]);
set(ax2,'XColor','r');
set(ax2,'YColor','r');
set(ax2,'YTick',[]);
hold on 
axis ij 
legend('RHOB')

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,4)%sonic and resis passeys

logDtc = LOG(logRange(1,1):logRange(1,2),logDtcIndex);
plot(logDtc, logdepth,'g')
xlim([logDtcXaxisRange(1,1) logDtcXaxisRange(1,2)])
ylim([logDepthRange(1,1) logDepthRange(1,2)])
axis ij
set(gca, 'XDir','reverse')
set(gca,'YTick',[]);
format long

ax1 = gca;
ax1_pos = get(ax1,'Position'); % position of first axes
ax2 = axes('Position',[ax1_pos(1,1) ax1_pos(1,2) 0.94*ax1_pos(1,3) ax1_pos(1,4)],...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
hold on 
logDresh = LOG(logRange(1,1):logRange(1,2),logDreshIndex);
semilogx(logDresh, logdepth, 'Color', 'r')%deep resistivity.
axis([logDreshXaxisRange(1,1) logDreshXaxisRange(1,2) logDepthRange(1,1) logDepthRange(1,2)])
set(ax2,'xscale','log');
set(ax2,'YTick',[]);
axis ij
format long
legend('sonic')
xlabel('Sonic and resistivity')
hold on 

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,5)%vshale

logVshale = (logECGR-logSandLineHist)./(logShaleLineHist-logSandLineHist);
logVshale(logVshale<=0) = 0.001; 
logVshale(logVshale>1) = 0.999; 
plot(logVshale,logdepth)
xlabel('Volume Percentages(%)')
ylabel('Depth (meters)')
xlim([logVshaleXaxisRange(1,1) logVshaleXaxisRange(1,2)])
ylim([logDepthRange(1,1) logDepthRange(1,2)])
%set(gca,'YTick',[]);
axis ij
logVclay = clayFactor.*logVshale;
hold on 
plot(logVclay,logdepth,'r')
for i=1:1:length(c)
    mineralVolumePercentage(i,:)= c(i,8)*c(i, 9:8+numberOfMinerals);
end
nonClayVolumePercentage = mineralVolumePercentage(:,1:size(xrdNonClayWeightPercent,2));
nonClayVolumePercentageSum = sum(nonClayVolumePercentage,2);
crystalsVolumePercentageSum  = sum(horzcat(mineralVolumePercentage(:,silicateIndex), mineralVolumePercentage(:,carbonateIndex)),2);
heaviesVolumePercentageSum  = sum(mineralVolumePercentage(:, heaviesIndex),2);
clayVolumePercentage = mineralVolumePercentage(:,(size(xrdNonClayWeightPercent,2)+1):size(mineralVolumePercentage,2));
clayVolumePercentageSum = sum(clayVolumePercentage,2);
clayVolumeSum = clayVolumePercentageSum/100;
hold on 
plot(clayVolumeSum, c(:,1),'ok')
legend('Vsh calculated from GR log','Vcl calculated from Vsh','Vclay converted from XRD wt %')
title('Vsh , Vclay Calculated and Vclay measured')

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,6)%TOC passeys and XRD;

deltaLogR = log(logDresh./dReshBaseHist)+0.02.*(logDtc - dtcBaseHist);
tocPassey = scalingFactor1.*deltaLogR.*10.^(0.297-0.1688.*levelOfMaturity)+scalingFactor2;
plot(tocPassey, logdepth,'r')
xlim([xrdTocPasseyXaxisRange(1,1) xrdTocPasseyXaxisRange(1,2)])
ylim([logDepthRange(1,1) logDepthRange(1,2)])
hold on 
plot(c(:,7),c(:,1),'ok')
set(gca,'YTick',[]);
axis ij 
hold on 
xlabel('TOC')
legend('TOC_Passey','TOC_XRD')

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,7)%bulk density

plot(logRhob,logdepth,'k')
xlabel('bulk density')
xlim([logRhobXaxisRange(1,1) logRhobXaxisRange(1,2)]);
ylim([logDepthRange(1,1) logDepthRange(1,2)]);
set(gca,'YTick',[]);
axis ij
hold on  
plot(c(:,8),c(:,1),'og')
format long
legend('Log blkd','GRI')

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,8)%porosity

porosityWithoutKerogen = (logRhob-(densityShale.*logVshale+densitySand.*(1-logVshale)))./(densityFluid-(densityShale.*logVshale+densitySand.*(1-logVshale)));
plot(porosityWithoutKerogen, logdepth,'k')
xlabel('phi')
xlim([porosityWithoutKeroXaxisRange(1,1) porosityWithoutKeroXaxisRange(1,2)])
ylim([logDepthRange(1,1) logDepthRange(1,2)])
axis ij
format long
hold on
plot(c(:,9+numberOfMinerals+1),c(:,1),'og')
set(gca,'YTick',[]);
format long
legend('phi-w/o-K','GRI')

%-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,9) %saturation from logs and core (gri)

plot(c(:,9+numberOfMinerals+2),c(:,1),'o')
xlim([saturationRange(1,1) saturationRange(1,2)])
ylim([logDepthRange(1,1) logDepthRange(1,2)])
xlabel('Saturation')
set(gca,'YTick',[]);
axis ij
format long
hold on 

% %-------------------------------------------------------------------------------------------------------------------------------------
subplot (1,10,10) %display clay non clays and kerogen

weightInGroups = zeros(length(griXrdCommonWeightPercentage),4);
for i=1:1:length(griXrdCommonWeightPercentage)
    for j=1:1:length(silicateIndex)
        weightInGroups(i,1) = weightInGroups(i,1) + griXrdCommonWeightPercentage(i, silicateIndex(1,j));
    end
    for j=1:1:length(carbonateIndex)
        weightInGroups(i,2) = weightInGroups(i,2) + griXrdCommonWeightPercentage(i, carbonateIndex(1,j));
    end
    for j=1:1:length(clayIndex)
        weightInGroups(i,3) = weightInGroups(i,3) + griXrdCommonWeightPercentage(i,clayIndex(1,j));
    end
    for j=1:1:length(heaviesIndex)
        weightInGroups(i,4) = weightInGroups(i,4) + griXrdCommonWeightPercentage(i,heaviesIndex(1,j));
    end
end

h=area(c(:,1), weightInGroups);
axis([logDepthRange(1,1) logDepthRange(1,2) 0 100])
view(90,90)
set(h(1),'FaceColor',[1 1 0]);
set(h(2),'FaceColor',[0.2 0.4 1]);
set(h(3),'FaceColor',[0.6 0.8 0.6]);
set(h(4),'FaceColor',[1 0 0]);
set(gca,'XTick',[]);
%[~,h_legend] = legend('Silicates','Carbonates','Clay','Heavies');
%PatchInLegend = findobj(h_legend, 'type', 'patch');
%set(PatchInLegend(1), 'FaceAlpha', 0.2);
%set(PatchInLegend(2), 'FaceAlpha', 0.4);
%set(PatchInLegend(3), 'FaceAlpha', 0.2);
%set(PatchInLegend(4), 'FaceAlpha', 0.4);   
%-------------------------------------------------------------------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,3,1)
format long g
hold on 
plot (c(:,2),c(:,3),'o')
xlim(griXrdPlotRange);
ylim(griXrdPlotRange);
xlabel('Measured Grain Density(gm/cc)')
ylabel('Calculated Grain Density(gm/cc)')

hold on 
%y=x line
x=griXrdPlotRange(1,1):0.1:griXrdPlotRange(1,2);
y=griXrdPlotRange(1,1):0.1:griXrdPlotRange(1,2);
plot(x,y)
title('Comparison b/w Measured and Calculated Grain Density ');

%-------------------------------------------------------------------------------------------------------------------------------------
subplot(2,3,2)

% griXrdCommonSilicatesCarbonates = sum(griXrdCommonWeightPercentage(:,horzcat(silicateIndex,carbonateIndex)),2);
% sumgriXrdClays = sum(griXrdCommonWeightPercentage(:,clayIndex),2);
plot (crystalsVolumePercentageSum, clayVolumePercentageSum, 'o');
xlabel('Core Crystals (Volume%)')
ylabel('Core Clays (Volume%)')
polyfitClaySilicateCarbonate = polyfit(crystalsVolumePercentageSum, clayVolumePercentageSum, 1);
func_1 = polyval(polyfitClaySilicateCarbonate,crystalsVolumePercentageSum);
hold on
plot(crystalsVolumePercentageSum,func_1,'--r')
str = strcat('y =  ',num2str(polyfitClaySilicateCarbonate(1)),'*x + ',num2str(polyfitClaySilicateCarbonate(2)));
title(str);
legend('Data at common Depth', 'Fitted Linear Trend')
hold off
constantA = polyfitClaySilicateCarbonate(2);
constantB = polyfitClaySilicateCarbonate(1);

%-------------------------------------------------------------------------------------------------------------------------------------
subplot(2,3,3)

% sumHeavies = sum(griXrdCommonWeightPercentage(:,heaviesIndex),2);
plot(heaviesVolumePercentageSum,c(:,9+numberOfMinerals),'o');  %heavies vs grain density plotted
xlabel('Core Heavies (Volume%)')
ylabel('Calculated GrainDensity w/o Kerogen')
polyfitheaviesGrainDensity = polyfit(heaviesVolumePercentageSum,c(:,9+numberOfMinerals), 1);
func_2 = polyval(polyfitheaviesGrainDensity,heaviesVolumePercentageSum);
hold on
plot(heaviesVolumePercentageSum,func_2,'--r')
str = strcat('y =  ',num2str(polyfitheaviesGrainDensity(1)),'*x + ',num2str(polyfitheaviesGrainDensity(2)));
title(str);
legend('Data at common Depth', 'Fitted Linear Trend')

hold off
constantC = polyfitheaviesGrainDensity(2);
constantD = polyfitheaviesGrainDensity(1);

%-------------------------------------------------------------------------------------------------------------------------------------
subplot(2,3,4)

blkdXrdCommonNorm = c(:,8);
tocCommonNorm = c(:,6)./1.1;
plot (blkdXrdCommonNorm, tocCommonNorm,'o') ;
xlabel('blkdXrdCommonNorm')
ylabel('tocCommonNorm')
polyfitblkdXrdToc = polyfit(blkdXrdCommonNorm,tocCommonNorm, 1);
func_3 = polyval(polyfitblkdXrdToc,blkdXrdCommonNorm);
hold on
plot(blkdXrdCommonNorm,func_3,'--r')
str = strcat('y =  ',num2str(polyfitblkdXrdToc(1)),'*x + ',num2str(polyfitblkdXrdToc(2)));
title(str);
hold off

figure(1)
subplot (1,10,6)
hold on
tocRohb = polyfitblkdXrdToc(1).*logRhob + polyfitblkdXrdToc(2);
%tocMix = (scalingFactor3/levelOfMaturity).*deltaLogR + scalingFactor4./logRhob;
plot(tocRohb, logdepth, 'b') %toc RHOB
hold on
%plot(tocMix, logdepth,'g')
legend('TOC_Passey','TOC_XRD','TOC_RHOB')
format long


%-------------------------------------------------------------------------------------------------------------------------------------
figure(2)
subplot(2,3,5)

xrdUrCommon = [];
index = 0;
for j= 1:length(UR)       
    for  k= 1:length(XRD)
        if  round(logdepth(j,1)*10)/10 == round(xrdDepth(k,1)*10)/10 %grain density xrd vs gri calibration at same depth data
            index = index+1;
            xrdUrCommon(index,1) = UR(j,1);
            xrdUrCommon(index,2) = xrdToc(k,1);
        end
    end 
end
plot (xrdUrCommon(:,2), xrdUrCommon(:,1), 'o');
xlabel('Core TOC (%)')
ylabel('Uranium')
title('CrossPlot : Uranium log with Core TOC')
polyfitXrdUr = polyfit(xrdUrCommon(:,2), xrdUrCommon(:,1), 1);
func_4 = polyval(polyfitXrdUr,xrdUrCommon(:,2));
hold on
plot(xrdUrCommon(:,2),func_4,'--r')
str = strcat('y =  ',num2str(polyfitXrdUr(1)),'*x + ',num2str(polyfitXrdUr(2)));
title(str);
hold off

%-------------------------------------------------------------------------------------------------------------------------------------
subplot(2,3,6)
plot(tocPassey, logdepth,'r')
xlim([xrdTocPasseyXaxisRange(1,1) xrdTocPasseyXaxisRange(1,2)])
ylim([logDepthRange(1,1) logDepthRange(1,2)])
hold on 
plot(c(:,7),c(:,1),'ok')
axis ij 
hold on 
plot(tocRohb, logdepth, 'b') %toc RHOB
hold on 
%plot(tocMix, logdepth,'g')
legend('TOC_Passey','TOC_XRD','TOC_RHOB')

%-------------------------------------------------------------------------------------------------------------------------------------
figure(1)
subplot (1,10,8)

averageGrainDensity = 2.72;
% clayWeightpercentUpscaled = (logVclay./logRhob).*averageGrainDensity.*100;
densityMA = (100 - (1.1/densityKerogen).*tocPassey.*logRhob + (constantC/constantD + constantA/constantB) - logVclay/constantB - logVclay)*constantD;

numeratorPhi = densityMA - logRhob.*((densityMA.*tocPassey./100)/densityKerogen - tocPassey./100 + 1);
denominatorPhi = densityMA - densityFluid + densityFluid.*tocPassey./100.*(1 - densityMA./densityKerogen);
porosityWithKerogen = numeratorPhi./denominatorPhi;

plot(porosityWithKerogen, logdepth,'r')
legend('phi-w/o-K','GRI','phi-w-K')
format long


%-------------------------------------------------------------------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1])
logDepthRangepicket=[1000 3100];
logdepthpicketvaluerange=[3872 17603];
logdepthpicket=LOG(logdepthpicketvaluerange(1):logdepthpicketvaluerange(2),logDepthIndex);
logECGRpicket=LOG(logdepthpicketvaluerange(1):logdepthpicketvaluerange(2),logEcgrIndex) ;
logNphipicket=LOG(logdepthpicketvaluerange(1):logdepthpicketvaluerange(2),logNphiIndex) ;
logRhobpicket=LOG(logdepthpicketvaluerange(1):logdepthpicketvaluerange(2),logRhobIndex) ;
logDreshpicket=LOG(logdepthpicketvaluerange(1):logdepthpicketvaluerange(2),logDreshIndex) ;


subplot(1,3,1)
plot(logECGRpicket,logdepthpicket) %ECGR
axis ij
xlim([logEcgrUrXaxisRange(1,1) logEcgrUrXaxisRange(1,2)]);
ylim([logDepthRangepicket(1,1) logDepthRangepicket(1,2)]);
xlabel('Gamma Ray (gAPI)')
ylabel('Depth (meters)')
legend('Gamma Ray')

subplot(1,3,2)
plot(logNphipicket,logdepthpicket,'g')%neutron
xlim([logNphiXaxisRange(1,1) logNphiXaxisRange(1,2)])
ylim([logDepthRangepicket(1,1) logDepthRangepicket(1,2)])
axis ij
set(gca,'XDir','reverse')
xlabel('density & neutron')
format long

axRhob = gca;
ax1_pos = get(axRhob,'Position'); 
ax2 = axes('Position',[ax1_pos(1,1) ax1_pos(1,2) ax1_pos(1,3) ax1_pos(1,4)],...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
hold on
plot(logRhobpicket,logdepthpicket,'parent',ax2,'color','r')%densitclosy
axis ij
xlim([logRhobXaxisRange(1,1) logRhobXaxisRange(1,2)]);
ylim([logDepthRangepicket(1,1) logDepthRangepicket(1,2)]);
set(ax2,'XColor','r');
set(ax2,'YColor','r');
hold on 
axis ij 
legend('RHOB')

subplot(1,3,3)
semilogx(logDreshpicket, logdepthpicket, 'Color', 'r')%deep resistivity.
axis([logDreshXaxisRange(1,1) logDreshXaxisRange(1,2) logDepthRange(1,1) logDepthRange(1,2)])
ylim([logDepthRangepicket(1,1) logDepthRangepicket(1,2)]);
axis ij
format long
xlabel('resistivity')

%-------------------------------------------------------------------------------------------------------------------------------------
[xValue,yValue] = ginput(3);

text(xValue(1,1),yValue(1,1),'---------------Start Depth---------------', ...
                'HorizontalAlignment','center', ...
                'Color', [0.9 0 1], ...
                'FontSize',8);
text(xValue(2,1),yValue(2,1),'----------------End Depth----------------', ...
                'HorizontalAlignment','center', ...
                'Color', [0.9 0 1], ...
                'FontSize',8);
text(xValue(3,1),yValue(3,1),'------Index for Resistivity-Porosity-----', ...
                'HorizontalAlignment','center', ...
                'Color', [0.9 0 1], ...
                'FontSize',8);

for i=1:size(logdepthpicket,1) 
   if logdepthpicket(i,1)<min(yValue(1,1), yValue(2,1))
       startRange = i;
   end  
   if logdepthpicket(i,1)<max(yValue(1,1), yValue(2,1))
       endRange = i;
   end 
   if logdepthpicket(i,1)<yValue(3,1)
       initalResPorIndex = i;
   end 
end

logVshalepicket = (logECGRpicket-logSandLineHist)./(logShaleLineHist-logSandLineHist);
logVshalepicket(logVshalepicket<=0) = 0.001; 
logVshalepicket(logVshalepicket>1) = 0.999; 
porosityjustmadeforpicket=(logRhobpicket-(densityShale.*logVshalepicket+densitySand.*(1-logVshalepicket)))./(densityFluid-(densityShale.*logVshalepicket+densitySand.*(1-logVshalepicket)));


initalResistivity = logDreshpicket(initalResPorIndex, 1);

initalPorosity = porosityjustmadeforpicket(initalResPorIndex, 1);
%initalPorosity = porosityWithKerogen(initalResPorIndex, 1);

%geothermal gradient 
surfaceTemp = 14;
tempGradient = 0.025;
temperature = surfaceTemp+tempGradient.*logdepthpicket;
resistivityArps = initalResistivity*(temperature(initalResPorIndex)+21.5)./(temperature+21.5);

%-------------------------------------------------------------------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1])
picketlogDresh=logDreshpicket(startRange:endRange);
picketPorosity=porosityjustmadeforpicket(startRange:endRange);
loglog(picketlogDresh,picketPorosity,'o')
hold on

%-------------------------------------------------------------------------------------------------------------------------------------
dlg_title = 'Cementation Parameter';
prompt = {'Enter cementation Value', 'Inital Porosity', 'Inital Resistivity'};
input = {'2', num2str(initalPorosity), num2str(initalResistivity)};
check=1;
while(check)
    input = inputdlg(prompt,dlg_title,[1 69],input);
    if(isnan(str2double(input{1})))
         prompt{1,1}='Enter Valid cement Value';
    elseif(isnan(str2double(input{2})))
         prompt{1,2}='Enter Valid Porosity Value';
    elseif(isnan(str2double(input{3})))
         prompt{1,3}='Enter Valid Resistivity Value';
    else
        cementationExponentM = str2double(input{1});
        initalPorosity = str2double(input{2});
        initalResistivity = str2double(input{3});

        resistivityArps = initalResistivity*(temperature(initalResPorIndex)+21.5)./(temperature+21.5);
        waterResistivity = initalResistivity*(initalPorosity^cementationExponentM);

        %-------------------------------------------------------------------------------------------------------------------------------------
        close figure 4
        figure('units','normalized','outerposition',[0 0 1 1])
        picketlogDresh=logDreshpicket(startRange:endRange);
        picketPorosity=porosityjustmadeforpicket(startRange:endRange);
        loglog(picketlogDresh,picketPorosity,'o')
        hold on

        %-------------------------------------------------------------------------------------------------------------------------------------
        x = picketplotXaxisRange(1,1) : picketplotXaxisRange(1,2);
        x = log10(x); 
        y=(-1./cementationExponentM).*x + log(waterResistivity)/cementationExponentM;
        x = 10.^x;
        y = 10.^y;
        loglog(x,y,'r');

        format long
        % xlim([picketplotXaxisRange(1,1) picketplotXaxisRange(1,2)])
        % ylim([picketplotYaxisRange(1,1) picketplotYaxisRange(1,2)])
        legend('porosity & Resistivity')
        xlabel('porosity and resistivity')
        hold on 
        str = strcat('log y =  ',num2str(-1./cementationExponentM),' * log x + ',num2str(log(waterResistivity)/cementationExponentM));
        title(str);

        %-------------------------------------------------------------------------------------------------------------------------------------
        figure(1)
        subplot(1,10,9)
        hold off
        plot(c(:,9+numberOfMinerals+2),c(:,1),'o')
        xlim([saturationRange(1,1) saturationRange(1,2)])
        ylim([logDepthRange(1,1) logDepthRange(1,2)])
        xlabel('Saturation')
        set(gca,'YTick',[]);
        axis ij
        hold on
        saturationWithKerogen = zeros(size(logdepth));
        for i=1:length(logdepth)
            saturationWithKerogen(i,1)=((TortuosityFactorA*resistivityArps(i,1))/((porosityWithKerogen(i,1)^cementationExponentM)*logDresh(i,1)))^(1/saturationExponentN);            
        end
        plot(saturationWithKerogen,logdepth,'r')
        legend('Saturation','Sat_GRI')
        format long
    end    
end
%-------------------------------------------------------------------------------------------------------------------------------------
