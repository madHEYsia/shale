clear all;
close all;
clc;

%mineral densities  ------------------------------------------------------
qtz=2.65; kfeld=2.52; pfeld=2.66;              %silicates
cal=2.71; dol=2.87;                              %carbonates
pyr=4.99; marc=4.87;                             %heavies
ill_smec=2.60; ill_mic=2.75; kaol=2.60; chl=2.94;  %clays
densities=[qtz, kfeld, pfeld, cal, dol, pyr, marc, ill_smec, ill_mic, kaol,chl];
kerogen = 1.35;
%-------------------------------------------------------------------------	

xrdFileName = 'xrd.xls';
griFileName = 'gri.xlsx';
logFileName = 'log.xlsx';

logRange = [15998 17603];
logDepthIndex = 1;
logCgrIndex = 4;
logRhobIndex = 6;
logUranIndex = 7;
xrdTocIndex = 3;
xrdDepthIndex = 1;
xrdClayRange = [4 10];
xrdNonClayRange = [11 14];
griDepthIndex = 2;
griGrainDensityIndex = 7;
griBulkDensityIndex = 2;

grXrdPlotRange = [2 3];

XRD = xlsread(xrdFileName);
GRI = xlsread(griFileName);
LOG = xlsread(logFileName);

logRhob = LOG(logRange(1,1):logRange(1,2),logRhobIndex);

xrdToc = XRD(:,xrdTocIndex);
xrdDepth = XRD(:,xrdDepthIndex);

xrdClayWeightPercent = XRD(:,xrdClayRange(1,1):xrdClayRange(1,2));
xrdNonClayWeightPercent = XRD(:,xrdNonClayRange(1,1):xrdNonClayRange(1,2));

weightPercentCombine = cat(2,xrdClayWeightPercent, xrdNonClayWeightPercent);
numberOfMinerals = size(weightPercentCombine,2);

weightPercentKerogen = 1.1*xrdToc;
weightPercentKerogenNormFactor = (1-weightPercentKerogen/100);

%normalize by adding kerogen weight percentage to xrd mineral weight
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
    XRDGrainDensityWithKerogen(i,1) = 100/(temp + weightPercentKerogen(i,1)/kerogen);
end
% -------------------------------------------------------------------------

griDepth = GRI(:,griDepthIndex);
griGrainDensity = GRI(:,griGrainDensityIndex);
griBulkDensity = GRI(:,griBulkDensityIndex);

%grain density xrd vs gri calibration at same depth data
c=[];
xrdNonClayNormSum = sum(weightPercentsNormalized(:,1:size(xrdClayWeightPercent,2)),2);
xrdClayNormSum=sum(weightPercentsNormalized(:,(size(xrdClayWeightPercent,2)+1):size(weightPercentCombine,2)),2);
cIndex = 0;
for j= 1:length(GRI)       
    for  k= 1:length(XRD)
        if griDepth(j,1)==xrdDepth(k,1)
        	cIndex = cIndex+1;
            c(cIndex,1) = griDepth(j,1);
            c(cIndex,2) = griGrainDensity(j,1);
            c(cIndex,3) = XRDGrainDensityWithKerogen(k,1);
            c(cIndex,4) = xrdNonClayNormSum(k,1);
            c(cIndex,5) = xrdClayNormSum(k,1);
            c(cIndex,6) = weightPercentKerogen(k,1);
            c(cIndex,8) = griBulkDensity(j,1);
            c(cIndex,9:8+numberOfMinerals) = weightByDensityNormalized(k,:);
            c(cIndex,9+numberOfMinerals) = XRDGrainDensityWithoutKerogen(k,1);  
        end
    end 
end

c(:,7) = c(cIndex,6)/1.1;
% rho_m_inverse_common=c(:,11:21);
% blkd_gri_common=c(:,8);
% nclays_common=c(:,4);
% heavies_sum=sum(c(:,28:29),2);
% crystals_sum=sum(c(:,23:27),2);
% clays_sum=c(:,5);

% %converting wt% of minerals to volume%
% v=zeros(size(c(:,11:18)));
% for i=1:length(c(:,1))
%     for j=1:11
% v(i,j) =rho_m_inverse_common(i,j).*blkd_gri_common(i,1);% volume percentages of all minerals.

%     end
% end
% v1=v(:,8:11);% volume percentage of all clay 
% volume_TCLAY=sum(v1,2)./100;





format long g
hold on 
plot (c(:,2),c(:,3),'o')
xlim(grXrdPlotRange);
ylim(grXrdPlotRange);

hold on 
%y=x line
x=grXrdPlotRange(1,1):0.1:grXrdPlotRange(1,2);
y=grXrdPlotRange(1,1):0.1:grXrdPlotRange(1,2);
plot(x,y)

figure; % start different figure

logdepth = LOG(logRange,logDepthIndex);
CGR = LOG(logRange,logCgrIndex);
UR = LOG(logRange,logUranIndex);

subplot (1,10,1) 	%GR
plot(CGR,logdepth,'c') %ECGR
axis ij
% xlim([min(CGR,[],'omitnan') max(CGR,[],'omitnan')])
xlim([0 400])
ylim([min(logdepth) max(logdepth)])
xlabel('GR & UR')

ax1 = gca;
ax1_pos = get(ax1,'Position'); % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
hold on 
plot(UR,logdepth,'Parent',ax2,'Color','r')%UR
set(ax2,'XColor','r');
set(ax2,'YColor','r');
axis ij
% xlim([min(UR,[],'omitnan') max(UR,[],'omitnan')])
xlim([0 400])
ylim([min(logdepth) max(logdepth)])
format long
hold on 	
 
%  ax3 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% %line(toc_xrd,dep_xrd,'Parent',ax3,'Color','b','LineStyle',':')%toc
% ax3.XColor = 'b';
% ax3.YColor = 'b';
% axis ij
% xlim([0 10])
% ylim([a b])
% format long

 
%  hold on 
% plot(toc_xrd,dep_xrd,'oc')
% xlim([0 10])
% ylim([a b])
% axis ij 
% hold on 
% %------------------------------------------------------------------------

% subplot (1,10,2)%grain density

% plot(c(:,2),c(:,1),'ok')
% hold on 
% plot(c(:,3),c(:,1),'or')
% axis ij
% xlim([2 3])
% ylim([a b])
% xlabel('grain den')
% format long 
% hold on
% legend('GRI','XRD')
% %------------------------------------------------------------------------
% subplot (1,10,3)%neutron poristy and bulk density 

% plot(welllogs_new(j,3),welllogs_new(j,1),'g')%neutron
% xlim([-0.15 0.45])
% ylim([a b])
% axis ij
% set(gca,'XDir','reverse')
% xlabel('density and neutron')
% format long

% ax1 = gca;
% ax1_pos = ax1.Position;
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% hold on
% plot(zRHOB(j),welllogs_new(j,1),'parent',ax2,'color','r')%density
% axis ij
% xlim([1.8 2.9])
% ylim([a b])
% hold on 
% axis ij 
% legend('Nphi')
% %------------------------------------------------------------------------
% %sonic and resis passeys:-
% subplot (1,10,4)%sonic
% plot(z_DTC(j),welllogs_new(j,1),'g')
% xlim([0 150])
% ylim([a b])
% axis ij
% set(gca, 'XDir','reverse')
% format long

% ax1 = gca;
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% hold on 
% semilogx(DRESHOHMM(j,1),welllogs_new(j,1),'Color','r')%deep resistivity.
% axis([0.25 2500 a b])
% set(ax2,'xscale','log');
% axis ij
% format long
% legend('sonic')
% hold on 
% %------------------------------------------------------------------------

% vsh=(CGR-80)./(190-80);
% subplot (1,10,5)%vshale
% plot(vsh(j),welllogs_new(j,1))
% xlabel('vsh')
% xlim([-1 2])
% ylim([a b])
% axis ij
% hold on 

% vclay=0.45.*vsh;
% plot(vclay(j),welllogs_new(j,1),'r')
% hold on 
% plot(volume_TCLAY,c(:,1),'ok')
% legend('Vsh','Vcl','XRDclay')
% %------------------------------------------------------------------------
% rho_sh=2.71; rho_sand=2.65;
% phi=(zRHOB(j)-(rho_sh.*vsh(j)+rho_sand.*(1-vsh(j))))./(1-(rho_sh.*vsh(j)+rho_sand.*(1-vsh(j))));
% subplot (1,10,6)%porosity
% plot(phi,welllogs_new(j,1),'k')
% xlabel('phi')
% xlim([-0.5 0.5])
% ylim([a b])
% axis ij
% format long
% hold on
% plot(POR_gri,DEPTH_gri,'og')
% xlim([-0.5 0.5])
% ylim([a b])
% axis ij
% format long
% legend('phi_Vsh','GRI')

% %-------    ----   ------       --------     ------    -------   ------   --------   --------
% hold on 
% grain_withouttoc=c(:,22);

% grain_withouttoc(grain_withouttoc==0)=NaN;
% grain_withouttoc(grain_withouttoc==0)=NaN;
% plot(grain_withouttoc,volume_TCLAY,'o')
% xlim([2.4 3])
% ylim([0 0.55])



% %{

% N=rho_m_star-rhob.*(rho_m_star.*toc./rho_ke-toc+1);
% D=rho_m_star-rho_fl+rho_m_star.*toc.*(1-rho_m_star./rho_ke);
% phi=N./D;
% plot()
% %}
% %------------------------------------------------------------------------
%  %Passeys TOC
%  DRESH=DRESHOHMM(j,1);
%  DRESH_base=35; %from histogram or graph
%  DTC=z_DTC(j,1);
%  DTC_base=60; %from histogram or graph
%  LOM=12;SF1=70;SF2=0.75;
% DLogR=log(DRESH./DRESH_base)+0.02.*(DTC-DTC_base);
% toc_passey=SF1.*DLogR.*10.^(0.297-0.1688.*LOM)+SF2;

% subplot (1,10,7)%TOC passeys and XRD;
% plot(toc_passey,welllogs_new(j,1),'r')
% xlim([0 10])
% ylim([a b])
% hold on 
% plot(c(:,7),c(:,1),'ok')
% xlim([-1 10])
% ylim([a b])
% axis ij 
% hold on 
% xlabel('TOC')
% legend('TOC_Passey','TOC_XRD')

% %------------------------------------------------------------------------

% subplot (1,10,8)%bulk density
% plot(zRHOB(j,1),welllogs_new(j,1),'k')
% xlabel('bulk density')
% xlim([2 3])
% ylim([a b])
% axis ij
% hold on  
% plot(c(:,8),c(:,1),'og')
% xlim([2 3])
% ylim([a b])
% axis ij
% format long
% legend('Log blkd','GRI')
% %------------------------------------------------------------------------
% a=1;m=2;
% subplot (1,10,9)%saturation from logs and core (gri)
% %{
% Sw=(a.*Rw./(phi.^m.*DRESHOHMM)).^0.5;


% plot(Sw(j,1),welllogs_new(j,1),'k')

% xlim([0 1])
% ylim([a b])
% axis ij
% hold on  
% %}

% plot(SW_gri,DEPTH_gri,'ok')
% xlim([0 1])
% ylim([2800 3100])
% axis ij
% xlabel('saturation')
% format long
% hold on 

% %------------------------------------------------------------------------
%  %display clay non clays and kerogen

    
% subplot (1,10,10)
% plot(nclay1,d1,'y')
% hold on 
% plot(clay1,d1,'k')
% hold on 
% plot(kerogen1,d1,'r')
% hold on 
% axis ij 
% xlim([0 100])
% ylim([2800 3100])
% legend('NClay','Clay','Kerogen')

% %------------------------------------------------------------------------
% %rhog_Vsh=rho_sh.*vsh(j)+rho_sand.*(1-vsh(j));
% %plot(rhog_Vsh,)