clear all;
close all;
clc;

xrdFileName = 'xrd.xls';
griFileName = 'gri.xls';
logFileName = 'log.xls';

XRD = xlsread(xrdFileName);
GRI = xlsread(griFileName);
LOG = xlsread(logFileName);





% load('26_11.mat')
% load('dissertation_toc.mat')
% clc
% %input
% rhob=RHOB(16000:17670);
% DEPTH=DEPTHM(16000:17670);
% Filename = 'xrd.xls';
% xrd = xlsread(Filename,1);
% toc_xrd=xrd(:,3);
% depth_xrd=xrd(:,2);
% wk=1.1.*toc_xrd;
% w_nclay=xrd(:,4:10);%(4th column to 10 column has nonclay mineral)
% w_clay=xrd(:,11:14);%(11th column to 14 column has clay mineral)
% wk1=(1-wk./100);%to normalize 
% w_combined=cat(2,w_nclay,w_clay);
% %-------------------------------------------------------------------------
% %mineral densities
% qtz=2.65;kfeld=2.52;pfeld=2.66;              %silicates
% cal=2.71;dol=2.87;                              %carbonates
% pyr=4.99;marc=4.87;                             %heavies
% ill_smec=2.60;ill_mic=2.75;kaol=2.60;chl=2.76;  %clays
% densities=[qtz,kfeld,pfeld,cal,dol,pyr,marc,ill_smec,ill_mic,kaol,chl];
% %-------------------------------------------------------------------------
% %sort toc vs bulkd for same depth 
% b=[];
% for j= 1:43
%     for  k= 1:71

%         if depth_leco(k,1)==depth_blkd(j,1)
%             b(j,1)=depth_blkd(j,1);%depth
%             b(j,2)=depth_blkd(j,2);%bulkdensity from GRI file
%             b(j,3)=depth_leco(k,2);%TOC_leco    from SRA file
%         end
%     end
% end
% format long g
% %-------------------------------------------------------------------------
% %crossplot 
% blkd=b(4:43,2);
% toc=b(4:43,3); 

% plot(blkd,toc,'or')
% %hold on
% %-------------------------------------------------------------------------
% %upscale TOC
% toc_up=-9.556.*RHOB+26.25; %RHOB is from well log data
% %-------------------------------------------------------------------------
% %normalize by adding kerogen weight percentage to xrd mineral weight
% Wi_combined = zeros(size(w_clay));
% for i=1:length(wk1) 
%     for j=1:11    %change 11 and make it generalized
%    Wi_combined(i,j)=wk1(i,1).*w_combined(i,j);
%     end
% end
% %[sum(w_combined,2),sum(Wi_combined,2)]
% %matrix density
% rho_m_inverse=zeros(size(w_clay));
% for i=1:length(wk1)
%     for j=1:11
% rho_m_inverse(i,j)=(Wi_combined(i,j)./densities(1,j));
%     end 
% end
% rho_m=100./sum(rho_m_inverse,2); % matrix density(rho ma) with depth;(,2in sum is to add columns)
% %-------------------------------------------------------------------------
% %crossplot CORE TOC with (1/rho_m)------1/rho_m=a*toc_xrd+b
% rho_m1=1./rho_m;
% plot(toc_xrd,rho_m1,'o')
% subplot(1,3,1)
% plot(xrd(:,3),xrd(:,2),'o')
% axis ij
% subplot(1,3,2)
% plot(depth_leco(:,2),depth_leco(:,1),'*')
% axis ij
% subplot(1,3,3)
% plot(rho_m,xrd(:,2),'*')
% axis ij
% %-------------------------------------------------------------------------
% %POROSITY 
% % rho_m2 is matrix density without kerogen or orgqnic matter=b  %
% % rho_ke is kerogen density when toc=100% in above equation     %
% % toc has to be renamed as toc well log                         %
% % rho_fl is fluid density usually 0.75 to 0.9                   %
% %N=rho_m2-rhob.*(rho_m2.*toc./rho_ke-toc+1);
% %D=rho_m2-rho_fl+rho_m2.*toc.*(1-rho_m2./rho_ke);
% %phi=N./D;
% %-------------------------------------------------------------------------
% %SATURATION
% a=1;m=2;
% %Sw=(a.*Rw./(phi.^m.*Rt)).^0.5;     %define Rt from well log

% %-------------------------------------------------------------------------
% %PLOTS

% %plot(depth_blkd(4:43,2),depth_blkd(4:43,1),'or')
% %axis ij 
% %hold on 
% %plot(RHOB(16000:17670),DEPTHM(16000:17670))
% %axis ij
% subplot(1,2,1)
% plot(toc_up(16000:17670),DEPTHM(16000:17670),'g')
% axis ij 
% hold on 
% plot(depth_leco(:,2),depth_leco(:,1),'or')
% axis ij
% hold on
% %plot(toc_xrd,depth_xrd,'-ob')
% %axis ij
% subplot(1,2,2)
% %plot(GR,DEPTH_GR)
% axis ij
