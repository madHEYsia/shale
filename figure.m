
figure
porosityWithoutKerogen = (logRhob-(densityShale.*logVshale+densitySand.*(1-logVshale)))./(densityFluid-(densityShale.*logVshale+densitySand.*(1-logVshale)));
plot(porosityWithoutKerogen, logdepth,'k')
xlabel('phi-w/o-K')
xlim([porosityWithoutKeroXaxisRange(1,1) porosityWithoutKeroXaxisRange(1,2)])
ylim([logDepthRange(1,1) logDepthRange(1,2)])
axis ij
format long
hold on
plot(c(:,9+numberOfMinerals+1),c(:,1),'og')
set(gca,'YTick',[]);
format long

averageGrainDensity = 2.72;
densityMA = (100 - (1.1/densityKerogen).*tocPassey.*logRhob + (constantC/constantD + constantA/constantB) - logVclay/constantB - logVclay)*constantD;

numeratorPhi = densityMA - logRhob.*((densityMA.*tocPassey./100)/densityKerogen - tocPassey./100 + 1);
denominatorPhi = densityMA - densityFluid + densityFluid.*tocPassey./100.*(1 - densityMA./densityKerogen);
porosityWithKerogen = numeratorPhi./denominatorPhi;
hold on
plot(porosityWithKerogen, logdepth,'r')
legend('phi-w/o-K','GRI','phi-w-K')
format long
