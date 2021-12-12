%plot VAR vs. model impulse responses

pnstep=horizon;
time=[0:1:pnstep-1];

conflevel=0.95;
factor=-norminv((1-conflevel)/2);

global select_idx
select_idx=1:12;

get_VAR_resp; 

%VAR responses and standard deviations

%Monetary Shock

policyresponse = -IRFFF(1:pnstep,:);
policyresponseHIGH=policyresponse+factor*IRFFFSE(1:pnstep,:);
policyresponseLOW=policyresponse-factor*IRFFFSE(1:pnstep,:);

policyresponseLOW(:,[2,3])=policyresponseLOW(:,[2,3])*4;
policyresponseHIGH(:,[2,3])=policyresponseHIGH(:,[2,3])*4;
policyresponse(:,[2,3])=policyresponse(:,[2,3])*4;

%Neutral Technology Shock

zresponse = IRFz(1:pnstep,:);
zresponseHIGH=zresponse+factor*IRFzSE(1:pnstep,:);
zresponseLOW=zresponse-factor*IRFzSE(1:pnstep,:);

zresponseLOW(:,[2,3])=zresponseLOW(:,[2,3])*4;
zresponseHIGH(:,[2,3])=zresponseHIGH(:,[2,3])*4;
zresponse(:,[2,3])=zresponse(:,[2,3])*4;

%Investment Technology Shock

uresponse = IRFu(1:pnstep,:);
uresponseHIGH=uresponse+factor*IRFuSE(1:pnstep,:);
uresponseLOW=uresponse-factor*IRFuSE(1:pnstep,:);

uresponseLOW(:,[2,3])=uresponseLOW(:,[2,3])*4;
uresponseHIGH(:,[2,3])=uresponseHIGH(:,[2,3])*4;
uresponse(:,[2,3])=uresponse(:,[2,3])*4;

reorderidx=[1 10 2 3 5 6 7 9 8 4 12 11];

%Investment Technology Shock Plot 

figure;
legend_string={};

for plotindx=1:1:size(policyresponse,2);

    subplot(3,4,plotindx)
    
    grpyat = [ [(0:pnstep-1)'  uresponseLOW(:,reorderidx(plotindx))] ; [(pnstep-1:-1:0)'  uresponseHIGH(pnstep:-1:1,reorderidx(plotindx))]];
    patch(grpyat(:,1),grpyat(:,2),[0.5 0.5 0.5]); hold on
    
    plot(time, uresponse(1:pnstep,reorderidx(plotindx)),'k'); hold on
    plot(time,model_resp_itech_baseline(reorderidx(plotindx),1:pnstep),'g-'); hold on
    
    if plotindx==1,
            legend_string={['VAR ',num2str(conflevel*100),'%'],'VAR Mean',' Nash (D/w=0.37,re-estimated)'}; 
    end
    hold off

    title(data_varnames(reorderidx(plotindx)),'FontSize',12);
    axis tight

end

legend1=legend(legend_string, 'Orientation','horizontal');
set(legend1,'Position',[0.3 0.95 0.4 0.03]);

%Neutral Technology Shock Plot  

figure;
legend_string={}; 

for plotindx=1:1:size(policyresponse,2);

    subplot(3,4,plotindx)
    
    grpyat = [ [(0:pnstep-1)'  zresponseLOW(:,reorderidx(plotindx))] ; [(pnstep-1:-1:0)'  zresponseHIGH(pnstep:-1:1,reorderidx(plotindx))]];
    patch(grpyat(:,1),grpyat(:,2),[0.5 0.5 0.5]); hold on
    
    plot(time, zresponse(1:pnstep,reorderidx(plotindx)),'k'); hold on
    plot(time,model_resp_ntech_baseline(reorderidx(plotindx),1:pnstep),'g-'); hold on
    legend_string={['VAR ',num2str(conflevel*100),'%'],'VAR Mean',' Nash(D/w=0.37,re-estimated'}; 
    hold off

    title(data_varnames(reorderidx(plotindx)),'FontSize',12);
    axis tight

end

legend1=legend(legend_string, 'Orientation','horizontal');
set(legend1,'Position',[0.3 0.95 0.4 0.03]);

%Monetary Shock Plot  

figure;
legend_string={};

for plotindx=1:1:size(policyresponse,2);

    subplot(3,4,plotindx)
    
    grpyat = [ [(0:pnstep-1)'  policyresponseLOW(:,reorderidx(plotindx))] ; [(pnstep-1:-1:0)'  policyresponseHIGH(pnstep:-1:1,reorderidx(plotindx))]];
    patch(grpyat(:,1),grpyat(:,2),[0.5 0.5 0.5]); hold on
    
    plot(time, policyresponse(1:pnstep,reorderidx(plotindx)),'k'); hold on
    plot(time,model_resp_mon_baseline(reorderidx(plotindx),1:pnstep),'g-'); hold on
    
    legend_string={['VAR ',num2str(conflevel*100),'%'],'VAR Mean',' Nash (D/w=0.37,re-estimated'}; 
    hold off
    
    title(data_varnames(reorderidx(plotindx)),'FontSize',12);
    axis tight

end

legend1=legend(legend_string, 'Orientation','horizontal');
set(legend1,'Position',[0.3 0.95 0.4 0.03]);
