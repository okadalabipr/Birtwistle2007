using PyPlot;

rc("figure",figsize = (12,15));
rc("font",size = 10);
rc("axes",linewidth = 1);
rc("lines",linewidth = 2);

subplots_adjust(wspace=0.4, hspace=0.5);

subplot(4,2,1)
plot(t,ERK_act[:,1]./589.5,"-",label="EGF=0 nM");
plot(t,ERK_act[:,4]./589.5,"--",label="EGF=0.5 nM");
plot(t,ERK_act[:,7]./589.5,":",label="EGF=10 nM");
title("EGF increasing; HRG=0.5 nM");
xlim(0,1900);
xticks([0,500,1000,1500]);
xlabel("Time (s)");
ylim(0,1.5);
yticks([0,0.5,1,1.5]);
ylabel("Normalized ERK*");
legend(loc="upper right",fontsize=8);

subplot(4,2,2)
plot(t,ERK_act[:,2]./589.5,"-",label="EGF=0 nM");
plot(t,ERK_act[:,5]./589.5,"--",label="EGF=0.5 nM");
plot(t,ERK_act[:,8]./589.5,":",label="EGF=10 nM");
title("EGF increasing; HRG=10 nM");
xlim(0,1900);
xticks([0,500,1000,1500]);
xlabel("Time (s)");
ylim(0,1.5);
yticks([0,0.5,1,1.5]);
ylabel("Normalized ERK*");
legend(loc="upper right",fontsize=8);

subplot(4,2,3)
plot(t,ERK_act[:,3]./589.5,"-",label="HRG=0 nM");
plot(t,ERK_act[:,4]./589.5,"--",label="HRG=0.5 nM");
plot(t,ERK_act[:,5]./589.5,":",label="HRG=10 nM");
title("EGF=0.5 nM; HRG increasing");
xlim(0,1900);
xticks([0,500,1000,1500]);
xlabel("Time (s)");
ylim(0,1.5);
yticks([0,0.5,1,1.5]);
ylabel("Normalized ERK*");
legend(loc="upper right",fontsize=8);

subplot(4,2,4)
plot(t,ERK_act[:,6]./589.5,"-",label="HRG=0 nM");
plot(t,ERK_act[:,7]./589.5,"--",label="HRG=0.5 nM");
plot(t,ERK_act[:,8]./589.5,":",label="HRG=10 nM");
title("EGF=10 nM; HRG increasing");
xlim(0,1900);
xticks([0,500,1000,1500]);
xlabel("Time (s)");
ylim(0,1.5);
yticks([0,0.5,1,1.5]);
ylabel("Normalized ERK*");
legend(loc="upper right",fontsize=8);

subplot(4,2,5)
plot(t,Akt_act[:,1]./18.8,"-",label="EGF=0 nM");
plot(t,Akt_act[:,4]./18.8,"--",label="EGF=0.5 nM");
plot(t,Akt_act[:,7]./18.8,":",label="EGF=10 nM");
title("EGF increasing; HRG=0.5 nM");
xlim(0,1900);
xticks([0,500,1000,1500]);
xlabel("Time (s)");
ylim(0,1.5);
yticks([0,0.5,1,1.5]);
ylabel("Normalized Akt*");
legend(loc="upper right",fontsize=8);

subplot(4,2,6)
plot(t,Akt_act[:,2]./18.8,"-",label="EGF=0 nM");
plot(t,Akt_act[:,5]./18.8,"--",label="EGF=0.5 nM");
plot(t,Akt_act[:,8]./18.8,":",label="EGF=10 nM");
title("EGF increasing; HRG=10 nM");
xlim(0,1900);
xticks([0,500,1000,1500]);
xlabel("Time (s)");
ylim(0,1.5);
yticks([0,0.5,1,1.5]);
ylabel("Normalized Akt*");
legend(loc="upper right",fontsize=8);

subplot(4,2,7)
plot(t,Akt_act[:,3]./18.8,"-",label="HRG=0 nM");
plot(t,Akt_act[:,4]./18.8,"--",label="HRG=0.5 nM");
plot(t,Akt_act[:,5]./18.8,":",label="HRG=10 nM");
title("EGF=0.5 nM; HRG increasing");
xlim(0,1900);
xticks([0,500,1000,1500]);
xlabel("Time (s)");
ylim(0,1.5);
yticks([0,0.5,1,1.5]);
ylabel("Normalized Akt*");
legend(loc="upper right",fontsize=8);

subplot(4,2,8)
plot(t,Akt_act[:,6]./18.8,"-",label="HRG=0 nM");
plot(t,Akt_act[:,7]./18.8,"--",label="HRG=0.5 nM");
plot(t,Akt_act[:,8]./18.8,":",label="HRG=10 nM");
title("EGF=10 nM; HRG increasing");
xlim(0,1900);
xticks([0,500,1000,1500]);
xlabel("Time (s)");
ylim(0,1.5);
yticks([0,0.5,1,1.5]);
ylabel("Normalized Akt*");
legend(loc="upper right",fontsize=8);