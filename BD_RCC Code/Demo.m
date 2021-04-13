%%
clc
clear all
close all
addpath('./minFunc');
Fs=10000;
N=Fs*1;
randn('seed',35);
N=5000;
fts=[0:1:N-1]*Fs/N;
defect = struct('type', 'outer', 'num', 1, 'interval', 5/19);
fgg=100;
param.RE = 1e-4;
param.iter = 100;
tep=floor(fgg./(Fs./N));
alpha=tep:tep:N-1;
alpha=alpha(1:3);
alpha=[alpha,alpha-1,alpha+1,alpha-2,alpha+2,alpha-3,alpha+3];

for i=1:101
    mparameters = struct('Fs', 10000, 'N', N, 'f_BPFO', fgg, 'f_RF',5.29,...
        'RSF', 2000, 'cyclostationary', 0, 'SDSF', 1/100,...
        'D', 1, 'maxtheta', 0.5*pi, 'qn', 10/9, 'SNR', -0.2*(i-1), 'PLOT', 0,'Beta',800);
    [odata, x, t, mparameters,~] = BearingSimulation(defect, mparameters);
    sx=odata;
    %%
    [optW2,rec2(i,:),fuv2(i)]=BD_RCC(sx,100,fgg,Fs);
    %%
    [h_est,s0(i,:),kappa(i),W,count,err] = MaxCycloBD(sx,100,alpha,Fs,param,2);
    %%
    [rec3(i,:), f_final, ckIter] = mckd(sx,100,400,floor(Fs/fgg),1,0);
    ck(i)=ckIter(end);
    %%
    [MKurt(i) f temp] = momeda(sx,100,ones(1,1),floor(Fs/fgg),0);
     temp2=zeros(1,4901);
     temp2(1:length(temp))=temp;
     rec4(i,:)=temp2;
    %%
    [optW,rec5(i,:),valiter]=min_lplq_Corr(sx,100,floor(Fs/fgg),1,2);
    corr(i)=valiter.trace.fval(end);
    %%
    blp1(i,:)=abs(fft(abs((rec2(i,:))).^2))/length(rec2(i,:));
    blp2(i,:)=abs(fft(abs((s0(i,:))).^2))/length(s0(i,:));
    blp3(i,:)=abs(fft(abs((rec3(i,:))).^2))/length(rec3(i,:));
    blp4(i,:)=abs(fft(abs((rec4(i,:))).^2))/length(rec4(i,:));
    blp5(i,:)=abs(fft(abs((rec5(i,:))).^2))/length(rec5(i,:));
    i
end
%%
fres=10000/4901;
close all
blp11=blp1;blp22=blp2;blp33=blp3;blp44=blp4;blp55=blp5;
blp11(:,2:end)=mapminmax(blp1(:,2:end),0,1);
blp22(:,2:end)=mapminmax(blp2(:,2:end),0,1);
blp33(:,2:end)=mapminmax(blp3(:,2:end),0,1);
blp44(:,2:end)=mapminmax(blp4(:,2:end),0,1);
blp55(:,2:end)=mapminmax(blp5(:,2:end),0,1);

ACC1=blp11(:,2:400);
ACC2=blp22(:,2:400);
ACC3=blp33(:,2:400);
ACC4=blp44(:,2:400);
ACC5=blp55(:,2:400);

ACC1=[zeros(size(ACC1(:,1))),ACC1];
ACC2=[zeros(size(ACC2(:,1))),ACC2];
ACC3=[zeros(size(ACC3(:,1))),ACC3];
ACC4=[zeros(size(ACC4(:,1))),ACC4];
ACC5=[zeros(size(ACC5(:,1))),ACC5];

aa=300;
bb=300;
xma=400;
figure
colormap(gray)
imagesc([0:1:400]*fres,[0:-0.2:-20],ACC1)
xlim([0,xma])
xticks([0:100:xma])
ylabel('SNR (dB)','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
colormap(gray)
imagesc([0:1:400]*fres,[0:-0.2:-20],ACC2)
% xticks([2:800])
ylabel('SNR (dB)','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
xlim([0,xma])

xticks([0:100:xma])

figure
colormap(gray)
imagesc([0:1:400]*fres,[0:-0.2:-20],ACC3)
% xticks([2:800])
ylabel('SNR (dB)','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
xlim([0,xma])
xticks([0:100:xma])

figure
colormap(gray)
imagesc([0:1:400]*fres,[0:-0.2:-20],ACC4)
% xticks([2:800])
ylabel('SNR (dB)','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
xlim([0,xma])
xticks([0:100:xma])

figure
colormap(gray)
imagesc([0:1:400]*fres,[0:-0.2:-20],ACC5)
% xticks([2:800])
ylabel('SNR (dB)','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
xlim([0,xma])
xticks([0:100:xma])
%%
close all
figure
plot([0:-0.2:-20],-fuv2,'LineWidth',1)
ylabel('ESHNR','fontsize',12)
xlabel('SNR (dB)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot([0:-0.2:-20],kappa,'LineWidth',1)
ylabel('ICS2','fontsize',12)
xlabel('SNR (dB)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot([0:-0.2:-20],ck,'LineWidth',1)
ylabel('CK','fontsize',12)
xlabel('SNR (dB)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot([0:-0.2:-20],MKurt,'LineWidth',1)
ylabel('MKurt','fontsize',12)
xlabel('SNR (dB)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot([0:-0.2:-20],corr.*size(rec5,2)^(1/2-1),'LineWidth',1)
ylabel('CG-L1/L2','fontsize',12)
xlabel('SNR (dB)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

%%
close all
indx=ceil(13/0.2);  %-13dB
aa=300;
bb=200;
figure
plot([1:2:200*2],ACC1(indx,1:200),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot([1:2:200*2],ACC2(indx,1:200),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot([1:2:200*2],ACC3(indx,1:200),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot([1:2:200*2],ACC4(indx,1:200),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);

figure
plot([1:2:200*2],ACC5(indx,1:200),'LineWidth',1)
ylabel('Amplitude','fontsize',12)
xlabel('Frequency (Hz)','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,aa,bb]);
