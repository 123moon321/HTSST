clc;clear;close all;
%%
N = 1000; 
fs = 100;
t =(0:N-1)/fs; 
w = (0:N/2)*fs/N;
f = (0:N/2)*fs/N;
%%
Amp1=exp(0.004*w);
Ph1=1.5*w+0.06/2*w.^2-0.0006/3*w.^3;
GD1=1.5+0.06*w-0.0006*w.^2;
x1=Amp1.*exp(-1i*2*pi*Ph1);    
x1(end) = -abs(x1(end));
x1 = [x1 conj(fliplr(x1(2:end-1)))];
y1=ifft(x1);

Amp2=exp(0.008*w);
Ph2=sin(28*pi*(w/100).^2)+6*w;
GD2=(56/100^2*pi*w).*cos(28*pi*(w/100).^2)+6;
x2=Amp2.*exp(-1i*2*pi*Ph2);
x2(end) = -abs(x2(end));
x2 = [x2 conj(fliplr(x2(2:end-1)))];
y2=ifft(x2);

ys = y1+y2;
sig = ys;
time = t;
fre = w;
%%
Q1=100;
dt=0.01;
timeinterval=dt*100;
sig=wavatten(sig',Q1,timeinterval,dt,1);
%%
s=0.2;
method = struct( 'type' , 'STFT' );
tfr = HTSST(sig , fs , s , method);

method = struct( 'type' , 'TSST' , 'order' , 1 , 'iteration' , 1 );
TSST = HTSST(sig', fs , s , method);

WindowOpt1 = struct('type','gauss','s',0.2);    %0.10
Parameter1 = struct('L',round(N/2)+1,'fmin',0,'fmax',fs/2);
[TSST2] = GHST(sig' , fs,  WindowOpt1, Parameter1, '2Ord');

method = struct( 'type' , 'TSST' , 'order' , 3 , 'iteration' , 1 );
TSST3 = HTSST(sig', fs , s , method);
%% 
minnum_tfr=min(min(tfr));maxnum_tfr=max(max(tfr));tfr=(tfr-minnum_tfr)./(maxnum_tfr-minnum_tfr);
minnum_TSST=min(min(TSST));maxnum_TSST=max(max(TSST));TSST=(TSST-minnum_TSST)./(maxnum_TSST-minnum_TSST);
minnum_TSST2=min(min(TSST2));maxnum_TSST2=max(max(TSST2));TSST2=(TSST2-minnum_TSST2)./(maxnum_TSST2-minnum_TSST2);
minnum_TSST3=min(min(TSST3));maxnum_TSST3=max(max(TSST3));TSST3=(TSST3-minnum_TSST3)./(maxnum_TSST3-minnum_TSST3);
%%
mymap=jet;

label=strings(9,1);
for k=1:9
    label(k)=['(',char(96+k),')'];
end

fig1 = figure('color',[1 1 1],'Units','characters','Position',[3 3 200 40]);
subplot('position',[0.05 0.10 0.38 0.85]);
imagesc(linspace(0,10,N),w,abs(tfr));axis xy;
xlabel('Time (s)','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ylabel('Frequency (Hz)','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
title(label(1),'Position',[-0.95 49],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.44 0.10 0.015 0.85]);colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
set(gca,'xtick',0:2:10) ;
set(gca,'ytick',0:10:50);
axis([0 10 0 50]);
m1=1.8;n1=18;width1=1.5;height1=5;
m2=5;n2=30;width2=1.5;height2=7;
rectangle('Position',[m1,n1,width1,height1],'LineWidth',1.5,'EdgeColor','y');
rectangle('Position',[m2,n2,width2,height2],'LineWidth',1.5,'EdgeColor','r');
%%
ha = subplot('position',[0.51 0.77 0.18 0.18]);
imagesc(t(fs*m1+1:fs*(m1+width1)+1),w(N/fs*n1+1:N/fs*(n1+height1)+1),abs(tfr(N/fs*n1+1:N/fs*(n1+height1)+1,fs*m1+1:fs*(m1+width1)+1)));
hold on;
plot(GD1(N/fs*n1+1:N/fs*(n1+height1)+1),w(N/fs*n1+1:N/fs*(n1+height1)+1),'y','linewidth',1.0);
title(label(2),'position',[1.69 18.7],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.695 0.77 0.011 0.18]);
colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ha.YColor = 'yellow';ha.XColor = 'yellow';

ha = subplot('position',[0.765 0.77 0.18 0.18]);
imagesc(t(fs*m2+1:fs*(m2+width2)+1),w(N/fs*n2+1:N/fs*(n2+height2)+1),abs(tfr(N/fs*n2+1:N/fs*(n2+height2)+1,fs*m2+1:fs*(m2+width2)+1)));
hold on;
plot(GD2(N/fs*n2+1:N/fs*(n2+height2)+1),w(N/fs*n2+1:N/fs*(n2+height2)+1),'r','linewidth',1.0);
title(label(6),'position',[4.89 28],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.95 0.77 0.011 0.18]);colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ha.YColor = 'red';ha.XColor = 'red';
%%
ha = subplot('position',[0.51 0.55 0.18 0.18]);
imagesc(t(fs*m1+1:fs*(m1+width1)+1),w(N/fs*n1+1:N/fs*(n1+height1)+1),abs(TSST(N/fs*n1+1:N/fs*(n1+height1)+1,fs*m1+1:fs*(m1+width1)+1)));
hold on;
plot(GD1(N/fs*n1+1:N/fs*(n1+height1)+1),w(N/fs*n1+1:N/fs*(n1+height1)+1),'y','linewidth',1.0);
title(label(3),'position',[1.69 18.7],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.695 0.55 0.011 0.18]);colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ha.YColor = 'yellow';ha.XColor = 'yellow';

ha = subplot('position',[0.765 0.55 0.18 0.18]);
imagesc(t(fs*m2+1:fs*(m2+width2)+1),w(N/fs*n2+1:N/fs*(n2+height2)+1),abs(TSST(N/fs*n2+1:N/fs*(n2+height2)+1,fs*m2+1:fs*(m2+width2)+1)));
hold on;
plot(GD2(N/fs*n2+1:N/fs*(n2+height2)+1),w(N/fs*n2+1:N/fs*(n2+height2)+1),'r','linewidth',1.0);
title(label(7),'position',[4.89 28],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.95 0.55 0.011 0.18]);colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ha.YColor = 'red';ha.XColor = 'red';
%%
ha = subplot('position',[0.51 0.33 0.18 0.18]);
imagesc(t(fs*m1+1:fs*(m1+width1)+1),w(N/fs*n1+1:N/fs*(n1+height1)+1),abs(TSST2(N/fs*n1+1:N/fs*(n1+height1)+1,fs*m1+1:fs*(m1+width1)+1)));
hold on;
plot(GD1(N/fs*n1+1:N/fs*(n1+height1)+1),w(N/fs*n1+1:N/fs*(n1+height1)+1),'y','linewidth',1.0);
title(label(4),'position',[1.69 18.7],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.695 0.33 0.011 0.18]);colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ha.YColor = 'yellow';ha.XColor = 'yellow';

ha = subplot('position',[0.765 0.33 0.18 0.18]);
imagesc(t(fs*m2+1:fs*(m2+width2)+1),w(N/fs*n2+1:N/fs*(n2+height2)+1),abs(TSST2(N/fs*n2+1:N/fs*(n2+height2)+1,fs*m2+1:fs*(m2+width2)+1)));
hold on;
plot(GD2(N/fs*n2+1:N/fs*(n2+height2)+1),w(N/fs*n2+1:N/fs*(n2+height2)+1),'r','linewidth',1.0);
title(label(8),'position',[4.89 28],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.95 0.33 0.011 0.18]);colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ha.YColor = 'red';ha.XColor = 'red';
%%
ha = subplot('position',[0.51 0.11 0.18 0.18]);
imagesc(t(fs*m1+1:fs*(m1+width1)+1),w(N/fs*n1+1:N/fs*(n1+height1)+1),abs(TSST3(N/fs*n1+1:N/fs*(n1+height1)+1,fs*m1+1:fs*(m1+width1)+1)));
hold on;
plot(GD1(N/fs*n1+1:N/fs*(n1+height1)+1),w(N/fs*n1+1:N/fs*(n1+height1)+1),'y','linewidth',1.0);
title(label(5),'position',[1.69 18.7],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.695 0.11 0.011 0.18]);colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ha.YColor = 'yellow';ha.XColor = 'yellow';

ha = subplot('position',[0.765 0.11 0.18 0.18]);
imagesc(t(fs*m2+1:fs*(m2+width2)+1),w(N/fs*n2+1:N/fs*(n2+height2)+1),abs(TSST3(N/fs*n2+1:N/fs*(n2+height2)+1,fs*m2+1:fs*(m2+width2)+1)));
hold on;
plot(GD2(N/fs*n2+1:N/fs*(n2+height2)+1),w(N/fs*n2+1:N/fs*(n2+height2)+1),'r','linewidth',1.0);
title(label(9),'position',[4.89 28],'units','normalized','FontName','Times New Roman','FontSize',15,'FontWeight','bold');
colorbar('position',[0.95 0.11 0.011 0.18]);colormap(mymap);
set(gca,'ydir','normal','xtick',[],'ytick',[],'linewidth',2,'FontName','Times New Roman','FontSize',15,'FontWeight','bold');
ha.YColor = 'red';ha.XColor = 'red';

r1 = renyi(abs(tfr));
r2 = renyi(abs(TSST));
r3 = renyi(abs(TSST2));
r4 = renyi(abs(TSST3));