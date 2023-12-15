clear all, close all, clc 

Folder   = '\\triton.meca.polymtl.ca\usagers\pemeh\profiles\Desktop\PhD\Rotatd Triangle\Ziada\XP2.08-0.0949-Re1870\SVD'; 

files = dir('*.csv');

for i=1:length(files);
    [folder, baseFileName, extension] = fileparts(files(i).name);
    [onlyFileNames{i}] = baseFileName;
    a=str2double(onlyFileNames);
    a=sort(a); 
end
for i=1:length(files);
    files(i).name=a(1,i);
    files(i).name=num2str(files(i).name);
    files(i).name = [files(i).name '.csv'];
end
%%

n=100;%number_of_sample_points
m=length(files);%number_of_snapshots
o=2;%number_of_modes 


X=zeros(n,length(files));
% Y=zeros(n,length(files));
Velocity=zeros(n,length(files));
% Velocity_v=zeros(n,length(files));
% VOL=zeros(n,length(files));
fu=zeros(n,length(files));
PSD=zeros(length(files),o);
F=zeros(length(files),o);

for i=1:length(files);
    data = csvread(files(i).name,5,0);
    X(:,i)=data(:,1);%x-direction
%     Y(:,i)=data(:,2);%y-direction
    Velocity(:,i)=data(:,2);%u-velocity
%     Velocity_v(:,i)=data(:,5);%v-velocity
%     VOL(:,i)=data(:,6);%vorticity
end


%%

% d=0.0127;
% [M N]=size(Velocity_v);
% Vavg=mean(Velocity_v,2);
% Vavg=Vavg*ones(1,N);
% fv=Velocity_v-Vavg;%fluctuation_velocity_v


d=0.0127;%mm
XP=1.6164;%mm
L1=0.8082*0.0254;%mm
T=0.7*0.0254;%mm


[M, N]=size(Velocity);
Uavg=mean(Velocity,2);
Uavg=Uavg*ones(1,N);
fu=Velocity-Uavg;%fluctuation_velocity_v
[U,S,V] = svd(fu);


figure
semilogy(diag(S)/sum(diag(S)),'*');
% title('Singular values of matrix')
xlabel('Mode')
ylabel('Energy Per Mode')
grid on

% Calculate the cumulative energy distribution for the modes
cumulative_energy = cumsum(diag(S).^2) / sum(diag(S).^2);

% Plot the cumulative energy distribution
figure;
plot(cumulative_energy, 'b*-');
title('Cumulative Energy Distribution');
xlabel('Mode');
ylabel('Cumulative Energy');
grid on;
%%
dt=0.005;
t=0:dt:15;

for i=1:4
    l1=length(V(:,i));
    freq=1/(dt*l1)*(0:l1);
    L=1:floor(l1/2);
    F(:,i)=fft(V(:,i));
    PSD(:,i)=F(:,i).*conj(F(:,i))/l1;
end

figure
 for i=1:4

    kk=(i)*ones(length(L),1);
    ROWS(:,1)=PSD(L,1);
    ROWS(:,2)=PSD(L,2);
    ROWS(:,3)=PSD(L,3);
    ROWS(:,4)=PSD(L,4);

    subplot(4,1,1)
    plot(freq(L),ROWS(L,1),'k')
    title('Mode 1')
    xlabel('Frequency [Hz]')
    ylabel({'Power Spectrum ', '(normalized)'})
    xlim([0 50])
    subplot(4,1,2)
    plot(freq(L),ROWS(L,2),'k')
    title('Mode 2')
    xlabel('Frequency [Hz]')
    ylabel({'Power Spectrum ', '(normalized)'})
    xlim([0 50])
    hold on
    subplot(4,1,3)
    plot(freq(L),ROWS(L,3),'k')
    title('Mode 3')
    xlabel('Frequency [Hz]')
    ylabel({'Power Spectrum ', '(normalized)'})
    xlim([0 50])
    hold on
    subplot(4,1,4)
    plot(freq(L),ROWS(L,4),'k')
    title('Mode 4')
    xlabel('Frequency [Hz]')
    ylabel({'Power Spectrum ', '(normalized)'})
    xlim([0 50])
    hold on
    end
%% geometry
% a1=(X(1,1)+0.0143)/XP;
% a2=(X(n,1)-0.0143)/XP;
% b1=(Y(1,1)-0.0111)/XP;
p=8;%number of tubes
XC=zeros(1,p);
YC=zeros(1,p);
for i=1:p;
    XC(i)=((T/2)*(2*i-1))/L1;
    if  i==5 | i==7;
        YC(i)=(L1/2)/L1;
    end
end

figure
for i=1:o
    subplot(1,o,i)
    plot(X(:,i)/L1,U(:,i),'k')
      ylim( [-1.5, 1])% check manually
      xlim([0 6.91])
        for j=1:length(XC)
            circle(XC(j),-1.5+YC(j),d/(2*L1));%put the min of line 117 as y
        end
    title(sprintf('POD Mode %i',i))
    xlabel('x/L')
    ylabel('Spatial eigen vector')
    grid on
    hold on
end
function circles = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
x_circle = r * cos(th) + x;
y_circle = r * sin(th) + y;
circles = plot(x_circle, y_circle,'k');
daspect([1,1,1])
hold off
end



   














% XC,min(U(:,i)),'ko','MarkerSize',48
%             set(gca,'xtick',[0:1:8]);
%              ylim([0 3])
%     set(gca,'xtick',[0:1:8]);