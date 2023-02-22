clear all
%close all
set(0,'defaulttextInterpreter','latex') %latex axis labels
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'DefaultFigureColormap',flipud(hot));


Time=0:0.1:100000;

%%%%%%%%%%%%% colormap
mymap=[1,1,1];
for u=1:100
    mymap=[mymap;[1-0.01*u,1-0.01*u,1]];
end

N=71; % hypercube size
rc=60; % radius cutoff 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flux 2.2916
%%%%%%%%%%%%%% derive quasicrystal coordinates %%%%%%%%%%%%%%
a=[1 0; 1/sqrt(2) 1/sqrt(2); 0 1; -1/sqrt(2) 1/sqrt(2)];%%%%%%%%%%%%%%%% 4 by 2
as=[1 0; -1/sqrt(2) 1/sqrt(2); 0 -1; 1/sqrt(2) 1/sqrt(2)];%%%%%%%%%%%%%% 4 by 2
a0=1+sqrt(2);

ins=[a0 1; 1 a0; -1 a0; -a0 1];
ins=[ins; -a0 -1; -1 -a0; 1 -a0; a0 -1];
ins=ins/1.99999999999999;

dim=4;
vec=de2bi((1:N^dim)-1,dim,N);
vec(:,1)=vec(:,1)-35+0;
vec(:,2)=vec(:,2)-35+0;
vec(:,3)=vec(:,3)-35+0.5;
vec(:,4)=vec(:,4)-35+0;

x0=vec*a;
xs=vec*as;
%1.8000,32.2000,922.4000,750.4000
ind=inpolygon(xs(:,1),xs(:,2),ins(:,1),ins(:,2));

x=x0(ind,1); % final coordinates
y=x0(ind,2);

n=vec(ind,1);
m=vec(ind,2);
l=vec(ind,3);
k=vec(ind,4);

xx=xs(ind,1);
yy=xs(ind,2);

%%%%%%%%% cut radius %%%%%%%%%%%%%%%%
r=(x.^2+y.^2).^(1/2);
indr=(r<rc);

x=x(indr);
y=y(indr);
n=n(indr);
m=m(indr);
l=l(indr);
k=k(indr);

xx=xx(indr);
yy=yy(indr);


z=x+1i*y;
theta=abs(angle(z));%%%%%%%%%%%%%%%%%%%%%%%% Angular position.

X=sparse(diag(x));
Y=sparse(diag(y));
pos=exp(-min(abs(x),abs(y)));
pos=pos./norm(pos);
% POS=diag(pos);
% POS=sparse(POS);
%POS=sparse(abs(Y));


%%%%%%%%%%%%%build hamiltonian%%%%%%%%%%%%%

ind1i=[];ind2i=[];hvec=[];

for i=1:length(x)
    ind1=1:length(x);
    ind2=i*ones(1,length(x));  
    
    dx=x-x(i);
    dy=y-y(i); 
    dxy=(dx.^2+dy.^2).^(1/2); % calculate distance

    tempind=ind1(logical((abs(dxy-1)<10^-6)));
    ind1i=[ind1i ind1(tempind)];
    ind2i=[ind2i ind2(tempind)];
    hvec=[hvec 1+ind1(tempind)*0];
    
end

%%%%%%%%% test draw %%%%%%%%%%%%%%%%%
%figure('Position',[-3,34,1180,961]);
% figure
% xvec=reshape([x(ind1i) x(ind2i) x(ind2i)+NaN]',[3*length(x(ind2i)) 1]); 
% yvec=reshape([y(ind1i) y(ind2i) y(ind2i)+NaN]',[3*length(x(ind2i)) 1]);
% line(xvec,yvec,'Color',[0.00,0.00,1.00]);
% hold on
% scatter(x,y,30,'filled');
% %axis off
% axis equal
% hold off

%%%%%%%%%%%%%add flux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=sparse(ind2i,ind1i,hvec);
h(1:1+size(h,1):end) = 0;
h0=h;
[ind1,ind2,value]=find(h0);

[V,D]=eig(full(h));

fraction=0.675;
flux=fraction*pi;
value2=transpose(exp(1i*flux*(x(ind1)-x(ind2)).*(y(ind1)+y(ind2))))'.*value;
h=sparse(ind1,ind2,value2);
s=length(x);
delta=zeros([s 1]);
[M,I]=min(abs((x-43.9203).^2+(y-1.207).^2));
%[M,I]=min(abs((x+41.01).^2+(y-1.207).^2));
%[M,I]=min(abs((x-24.0208).^2+(y-26.935).^2));

delta(I)=1;
delta=sparse(delta);
delta0=delta;
%close all
%fig1=figure('Position',[1.8,50.6,867.2,724]);
fig2=figure('Position',[1.8,50.6,867.2,724]);
step=-sqrt(-1)*0.1;
% v = VideoWriter('dynamics.avi');
% open(v);

fidelity=[];


for j=0:10*Time(length(Time))
    fig2;
    fidelity=[fidelity,abs(delta0'*delta)];
    plot(0:0.1:Time(length(fidelity)),fidelity);
    xlabel('Time');
    ylabel('Fidelity');
    drawnow
%     if j==0
%         fig1;
%         scatter(x,y,100,abs(delta).^2,'filled');
%         colormap(mymap);
%         colorbar;
%         hold on
%         line(xvec,yvec,'Color',[0.8,0.8,0.8]);
%         hold off
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%         continue
%     end
    k1=step*h*delta;
    k2=step*h*(delta+k1./2);
    k3=step*h*(delta+k2./2);
    k4=step*h*(delta+k3);
    delta=delta+k1./6+k2./3+k3./3+k4./6;
    
    delta=delta/norm(delta);%%%%%%%%%%%%%% Normalization
    
%     fig1;
%     title("Time is "+Time(j))
%     scatter(x,y,100,abs(delta).^2,'filled');
%     colormap(mymap);
%     colorbar;
%     hold on
%     line(xvec,yvec,'Color',[0.8,0.8,0.8]);
%     hold off
%     frame=getframe(gcf);
%     writeVideo(v,frame);
end

%close(v);


% 
% hold on;angles = linspace(0, 2*pi, 500);
% radius = 14.4852;
% CenterX = 43.9203;
% CenterY = -3.6213;
% x = radius * cos(angles) + CenterX;
% y = radius * sin(angles) + CenterY;
% plot(x, y, 'r-', 'LineWidth', 2);
%9.94975,

    %print(strcat('C:\Users\junmo\OneDrive - kaist.ac.kr\NewQC\MP_theta\code','\',NUM,'\','0'),'-djpeg')
    %saveas(gcf,strcat('C:\Users\junmo\OneDrive - kaist.ac.kr\NewQC\MP_theta\code','\',NUM,'\0.fig'))
    


