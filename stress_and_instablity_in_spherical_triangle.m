function stress_and_instablity_in_spherical_triangle(R,friction,plot_variable,varargin)
% plot_variable 指定绘图变量
%    'ns': normal stress
%    'ss': shear stress
%    'i' : instability
% 可选的输入变量：'isoline'，是否绘制等值线，该选项仅对 plot_variable 为 'i' 有效
%    'on' : 绘制等值线（默认）
%    'off': 不绘制

[x,y,ns,ss,instability]=NormalShearStress_and_instability(R,friction);

figure('color','w','position',[500,300,700,600])
hold on
axis equal
axis off
set(gca,'position',[0.15,0.2,0.6,0.6])
text(0.7,0.8,['R=',num2str(R)],'horizontalalignment','center','FontSize',14,'Fontname','TimesNewRoman')
h=colorbar('position',[0.8,0.25,0.025,0.4],'TickLength',0.02,'FontSize',10,'FontName','TimesNewRoman','LineWidth',0.8);

switch plot_variable
    case 'ns'
        pcolor(x,y,ns)
        caxis([-1,1])
        isocline % 绘制等倾线
        ylabel(h,'normal stress','FontSize',10,'Fontname','TimesNewRoman')

    case 'ss'
        pcolor(x,y,ss)
        caxis([0,1])
        isocline % 绘制等倾线
        ylabel(h,'shear stress','FontSize',10,'Fontname','TimesNewRoman')
        
    case 'i'
        pcolor(x,y,instability)
        caxis([0,1])
        isocline % 绘制等倾线
        ylabel(h,'instability','FontSize',10,'Fontname','TimesNewRoman')
        if isempty(varargin) || strcmp(varargin{2},'on')
            level_list=[0.7,0.8,0.9];
            contour(x,y,instability,level_list,'LineColor','k','ShowText','on','LineWidth',0.8)
        end
        
end

shading interp
add_coordinate



end




% ----------------- functions ------------------------------    

function [h,v]=kaverina(L,M,K,radius)
% 主轴坐标系下，断层面（法向）在球面三角形上的等面积投影
%   input:  L, K 分别为法向与压轴和张轴夹角余弦的绝对值
%   output: h, v 分别为投影后的横纵坐标
%   图的上，坐，右顶点分别为 sigma1,sigma2 和 sigma3 轴


if nargin==3
    radius=2;
end

phi=acos( 1/sqrt(3) * (L+M+K) );
r=radius*sin(phi/2);

if L==M && M==K && L==K
    nm=1;
else
    nm=sqrt( 2 * ( (K-M)^2 + (L-K)^2 + (M-L)^2 ) );
end

cos_theta=(2*L-M-K) / nm;
sin_theta=(K-M)*sqrt(3) / nm;


h=r*sin_theta;
v=r*cos_theta;


end


function isocline

np=100;
inter=10;

color=[240,255,255]/255;

% 法向与 sigma1 轴的等夹角线
d1=0:inter:90;
for i=1:length(d1)
    n1=cos(d1(i)*pi/180);
    n2_range=linspace(0,sqrt(1-n1^2),np);
    x=nan(1,np); y=nan(1,np);
    for j=1:np
        n2=n2_range(j);
        n3=sqrt(1-n1^2-n2^2);
        [x(j),y(j)]=kaverina(n1,n2,n3);
    end
    if d1(i)==90
        plot(real(x),real(y),'k','LineWidth',1.5)
    else
        plot(real(x),real(y),'-','color',color,'LineWidth',0.5)
    end
end

% 法向与 sigma2 轴的等夹角线
d2=0:inter:90;
for i=1:length(d2)
    n2=cos(d2(i)*pi/180);
    n3_range=linspace(0,sqrt(1-n2^2),np);
    x=nan(1,np); y=nan(1,np);
    for j=1:np
        n3=n3_range(j);
        n1=sqrt(1-n2^2-n3^2);
        [x(j),y(j)]=kaverina(n1,n2,n3);
    end
    if d2(i)==90
        plot(real(x),real(y),'k','LineWidth',1.5)
    else
        plot(real(x),real(y),'-','color',color,'LineWidth',0.5)
    end
end

% 法向与 sigma3 轴的等夹角线
d3=0:inter:90;
for i=1:length(d3)
    n3=cos(d3(i)*pi/180);
    n1_range=linspace(0,sqrt(1-n3^2),np);
    x=nan(1,np); y=nan(1,np);
    for j=1:np
        n1=n1_range(j);
        n2=sqrt(1-n1^2-n3^2);
        [x(j),y(j)]=kaverina(n1,n2,n3);
    end
    if d3(i)==90
        plot(real(x),real(y),'k','LineWidth',1.5)
    else
        plot(real(x),real(y),'-','color',color,'LineWidth',0.5)
    end
end


end


function [h,v,ns,ss,instability]=NormalShearStress_and_instability(R,friction)

np=200;
n1_list=linspace(0,1,np);
n2_max=sqrt(1-n1_list.^2);
n1=repmat(n1_list,np,1);
n2=nan(np,np);
for i=1:np
    n2(:,i)=reshape(linspace(0,n2_max(i),np), np, 1);
end

n3=sqrt(1-n1.^2-n2.^2);

h=nan(np,np);
v=nan(np,np);
for i=1:np
    for j=1:np
        [h(i,j),v(i,j)]=kaverina(n1(i,j),n2(i,j),n3(i,j));
    end
end
h=real(h);
v=real(v);

% 引自 Vavrycuk(2014) 的计算方法
ns=n1.^2+(1-2*R)*n2.^2-n3.^2;
ss=sqrt( ( n1.^2+(1-2*R)^2*n2.^2+n3.^2 ) - ns.^2 );
ss=real(ss);
instability=( ss-friction*(ns-1) ) / ( friction+sqrt(1+friction^2) );


end


function add_coordinate

inter=10;
n1=0;
d2=0:inter:90;
num=length(d2);

n2=cos( d2*pi/180 );
n3=sqrt(1-n1.^2-n2.^2);
x=nan(1,num);
y=nan(1,num);
for i=1:num
    [x(i),y(i)]=kaverina(n1,n2(i),n3(i));
end
text(x,y-0.05,[num2str(d2'),repmat('\circ',num,1)],'horizontalalignment','center','FontSize',11,'Fontname','TimesNewRoman')
% text(x(1)-0.08,y(1)-0.05,'\sigma_2','horizontalalignment','center','FontSize',18,'Fontname','TimesNewRoman')


n2=0;
d3=0:inter:90;
n3=cos( d3*pi/180 );
n1=sqrt(1-n2.^2-n3.^2);
x=nan(1,num);
y=nan(1,num);
for i=1:num
    [x(i),y(i)]=kaverina(n1(i),n2,n3(i));
end
text(x+0.07,y+0.02,[num2str(d3'),repmat('\circ',num,1)],'horizontalalignment','center','FontSize',11,'Fontname','TimesNewRoman')
% text(x(1)+0.12,y(1)-0.04,'\sigma_3','horizontalalignment','center','FontSize',18,'Fontname','TimesNewRoman')

n3=0;
d1=0:inter:90;
n1=cos( d1*pi/180 );
n2=sqrt(1-n1.^2-n3.^2);
x=nan(1,num);
y=nan(1,num);
for i=1:num
    [x(i),y(i)]=kaverina(n1(i),n2(i),n3);
end
text(x-0.05,y+0.02,[num2str(d1'),repmat('\circ',num,1)],'horizontalalignment','center','FontSize',11,'Fontname','TimesNewRoman')
% text(x(1),y(1)+0.12,'\sigma_1','horizontalalignment','center','FontSize',18,'Fontname','TimesNewRoman')

text(0,-0.75,'angle with \sigma_2','fontsize',11,'horizontalalignment','center','rotation',0,'Fontname','TimesNewRoman');
text(0.68,0.37,'angle with \sigma_3','fontsize',11,'horizontalalignment','center','rotation',-60,'Fontname','TimesNewRoman');
text(-0.67,0.37,'angle with \sigma_1','fontsize',11,'horizontalalignment','center','rotation',60,'Fontname','TimesNewRoman');


end





