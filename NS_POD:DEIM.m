% 完全降维模型
% 用F（rhs）做deim 且不直接计算P
% F的采样来自半降阶
%% 初始化
clear
close all;
%-----------------网格建立------------------------
Re = 1e2;      % Reynolds numbe 
max_it=2000;          % 最大演算步数
dt = 1e-2;            % stepwise of time
nx = 60;      ny = 60;      % number of xy-gridpoints
sample.uv =1000;       dimension.uv =20;    % 采样数与维度
sample.f  =1000;       dimension.f  =20;    % 对流非线性项
sample.rhs=1000;       dimension.rhs=30;    % poisson方程

lx = 1;       ly = 1;       % height of domin
dx = lx/nx;   dy = ly/ny;   % stepwise of u,v
fprintf('-------------------------------------------------------\n')
fprintf('Rayleigh Nnmber:%d; Domain:(0,%d)*(0,%d);Grid:%d*%d\n',Re,lx,ly,nx,ny)
fprintf('stepwise of time:%f; max_iteration:%d\n',dt,max_it)
%---------------------网格建立---------------------
t=0:dt:max_it*dt;
xu=linspace(dx,1-dx,nx-1);  yu=linspace(dy/2,1-dy/2,ny);
xv=linspace(dx/2,1-dx/2,nx);yv=linspace(dy,1-dy,ny-1);
% 初始条件
u=zeros(ny,nx-1);v=zeros(ny-1,nx);p=zeros(ny,nx);    

% 生成数值计算所需矩阵等
Kx=dt/Re/dx/dx; Ky=dt/Re/dy/dy;
%------------对流------------------------
[GammaLu,GammaUu]=C_Gamma(nx-1,ny,1);   [GammaLv,GammaUv]=C_Gamma(ny-1,nx,1);
[GammaLut,GammaUut]=C_Gamma(nx-1,ny,-1);[GammaLvt,GammaUvt]=C_Gamma(ny-1,nx,-1);
[PhiLu,PhiUu]=C_Phi(nx-1,ny,1);   [PhiLv,PhiUv]=C_Phi(ny-1,nx,1);
[PhiLut,PhiUut]=C_Phi(nx-1,ny,-1);[PhiLvt,PhiUvt]=C_Phi(ny-1,nx,-1);
[PsiL1u,PsiL2u,PsiU1u,PsiU2u]=C_Psi(nx,ny);
[PsiL1v,PsiL2v,PsiU1v,PsiU2v]=C_Psi(ny,nx);
%------------扩散------------------------
Hu=speye((nx-1)*ny)-D_J(Kx,Ky,nx-1,ny);
Hv=speye((ny-1)*nx)-D_J(Ky,Kx,ny-1,nx);
%-------------solve Poisson---------------
Lp = kron(speye(ny),K1(nx,dx,1))+kron(K1(ny,dy,1),speye(nx));Lp0=Lp;
Lp(1,1) = 3/2*Lp(1,1); %为了保证正定，而不改变结果
perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
Lq = kron(speye(ny-1),K1(nx-1,dx,2))+kron(K1(ny-1,dy,2),speye(nx-1));
perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
%------------采样信息初始化----------------
n_sample1=floor(max_it/sample.uv);j1=1;
n_sample2=floor(max_it/sample.f);j2=1;
n_sample3=floor(max_it/sample.rhs);j3=1;
podur=zeros((nx-1)*ny,sample.uv);podvr=zeros((ny-1)*nx,sample.uv);
deimf1=zeros((nx-1)*ny,sample.f);deimf2=zeros((ny-1)*nx,sample.f);
podrhs=zeros(nx*ny,sample.rhs);
maxu=zeros((nx-1)*ny,1);maxv=zeros((ny-1)*nx,1);
%% 数值解 层级递推
% 边界条件 _Tilde
u_tilde1=zeros(ny,1);    u_tilde2=zeros(ny,1);
u_tilde=reshape(transpose([u_tilde1,zeros(ny,nx-3),u_tilde2]),[],1);
v_tilde1=zeros(1,nx);    v_tilde2=zeros(1,nx);
v_tilde=reshape([v_tilde1;zeros(ny-3,nx);v_tilde2],[],1);
BC.u2R=reshape(transpose([zeros(ny,nx-2),u_tilde2]),[],1);
BC.u2L=reshape(transpose([u_tilde1,zeros(ny,nx-2)]),[],1);
BC.v2R=reshape([zeros(ny-2,nx);v_tilde2],[],1);
BC.v2L=reshape([v_tilde1;zeros(ny-2,nx)],[],1);
BC.uvyR=PsiU2u*v_tilde;BC.uvyL=PsiL2u*v_tilde;
BC.uvxR=PsiU2v*u_tilde;BC.uvxL=PsiL2v*u_tilde;
tic
for i=1:max_it
% 虚点（插值合成）   _Hat
u_hat1=2*zeros(1,nx-1)-u(1,:);u_hat2=2*ones(1,nx-1)-u(end,:);
u_hat=[u_hat1';zeros((nx-2)*(ny-1),1);u_hat2'];
v_hat1=2*zeros(ny-1,1)-v(:,1);v_hat2=2*zeros(ny-1,1)-v(:,end);
v_hat=[v_hat1;zeros((nx-1)*(ny-2),1);v_hat2];
HAT.uR=[zeros((ny-1)*(nx-1),1);u_hat2'];HAT.uL=[u_hat1';zeros((ny-1)*(nx-1),1)];
HAT.vR=[zeros((nx-1)*(ny-1),1);v_hat2]; HAT.vL=[v_hat1;zeros((nx-1)*(ny-1),1)];

ur=reshape(u',[(nx-1)*ny,1]);vr=reshape(v,[nx*(ny-1),1]);
gamma = min(1.2*dt*max(max(max(abs(u)))/dx,max(max(abs(v)))/dy),1);
% ---------------u--------------------------------------
a=(GammaUu*ur+BC.u2R).^2-(GammaLu*ur+BC.u2L).^2;
a1=abs(GammaUu*ur+BC.u2R).*(GammaUut*ur+BC.u2R)-abs(GammaLu*ur+BC.u2L).*(GammaLut*ur-BC.u2L);
b=PhiUu*ur+HAT.uR;  b1=PhiUut*ur+HAT.uR;  c=PsiU1u*vr+BC.uvyR;
d=PhiLu*ur+HAT.uL;  d1=PhiLut*ur-HAT.uL;  e=PsiL1u*vr+BC.uvyL;
f1=(a-gamma*a1)/dx/4+((b.*c-d.*e)-gamma*(b1.*abs(c)-d1.*abs(e)))/dy/4;
ur=Hu\(ur-dt*f1+Kx*u_tilde+Ky*u_hat);
u=transpose(reshape(ur,[nx-1,ny]));
% ---------------v--------------------------------------
a=(GammaUv*vr+BC.v2R).^2-(GammaLv*vr+BC.v2L).^2;
a1=abs(GammaUv*vr+BC.v2R).*(GammaUvt*vr+BC.v2R)-abs(GammaLv*vr+BC.v2L).*(GammaLvt*vr-BC.v2L);
b=PhiUv*vr+HAT.vR;  b1=PhiUvt*vr+HAT.vR;  c=PsiU1v*ur+BC.uvxR;
d=PhiLv*vr+HAT.vL;  d1=PhiLvt*vr-HAT.vL;  e=PsiL1v*ur+BC.uvxL;
f2=(a-gamma*a1)/dy/4+((b.*c-d.*e)-gamma*(b1.*abs(c)-d1.*abs(e)))/dx/4;
vr=Hv\(vr-dt*f2+Ky*v_tilde+Kx*v_hat);
v=reshape(vr,[ny-1,nx]);
% 非线性项deim采样
if mod(i,n_sample2)==0
    deimf1(:,j2)=f1;deimf2(:,j2)=f2;j2=j2+1;
end
% ---------------pressure update---------------------------
rhs = reshape((diff([u_tilde1,u,u_tilde2],1,2)/dx+diff([v_tilde1;v;v_tilde2])/dy),[],1);
p(perp) = -Rp\(Rpt\rhs(perp));
u=u-diff(p,1,2)/dx;v=v-diff(p)/dy;
% POD采样 (记录最大点)
if mod(i,n_sample1)==0
    podur(:,j1)=reshape(u',[],1);
    [~,h1]=max(podur(:,j1));
    maxu(h1)=h1;
    podvr(:,j1)=reshape(v,[],1);
    [~,h2]=max(podvr(:,j1));
    maxv(h2)=h2;
    j1=j1+1;
end
end
disp('1.完整模型建立：')
toc
maxu(maxu==0)=[];maxv(maxv==0)=[];
u_full=u;v_full=v;  % 保存数据作为对比

% 可视化：完整模型
us=([u_tilde1,u_full]+[u_full,u_tilde2])/2;
vs=([v_tilde1;v_full]+[v_full;v_tilde2])/2;
rhsq = reshape((diff(u_full)/dy-diff(v_full')'/dx),[],1);
q(perq) = Rq\(Rqt\rhsq(perq));
Q=zeros(ny+1,nx+1);Q(2:ny,2:nx)=reshape(q,nx-1,ny-1);
[xx,yy]=meshgrid(linspace(0.5*dx,lx-0.5*dx,nx),linspace(0.5*dy,ly-0.5*dy,ny));
[xx1,yy1]=meshgrid(linspace(0,lx,nx+1),linspace(0,lx,ny+1));
figure("Name",'完整模型结果')
Len = sqrt(us.^2+vs.^2+eps);
quiver(xx,yy,us./Len,vs./Len,.4,'k-')
hold on
contour(xx1,yy1,Q,12,'b','LineWidth',1.5);
axis image
set(gca,'xtick',[]); set(gca,'ytick',[]);
hold off
drawnow
%% 建立ROM
% POD 对 u,v降维
[U,spu,~]=svd(podur);PSIu=U(:,1:dimension.uv);spu=diag(spu);I_u=sum(spu(1:dimension.uv))/sum(spu);
[U,spv,~]=svd(podvr);PSIv=U(:,1:dimension.uv);spv=diag(spv);I_v=sum(spv(1:dimension.uv))/sum(spv);
 fprintf('对U进行SVD降至%d维，能量占比为%f \n',dimension.uv,I_u)
% fprintf('对V进行SVD降至%d维，能量占比为%f \n',dimension.uv,I_v)
% deim 对f1,f2
[U,spf1,~]=svd(deimf1);PHIu=U(:,1:dimension.f);spf1=diag(spf1);
[U,spf2,~]=svd(deimf2);PHIv=U(:,1:dimension.f);spf2=diag(spf2);
[p1,P1]=deim(PHIu);
[p2,P2]=deim(PHIv);
% fprintf('对F进行DEIM维数%d\n',dimension.f)
% fprintf('对P进行SVD降至%d维，能量占比为%f \n',dimension.p,I_p)
disp("2.POD-DEIM正交基计算完成")
%% 半降阶中对Poisson所需数据重采样(为适应降阶Burgers带来的误差)
% 重新初始化
u=zeros(ny,nx-1);   v=zeros(ny-1,nx);   p=zeros(ny,nx);   
% 边界调整
u_tilde=reshape(transpose([u_tilde1,zeros(ny,nx-3),u_tilde2]),[],1);
v_tilde=reshape([v_tilde1;zeros(ny-3,nx);v_tilde2],[],1);
BCr.u2R=P1'*reshape(transpose([zeros(ny,nx-2),u_tilde2]),[],1);%f中的边界乘P
BCr.u2L=P1'*reshape(transpose([u_tilde1,zeros(ny,nx-2)]),[],1);
BCr.v2R=P2'*reshape([zeros(ny-2,nx);v_tilde2],[],1);
BCr.v2L=P2'*reshape([v_tilde1;zeros(ny-2,nx)],[],1);
BCr.uvyR=P1'*PsiU2u*v_tilde;BCr.uvyL=P1'*PsiL2u*v_tilde;
BCr.uvxR=P2'*PsiU2v*u_tilde;BCr.uvxL=P2'*PsiL2v*u_tilde;
BCr.lapu=Kx*PSIu'*u_tilde;BCr.lapv=Ky*PSIv'*v_tilde;
% 矩阵重定义
U2.R =P1'*GammaUu* PSIu; U2.L= P1'*GammaLu* PSIu;
U2.Rt=P1'*GammaUut*PSIu; U2.Lt=P1'*GammaLut*PSIu;
V2.R =P2'*GammaUv* PSIv; V2.L= P2'*GammaLv* PSIv;
V2.Rt=P2'*GammaUvt*PSIv; V2.Lt=P2'*GammaLvt*PSIv;
UVy.uR =P1'*PhiUu *PSIu; UVy.uL =P1'*PhiLu *PSIu;
UVy.uRt=P1'*PhiUut*PSIu; UVy.uLt=P1'*PhiLut*PSIu;
UVy.vR= P1'*PsiU1u*PSIv; UVy.vL= P1'*PsiL1u*PSIv;
UVx.vR =P2'*PhiUv *PSIv; UVx.vL =P2'*PhiLv *PSIv;
UVx.vRt=P2'*PhiUvt*PSIv; UVx.vLt=P2'*PhiLvt*PSIv;
UVx.uR= P2'*PsiU1v*PSIu; UVx.uL= P2'*PsiL1v*PSIu;
DEIM.u1=dt*PSIu'*PHIu;   DEIM.u2=P1'*PHIu;
DEIM.v1=dt*PSIv'*PHIv;   DEIM.v2=P2'*PHIv;
H1=PSIu'*(Hu\PSIu); H2=PSIv'*(Hv\PSIv);
tic
for i=1:max_it
% 虚点（插值合成）   _Hat
u_hat1=2*zeros(1,nx-1)-u(1,:);u_hat2=2*ones(1,nx-1)-u(end,:);
u_hat=[u_hat1';zeros((nx-2)*(ny-1),1);u_hat2'];
v_hat1=2*zeros(ny-1,1)-v(:,1);v_hat2=2*zeros(ny-1,1)-v(:,end);
v_hat=[v_hat1;zeros((nx-1)*(ny-2),1);v_hat2];
HATr.uR=P1'*[zeros((ny-1)*(nx-1),1);u_hat2'];HATr.uL=P1'*[u_hat1';zeros((ny-1)*(nx-1),1)];
HATr.vR=P2'*[zeros((nx-1)*(ny-1),1);v_hat2]; HATr.vL=P2'*[v_hat1;zeros((nx-1)*(ny-1),1)];

ur=reshape(u',[(nx-1)*ny,1]);vr=reshape(v,[nx*(ny-1),1]);
alpha=PSIu'*ur;     beta=PSIv'*vr;
gamma = min(1.2*dt*max(max(max(abs(u)))/dx,max(max(abs(v)))/dy),1);
% ---------------u--------------------------------------
a=(U2.R*alpha+BCr.u2R).^2-(U2.L*alpha+BCr.u2L).^2;
a1=abs(U2.R*alpha+BCr.u2R).*(U2.Rt*alpha+BCr.u2R)-abs(U2.L*alpha+BCr.u2L).*(U2.Lt*alpha-BCr.u2L);
b=UVy.uR*alpha+HATr.uR;  b1=UVy.uRt*alpha+HATr.uR;  c=UVy.vR*beta+BCr.uvyR;
d=UVy.uL*alpha+HATr.uL;  d1=UVy.uLt*alpha-HATr.uL;  e=UVy.vL*beta+BCr.uvyL;
f1=(a-gamma*a1)/dx/4+((b.*c-d.*e)-gamma*(b1.*abs(c)-d1.*abs(e)))/dy/4;
alpha=H1*(alpha+BCr.lapu+Ky*PSIu'*u_hat)-DEIM.u1*(DEIM.u2\f1);
u=transpose(reshape(PSIu*alpha,[nx-1,ny]));
% ---------------v--------------------------------------
a=(V2.R*beta+BCr.v2R).^2-(V2.L*beta+BCr.v2L).^2;
a1=abs(V2.R*beta+BCr.v2R).*(V2.Rt*beta+BCr.v2R)-abs(V2.L*beta+BCr.v2L).*(V2.Lt*beta-BCr.v2L);
b=UVx.vR*beta+HATr.vR;   b1=UVx.vRt*beta+HATr.vR;   c=UVx.uR*alpha+BCr.uvxR;
d=UVx.vL*beta+HATr.vL;   d1=UVx.vLt*beta-HATr.vL;   e=UVx.uL*alpha+BCr.uvxL;
f2=(a-gamma*a1)/dy/4+((b.*c-d.*e)-gamma*(b1.*abs(c)-d1.*abs(e)))/dx/4;
beta=H2*(beta+BCr.lapv+Kx*PSIv'*v_hat)-DEIM.v1*(DEIM.v2\f2);
v=reshape(PSIv*beta,[ny-1,nx]);
% ---------------pressure update---------------------------
rhs = reshape((diff([u_tilde1,u,u_tilde2],1,2)/dx+diff([v_tilde1;v;v_tilde2])/dy),[],1);
p(perp) = -Rp\(Rpt\rhs(perp));
u=u-diff(p,1,2)/dx;v=v-diff(p)/dy;
if mod(i,n_sample3)==0
    podrhs(:,j3)=reshape(rhs,[],1);
    j3=j3+1;
end
end
[U,~,~]=svd(podrhs); PSIrhs=U(:,1:dimension.rhs);
% 对nablaUV做DEIM
% [sf,S]=deim(PSIrhs);
S=q_deim(PSIrhs);
disp("3.对poisson方程非齐次项F采样与降维完成")
%% 降维模型
% 重新初始化
u=zeros(ny,nx-1);v=zeros(ny-1,nx);p=zeros(ny,nx);
alpha=zeros(dimension.uv,1);  beta=zeros(dimension.uv,1);  

nablaU=S'*acute(nx,ny)*Phi(nx-1,ny,-1)*PSIu; nablaV=S'*Phi(ny-1,nx,-1)*PSIv;
BCr.nu=S'*acute(nx,ny)*reshape(transpose([-1*u_tilde1,zeros(ny,nx-2),u_tilde2]),[],1);
BCr.nv=S'*reshape([-1*v_tilde1;zeros(ny-2,nx);v_tilde2],[],1);
% lapP=PSIp'*inv(Lp)*PSIrhs*inv(S'*PSIrhs);
nabla1=Psi1(nx,ny,-1);  nabla2=Phi1(ny,nx,-1);
Pois1=PSIu'*acute(ny,nx-1)*nabla1*inv(Lp)*  PSIrhs*inv(S'*PSIrhs);
Pois2=PSIv'*nabla2               *inv(Lp)*  PSIrhs*inv(S'*PSIrhs);

% 虚点处理
hatN=2*[zeros((nx-1)*(ny-1),1);ones(nx-1,1)];
hatS=2*zeros(ny*(nx-1),1);
hatE=2*zeros(nx*(ny-1),1);
hatW=2*zeros(nx*(ny-1),1);

HaU=Ky*PSIu'*(hatN+hatS);
HaV=Kx*PSIv'*(hatE+hatW);
HaNf1=P1'*hatN;    HaSf1=P1'*hatS;
HaWf2=P2'*hatW;    HaEf2=P2'*hatE;

% 筛选矩阵
north=spdiags([zeros((nx-1)*(ny-1),1);ones(nx-1,1)],0,(nx-1)*ny,(nx-1)*ny);
south=spdiags([ones(nx-1,1);zeros((nx-1)*(ny-1),1)],0,(nx-1)*ny,(nx-1)*ny);
west =spdiags([zeros((nx-1)*(ny-1),1);ones(ny-1,1)],0,(ny-1)*nx,(ny-1)*nx);
east =spdiags([ones(ny-1,1);zeros((nx-1)*(ny-1),1)],0,(ny-1)*nx,(ny-1)*nx);

hat_NS_u =Ky*PSIu'*(north+south)*PSIu;
hat_WE_v =Kx*PSIv'*( west+ east)*PSIv;
hat_N_f1=P1'*north*PSIu;  hat_S_f1=P1'*south*PSIu;
hat_W_f2=P2'* west*PSIv;  hat_E_f2=P2'* east*PSIv;

pmu=PSIu(maxu,:);pmv=PSIv(maxv,:);% 低维gamma用正交基
tic
for i=1:max_it
gamma = min(1.2*dt*max(max(abs(pmu*alpha))/dx,max(abs(pmv*beta))/dy),1);
% gamma =0;
HATr.u=HaU-hat_NS_u*alpha;      HATr.v=HaV-hat_WE_v* beta;
HATr.uR=HaNf1-hat_N_f1*alpha;   HATr.uL=HaSf1-hat_S_f1*alpha;
HATr.vR=HaWf2-hat_W_f2* beta;   HATr.vL=HaEf2-hat_E_f2* beta;
% ---------------u--------------------------------------
a=(U2.R*alpha+BCr.u2R).^2-(U2.L*alpha+BCr.u2L).^2;
a1=abs(U2.R*alpha+BCr.u2R).*(U2.Rt*alpha+BCr.u2R)-abs(U2.L*alpha+BCr.u2L).*(U2.Lt*alpha-BCr.u2L);
b=UVy.uR*alpha+HATr.uR;  b1=UVy.uRt*alpha+HATr.uR;  c=UVy.vR*beta+BCr.uvyR;
d=UVy.uL*alpha+HATr.uL;  d1=UVy.uLt*alpha-HATr.uL;  e=UVy.vL*beta+BCr.uvyL;
f1=(a-gamma*a1)/dx/4+((b.*c-d.*e)-gamma*(b1.*abs(c)-d1.*abs(e)))/dy/4;
alpha=H1*(alpha+BCr.lapu+HATr.u)-DEIM.u1*(DEIM.u2\f1);
% ---------------v--------------------------------------
a=(V2.R*beta+BCr.v2R).^2-(V2.L*beta+BCr.v2L).^2;
a1=abs(V2.R*beta+BCr.v2R).*(V2.Rt*beta+BCr.v2R)-abs(V2.L*beta+BCr.v2L).*(V2.Lt*beta-BCr.v2L);
b=UVx.vR*beta+HATr.vR;   b1=UVx.vRt*beta+HATr.vR;   c=UVx.uR*alpha+BCr.uvxR;
d=UVx.vL*beta+HATr.vL;   d1=UVx.vLt*beta-HATr.vL;   e=UVx.uL*alpha+BCr.uvxL;
f2=(a-gamma*a1)/dy/4+((b.*c-d.*e)-gamma*(b1.*abs(c)-d1.*abs(e)))/dx/4;
beta=H2*(beta+BCr.lapv+HATr.v)-DEIM.v1*(DEIM.v2\f2);
% ---------------pressure update---------------------------
rhs=(nablaU*alpha+BCr.nu)/dx+(nablaV*beta+BCr.nv)/dy; % 次出rhs是低维的
alpha=alpha+(Pois1*rhs)/dx;
beta =beta +(Pois2*rhs)/dy;
end
toc
disp("4.降维模型运算完成")

u_deim=transpose(reshape(PSIu*alpha,[nx-1,ny]));
v_deim=reshape(PSIv*beta,[ny-1,nx]);

% 可视化：降维模型 误差显示
% 误差分析
Eu=norm(u_full-u_deim)/norm(u_full);
RMSE=sqrt(sum(sum((u_full-u_deim).^2))/(nx-1)/ny);
fprintf('U的相对误差%f\n',Eu)
fprintf('U的均方误差%f\n',RMSE)
% visualization
us=([u_tilde1,u_deim]+[u_deim,u_tilde2])/2;
vs=([v_tilde1;v_deim]+[v_deim;v_tilde2])/2;
rhsq = reshape((diff(u_deim)/dy-diff(v_deim')'/dx),[],1);
q(perq) = Rq\(Rqt\rhsq(perq));
Q=zeros(ny+1,nx+1);Q(2:ny,2:nx)=reshape(q,nx-1,ny-1);

figure("Name",'降维模型结果')
Len = sqrt(us.^2+vs.^2+eps);
quiver(xx,yy,us./Len,vs./Len,.4,'k-')
hold on
contour(xx1,yy1,Q,12,'b','LineWidth',1.5);
axis image
set(gca,'xtick',[]); set(gca,'ytick',[]);
hold off

%  对比与误差
% figure
% subplot(121),surf(xu,yu,u_full),title('full')
% subplot(122),surf(xu,yu,u_deim),title('deim')
% figure
% [C,h]=contourf(xu,yu,abs(u_full-u_deim),50);
% set(h,'Color','none'),title('error')
% axis square
% colorbar

%% 函数 生成矩阵
% n对应一行几个元素  m对应整体几行
%---------------对流（平方）--------------------
function [GammaL,GammaU]=C_Gamma(n,m,a)
U=spdiags([a*ones(n,1),ones(n,1)],0:1,n,n);
L=spdiags([a*ones(n,1),ones(n,1)],-1:0,n,n);
E=speye(m);
GammaL=kron(E,L);GammaU=kron(E,U);
end
function A1=Phi(n,m,a)
A=spdiags([a*ones(n,1),ones(n,1)],-1:0,n+1,n);
A1=kron(spdiags(ones(m,1),0,m,m),A);
end
function A1=Phi1(n,m,a)
A=spdiags([a*ones(n-1,1),ones(n-1,1)],0:1,n-1,n);
A1=kron(spdiags(ones(m,1),0,m,m),A);
end
function B=Psi1(n,m,a)
B=spdiags([a*ones(n*(m-1),1),ones(n*(m-1),1)],[0,n],n*(m-1),m*n);
end
%---------------对流（混合）--------------------
function [PhiL,PhiU]=C_Phi(n,m,a)
% 对流：混合项同名
U=spdiags([a*ones(m,1),ones(m,1)],0:1,m,m);
L=spdiags([a*ones(m,1),ones(m,1)],-1:0,m,m);
E=speye(n);
PhiL=kron(L,E);PhiU=kron(U,E);    
end

function [PsiL1,PsiL2,PsiU1,PsiU2]=C_Psi(n,m)
% 对流：混合项不同名
zero=sparse(n-1,n*(m-1));
P=[];E=speye(n-1,n)+sparse(1:n-1,2:n,1,n-1,n);
for i=1:m-1
    row=sparse(1,i,1,1,m-1);
    P=[P;kron(E,row)];
end
PsiL1=[zero;P];PsiU1=[P;zero];
PsiL2=[kron(E,sparse(1,1,1,1,m-1));zeros(size(P))];
PsiU2=[zeros(size(P));kron(E,sparse(1,m-1,1,1,m-1))];
end
%---------------扩散----------------------------
function J=D_J(a,b,n,m)
x=a*ones(n,1);y=b*ones(m,1);
A = spdiags([x -2*x x],-1:1,n,n);
B = spdiags([y -2*y y],-1:1,m,m);
J=kron(speye(m),A)+kron(B,speye(n));
end
function ua=acute(n,m)
t=reshape(transpose(reshape(1:n*m,n,m)),[],1);
ua=sparse(1:n*m,t,ones(n*m,1));
end
%---------------Chol-------------------------
function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;
end
%-----------------deim------------
function [ S, M ] = q_deim( U )
% M为正交基变体
% 内置函数qr出来的P为列主元
% 这里将S的输出改为降维矩阵P
[n,m] = size(U) ;
if nargout == 1
    [~,~,P] = qr(U','vector') ;p=P(1:m);
else
    [Q,R,P] = qr(U','vect') ;p=P(1:m);
    M = [eye(m) ; (R(:,1:m)\R(:,m+1:n))'] ;
    Pinverse(P) = 1 : n ; M = M(Pinverse,:) ;
end
S=zeros(n,m);
for i=1:m
    S(p(i),i)=1;
end
end

