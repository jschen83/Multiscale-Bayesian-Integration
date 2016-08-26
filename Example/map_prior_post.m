clear all;
clf;
%
fid=fopen('./Output/mytest_run1.est','r');
SS=fgetl(fid);
tmpint=fscanf(fid,'%d');
nx=tmpint(1);
ny=tmpint(2);
SS=fgetl(fid);
A=fscanf(fid,'%f',[6 inf]);
fclose(fid);
primu=reshape(A(3,:),nx,ny);
postmu=reshape(A(4,:),nx,ny);
pristd=reshape(A(5,:),nx,ny);
poststd=reshape(A(6,:),nx,ny);
%
xx=A(1,1:nx);
yy=A(2,1:nx:(nx*ny));
%
subplot(2,2,1);
H1=pcolor(xx,yy,primu);
set(H1,'LineStyle','none');
colorbar('vert');
title('(a) Prior mean');
subplot(2,2,2);
H2=pcolor(xx,yy,pristd);
set(H2,'LineStyle','none');
colorbar('vert');
title('(b) Prior std');
subplot(2,2,3);
H3=pcolor(xx,yy,postmu);
set(H3,'LineStyle','none');
colorbar('vert');
title('(c) Posterior mean');
subplot(2,2,4);
H4=pcolor(xx,yy,poststd);
set(H4,'LineStyle','none');
colorbar('vert');
title('(d) Posterior std');
