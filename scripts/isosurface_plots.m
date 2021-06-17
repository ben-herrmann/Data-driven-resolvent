clear all, close all, clc
% Project/scripts should be set as the working directory

addpath('../src')
red = [1 0 0];
green = [0 0.8 0];
yellow = [1 1 0];
blue = [0 0 0.7];
white = [0.9 0.9 0.9];
black = [0.1 0.1 0.1];

%% DATASET 1

%% Snapshots

dir_load = '../data/sims/channel_flow/random_snapshots.mat';
load(dir_load)

c = 0.5;
l = 2;
k = [2];%,102,202];
scale = max(abs(X{1}(:)));
for i=1:length(k)
u = get_component(X{k(i)},l,Nx,Ny,Nz);
u = u/max(abs(u(:)));
fig(k) = isosurf(x,y(end:-1:1),z,u,c,red,green);
print(fig(k),'../plots/random_x0','-depsc');
end

%% Data-driven resolvent modes

% w=0, mode 1

dir_load = '../data/sims/channel_flow/resolvent_rand_w=0.0_mode1.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y,z,V,0.5,red,green);
print(f1,'../plots/phi_random_w=0_mode1','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y,z,U,0.5,yellow,blue);
print(f2,'../plots/psi_random_w=0_mode1','-depsc');

% w=0, mode 2

dir_load = '../data/sims/channel_flow/resolvent_rand_w=0.0_mode2.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y,z,V,0.5,red,green);
print(f1,'../plots/phi_random_w=0_mode2','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y,z,U,0.5,yellow,blue);
print(f2,'../plots/psi_random_w=0_mode2','-depsc');

% w=0, mode 3

dir_load = '../data/sims/channel_flow/resolvent_rand_w=0.0_mode3.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y,z,V,0.5,red,green);
print(f1,'../plots/phi_random_w=0_mode3','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y,z,U,0.5,yellow,blue);
print(f2,'../plots/psi_random_w=0_mode3','-depsc');

% w=0.4, mode 1

dir_load = '../data/sims/channel_flow/resolvent_rand_w=0.4_mode1.mat';
load(dir_load)

V = get_component(phi,1,Nx,Ny,Nz);
f1 = isosurf(x,y,z,V,0.5,yellow,blue);
print(f1,'../plots/phi_random_w=04_mode1','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y,z,U,0.5,yellow,blue);
print(f2,'../plots/psi_random_w=04_mode1','-depsc');

%% DATASET 2

dir_load = '../data/sims/channel_flow/optimal_snapshots.mat';
load(dir_load)

c = 0.5;
l = 2;
k = [1];%,102,202];
scale = max(abs(X{1}(:)));
for i=1:length(k)
u = get_component(X{k(i)},l,Nx,Ny,Nz);
u = u/max(abs(u(:)));
fig(k) = isosurf(x,y(end:-1:1),z,u,c,red,green);
print(fig(k),'../plots/optimal_x0','-depsc');
end

%% Data-driven resolvent modes

% w=0, mode 1

dir_load = '../data/sims/channel_flow/resolvent_optimal_w=0.0_mode1.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y(end:-1:1),z,V,0.5,red,green);
print(f1,'../plots/phi_optimal_w=0_mode1','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y(end:-1:1),z,U,0.5,yellow,blue);
print(f2,'../plots/psi_optimal_w=0_mode1','-depsc');

%% DATASET 2

dir_load = '../data/sims/channel_flow/localized_snapshots.mat';
load(dir_load)

c = 0.05;
l = 2;
k = [1];%,102,202];
scale = max(abs(X{1}(:)));
for i=1:length(k)
u = get_component(X{k(i)},l,Nx,Ny,Nz);
u = u/max(abs(u(:)));
fig(k) = isosurf(x,y(end:-1:1),z,u,c,red,green);
print(fig(k),'../plots/localized_x0','-depsc');
end

%% Data-driven resolvent modes

% w=0, mode 1

dir_load = '../data/sims/channel_flow/resolvent_localized_w=0.0_mode1.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y(end:-1:1),z,V,0.3,red,green);
print(f1,'../plots/phi_localized_w=0_mode1','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y(end:-1:1),z,U,0.3,yellow,blue);
print(f2,'../plots/psi_localized_w=0_mode1','-depsc');

%% Operator-based projected resolvent modes

% w=0, mode 1

dir_load = '../data/sims/channel_flow/resolvent_localized_projected.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y(end:-1:1),z,V,0.3,red,green);
print(f1,'../plots/phi_localized_projected','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y(end:-1:1),z,U,0.3,yellow,blue);
print(f2,'../plots/psi_localized_projected','-depsc');

%% OPERATOR-BASED

%% True resolvent modes

% w=0, mode 1

dir_load = '../data/sims/channel_flow/resolvent_truth_w=0.0_mode1.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y,z,V,0.5,red,green);
print(f1,'../plots/phi_truth_w=0_mode1','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y,z,U,0.5,yellow,blue);
print(f2,'../plots/psi_truth_w=0_mode1','-depsc');

% w=0, mode 2

dir_load = '../data/sims/channel_flow/resolvent_truth_w=0.0_mode2.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y,z,V,0.5,red,green);
print(f1,'../plots/phi_truth_w=0_mode2','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y,z,U,0.5,yellow,blue);
print(f2,'../plots/psi_truth_w=0_mode2','-depsc');

% w=0, mode 3

dir_load = '../data/sims/channel_flow/resolvent_truth_w=0.0_mode3.mat';
load(dir_load)

V = get_component(phi,2,Nx,Ny,Nz);
f1 = isosurf(x,y,z,V,0.5,red,green);
print(f1,'../plots/phi_truth_w=0_mode3','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y,z,U,0.5,yellow,blue);
print(f2,'../plots/psi_truth_w=0_mode3','-depsc');

% w=0.4, mode 1

dir_load = '../data/sims/channel_flow/resolvent_truth_w=0.4_mode1.mat';
load(dir_load)

V = get_component(phi,1,Nx,Ny,Nz);
f1 = isosurf(x,y,z,V,0.5,yellow,blue);
print(f1,'../plots/phi_truth_w=04_mode1','-depsc');

U = get_component(psi,1,Nx,Ny,Nz);
f2 = isosurf(x,y,z,U,0.5,yellow,blue);
print(f2,'../plots/psi_truth_w=04_mode1','-depsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
dir_load = '../data/sims/channel_flow/resolvent_modes.mat';
load(dir_load)

%% Load operator-based resolvent modes
%% w=0, mode 1
dir_load = '../data/sims/channel_flow/resolvent_truth_w=0.0_mode1.mat';
load(dir_load)
%% Get components
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = -V/max(abs(u(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,Ny+1-j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(u(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.5));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.5));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
% xlab = xlabel('x', 'interpreter', 'tex');
% ylab = ylabel('z', 'interpreter', 'tex');
% zlab = zlabel('y', 'interpreter', 'tex');
% set(xlab, 'FontName', 'cmmi10', 'FontSize', 25);
% set(ylab, 'FontName', 'cmmi10', 'FontSize', 25);
% set(zlab, 'FontName', 'cmmi10', 'FontSize', 25);
% h = get(gcf,'CurrentAxes');
% set(h,'FontName','cmr10','FontSize',15,'xscale','lin','yscale','lin');
hold off
% print(f1,'../plots/phi_truth_w=0_mode1','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.5));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.5));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
% print(f2,'../plots/psi_truth_w=0_mode1','-depsc');

%% w=0, mode 2
dir_load = '../data/sims/channel_flow/resolvent_truth_w=0.0_mode2.mat';
load(dir_load)
%% Get components
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = -V/max(abs(u(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,Ny+1-j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(u(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.5));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.5));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
% xlab = xlabel('x', 'interpreter', 'tex');
% ylab = ylabel('z', 'interpreter', 'tex');
% zlab = zlabel('y', 'interpreter', 'tex');
% set(xlab, 'FontName', 'cmmi10', 'FontSize', 25);
% set(ylab, 'FontName', 'cmmi10', 'FontSize', 25);
% set(zlab, 'FontName', 'cmmi10', 'FontSize', 25);
% h = get(gcf,'CurrentAxes');
% set(h,'FontName','cmr10','FontSize',15,'xscale','lin','yscale','lin');
hold off
print(f1,'../plots/phi_truth_w=0_mode2','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.5));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.5));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_truth_w=0_mode2','-depsc');

%% w=0, mode 3
dir_load = '../data/sims/channel_flow/resolvent_truth_w=0.0_mode3.mat';
load(dir_load)
%% Get components
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = -V/max(abs(u(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,Ny+1-j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(u(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.5));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.5));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/phi_truth_w=0_mode3','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.5));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.5));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_truth_w=0_mode3','-depsc');

%% Load operator-based resolvent modes
dir_load = '../data/sims/channel_flow/resolvent_truth.mat';
load(dir_load)
%% Get components
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = V/max(abs(u(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(u(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.5));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.5));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
% xlab = xlabel('x', 'interpreter', 'tex');
% ylab = ylabel('z', 'interpreter', 'tex');
% zlab = zlabel('y', 'interpreter', 'tex');
% set(xlab, 'FontName', 'cmmi10', 'FontSize', 25);
% set(ylab, 'FontName', 'cmmi10', 'FontSize', 25);
% set(zlab, 'FontName', 'cmmi10', 'FontSize', 25);
% h = get(gcf,'CurrentAxes');
% set(h,'FontName','cmr10','FontSize',15,'xscale','lin','yscale','lin');
hold off
print(f1,'../plots/phi_truth','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.5));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.5));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_truth','-depsc');

%% Load small-wavenumber random disturbances mode 1
dir_load = '../data/sims/channel_flow/resolvent_rand_w=0.0_mode1.mat';
load(dir_load)
XX = X;

%% Resolvent modes
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = -V/max(abs(u(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,Ny+1-j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(u(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.5));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.5));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/phi_randomw=0_mode1','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.5));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.5));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_random_w=0_mode1','-depsc');

%% Load small-wavenumber random disturbances mode 2
dir_load = '../data/sims/channel_flow/resolvent_rand_w=0.0_mode2.mat';
load(dir_load)
XX = X;

%% Resolvent modes
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = -V/max(abs(u(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,Ny+1-j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(u(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.5));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.5));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/phi_randomw=0_mode2','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.5));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.5));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_random_w=0_mode2','-depsc');

%% Load small-wavenumber random disturbances mode 2
dir_load = '../data/sims/channel_flow/resolvent_rand_w=0.0_mode3.mat';
load(dir_load)
XX = X;

%% Resolvent modes
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = -V/max(abs(u(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,Ny+1-j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(u(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.5));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.5));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/phi_randomw=0_mode3','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.5));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.5));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_random_w=0_mode3','-depsc');


%% Load small-wavenumber random disturbances
dir_load = '../data/sims/channel_flow/resolvent_random.mat';
load(dir_load)
XX = X;

%% Resolvent modes
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = -V/max(abs(u(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,Ny+1-j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(u(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.5));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.5));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/phi_random','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.5));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.5));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_random','-depsc');

%% Initial conditions
d = 999;
u = reshape(XX(:,1),Nz,Ny,Nx,3);
V1 = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V1(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V1 = V1/max(abs(V1(:)));

u = reshape(XX(:,1+d),Nz,Ny,Nx,3);
V2 = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V2(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
V2 = V2/max(abs(V2(:)));

u = reshape(XX(:,1+2*d),Nz,Ny,Nx,3);
V3 = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V3(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
V3 = V3/max(abs(V3(:)));

%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V1,0.5));
isonormals(X,Z,Y,V1,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V1,-0.5));
isonormals(X,Z,Y,V1,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/x0_1_random','-depsc');

f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V2,0.5));
isonormals(X,Z,Y,V2,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V2,-0.5));
isonormals(X,Z,Y,V2,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/x0_2_random','-depsc');

f3 = figure(3);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V3,0.5));
isonormals(X,Z,Y,V3,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V3,-0.5));
isonormals(X,Z,Y,V3,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f3,'../plots/x0_3_random','-depsc');

%% Transients
d = 2*999;
u = reshape(XX(:,1+d),Nz,Ny,Nx,3);
V1 = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V1(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V1 = V1/max(abs(V1(:)));

u = reshape(XX(:,14+d),Nz,Ny,Nx,3);
V2 = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V2(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
V2 = V2/max(abs(V2(:)));

u = reshape(XX(:,27+d),Nz,Ny,Nx,3);
V3 = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V3(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
V3 = V3/max(abs(V3(:)));

u = reshape(XX(:,40+d),Nz,Ny,Nx,3);
V4 = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V4(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
V4 = V4/max(abs(V4(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V1,0.5));
isonormals(X,Z,Y,V1,p)
% red
set(p,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V1,-0.5));
isonormals(X,Z,Y,V1,p)
% green
set(p,'FaceColor',[0.1 0.1 0.1],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/x1_3_random','-depsc');

f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V2,0.5));
isonormals(X,Z,Y,V2,p)
% red
set(p,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V2,-0.5));
isonormals(X,Z,Y,V2,p)
% green
set(p,'FaceColor',[0.1 0.1 0.1],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/x2_3_random','-depsc');

f3 = figure(3);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V3,0.5));
isonormals(X,Z,Y,V3,p)
% red
set(p,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V3,-0.5));
isonormals(X,Z,Y,V3,p)
% green
set(p,'FaceColor',[0.1 0.1 0.1],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f3,'../plots/x3_3_random','-depsc');

f4 = figure(4);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V4,0.5));
isonormals(X,Z,Y,V4,p)
% red
set(p,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V4,-0.5));
isonormals(X,Z,Y,V4,p)
% green
set(p,'FaceColor',[0.1 0.1 0.1],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f4,'../plots/x4_3_random','-depsc');

%% Load localized disturbances
dir_load = '../data/sims/channel_flow/resolvent_localized.mat';
load(dir_load)
XX = X;
%% Initial condition
u = reshape(XX(:,1),Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = V/max(u(:));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.05));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.05));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/x1_localized','-depsc');
%% Resolvent modes
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,Ny+1-j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = -V/max(abs(V(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,Ny+1-j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, sort(y));
U = U/max(abs(U(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.3));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.3));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/phi_localized','-depsc');

close all
f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.3));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
% camlight
% lighting gouraud
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.3));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_localized','-depsc');
%% Load resolvent projected disturbances
dir_load = '../data/sims/channel_flow/resolvent_projected.mat';
load(dir_load)
%% Resolvent modes
u = reshape(phi,Nz,Ny,Nx,3);
V = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            V(i,k,j) = u(i,j,k,2);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
V = V/max(abs(V(:)));

u = reshape(psi,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,j,k,1);
        end
    end
end
[X Z Y] = meshgrid(x-pi, z-pi, y);
U = U/max(abs(U(:)));
%% Figures
close all
f1 = figure(1);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,V,0.3));
isonormals(X,Z,Y,V,p)
% red
set(p,'FaceColor',[1 0 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,V,-0.3));
isonormals(X,Z,Y,V,p)
% green
set(p,'FaceColor',[0 0.8 0],'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
% camlight
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
print(f1,'../plots/phi_projected','-depsc');

f2 = figure(2);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,0.3));
isonormals(X,Z,Y,U,p)
% yellow
set(p,'FaceColor',[1 1 0],'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
% camlight
% lighting gouraud
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-0.3));
isonormals(X,Z,Y,U,p)
% blue
set(p,'FaceColor',[0 0 0.7],'EdgeColor','none');
% camlight
lighting gouraud

xticks([])
yticks([])
zticks([])
hold off
print(f2,'../plots/psi_projected','-depsc');