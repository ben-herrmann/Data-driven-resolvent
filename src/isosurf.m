
function fig = isosurf(x,y,z,U,c,color1,color2)
fig = figure;
[X,Z,Y] = meshgrid(x-pi, z-pi, y);
ax = gca;
ax.BoxStyle = 'full';
box on
p = patch(isosurface(X,Z,Y,U,c));
isonormals(X,Z,Y,U,p)
set(p,'FaceColor',color1,'EdgeColor','none');
daspect([1 1 1])
view([20,25]);
camlight
lighting gouraud
view([-70,25]);
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);

hold on
p = patch(isosurface(X,Z,Y,U,-c));
isonormals(X,Z,Y,U,p)
set(p,'FaceColor',color2,'EdgeColor','none');
daspect([1 1 1])
axis([X(1) X(end) Z(1) Z(end) Y(1) Y(end)]);
lighting gouraud
xticks([])
yticks([])
zticks([])
hold off
end