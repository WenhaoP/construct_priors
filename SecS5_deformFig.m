% Part 1
x1 = linspace(0, 3, 21);
y1 = linspace(0, 3, 21);

[x,y] = meshgrid(x1, y1);

X = 3/pi*cos(y*pi/3);
Z = 3/pi*sin(y*pi/3);
Y = x;

surf(X,Y,Z, 0*Z);
hold on;
% 

x1 = linspace(1, 4, 21);
y1 = linspace(0, pi, 21);

[x,y] = meshgrid(x1, y1);

X = 3/pi*x.*cos(y);
Z = 3/pi*x.*sin(y);
Y = 0.5*(x.*sqrt(x.^2-1)-log(sqrt(x.^2-1)+x))+3;

surf(X,Y,Z, 0*Z)

% Part 3
x1 = linspace(0, 3, 20/4+1);
y1 = linspace(0, 3, 21);

[x,y] = meshgrid(x1, y1);

X = 3/pi*4*cos(y*pi/3);
Z = 3/pi*4*sin(y*pi/3);
Y = x+0.5*(4*sqrt(4^2-1)-log(sqrt(4^2-1)+4))+3;

surf(X,Y,Z, 0*Z);
hold on;

axis equal
xlim([-4 4])
view(63, 20)
set(gca, 'FontSize', 18)
print('-depsc2', 'Figures/VaryingR.eps')
