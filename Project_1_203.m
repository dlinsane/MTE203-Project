clear all;
clc;
%1A
[x, y] = meshgrid(1:0.045:5.5, -4:0.07:3);
z = 0.00125 * exp(-((0.5*y.^2)+(x-3).^2)).*(sin(2*x)+2*sin(0.75*(0.5*y-2).^2)).*(16*x+64*x.^2 + y.^2);

figure(1)
surf(x,y,z)
xlabel('East-West Distance (km)')
ylabel('North-South Distance (km)')
zlabel('Elevation (km)')
title('Terrain Plot')

figure(2)
[c, h]=contour(x, y, z, 25);
legend = colorbar;
clabel(c,h, 'FontSize', 6);
xlabel('East-West Distance (km)')
ylabel('North-South Distance (km)')
ylabel(legend, 'Elevation (km)')
title('Contour Plot of Terrain')


%1B
syms x y z T zz l zzz
z(x,y)= 0.00125 * exp(-((0.5*y^2)+(x-3)^2))*(sin(2*x)+2*sin(0.75*(0.5*y-2)^2))*(16*x+64*x^2 + y^2);
dzdx(x,y) = diff(z, x);
dzdy(x,y) = diff(z, y);
slope = sqrt(dzdx.^2 + dzdy.^2);
gradientSlope = gradient(slope);
gradientSlopeMatlab = matlabFunction(gradientSlope, 'Vars', {[x y]});
guess = [3,0];
pointMaxSlope = fsolve(gradientSlopeMatlab, guess);

%1C
gradientZ = gradient(z);
gradientZMatlab = matlabFunction(gradientZ, 'Vars', {[x y]});
guess1 = [3.5,1];
guess2 = [3,-1];
guess3 = [3.5,-2];
cp1=fsolve(gradientZMatlab, guess1);
cp2=fsolve(gradientZMatlab, guess2);
cp3=fsolve(gradientZMatlab, guess3);
A = diff(dzdx, x);
B = diff(dzdx, y);
C = diff(dzdy, y);
Astatus = deriv2Test(A, B, C, cp1, x, y);
Bstatus = deriv2Test(A, B, C, cp2, x, y);
Cstatus = deriv2Test(A, B, C, cp3, x, y);
maxHeight = double(z(cp1(1), cp1(2)));
minHeight = double(z(cp2(1), cp2(2)));
localMax = double(z(cp3(1), cp3(2)));

%2A
T(x,y,zz) = (-0.1*(zz^2))+17*exp(-0.1*((0.1*x-2)-(0.05*y-1)^2-(zz-1)^2))-10;
tempAtPeak = double(T(cp1(1), cp1(2), maxHeight));
tempAtValley = double(T(cp2(1), cp2(2), minHeight));

%2B
heightAtPoint = z(4,-0.3);
tempAtPoint = T(4,-0.3,heightAtPoint);
[xx, yy] = meshgrid(1:0.09:5.5,-4:0.14:3);
isotherm = (-0.1*(heightAtPoint^2))+17*exp(-0.1*((0.1*xx-2)-(0.05*yy-1).^2-(heightAtPoint-1)^2))-10;

figure(3)
[c, h] = contour(xx,yy,isotherm,30);
axis([1 5.5 -4 3]);
legend = colorbar;
clabel(c,h);
xlabel('East-West Distance (km)');
ylabel('North-South Distance (km)');
ylabel(legend, 'Temperature (Degrees Celsius)');
title('Isotherms at Same Elevation as Point(4,-0.3)');

%2C
NWUnit2D = [-1/sqrt(2), 1/sqrt(2)];
SWUnit2D = [-1/sqrt(2), -1/sqrt(2)];
gradientAtPoint = [dzdx(4,-0.3), dzdy(4,-0.3)];
heightChangeNW = double(dot(NWUnit2D, gradientAtPoint));
heightChangeSW = double(dot(SWUnit2D, gradientAtPoint));

tangentPlaneZ(x,y) = dzdx(4,-0.3)*(x-4)+dzdy(4,-0.3)*(y+0.3);
NWZ = tangentPlaneZ(-1, 1);
SWZ = tangentPlaneZ(-1, -1);
NWUnit3D = [-1, 1,NWZ]/sqrt(1+1+NWZ^2);
SWUnit3D = [-1, -1,SWZ]/sqrt(1+1+SWZ^2);
gradientT = gradient(T);
gradientTAtPoint = gradientT(4,-0.3,heightAtPoint);
tempChangeNW = double(dot(gradientTAtPoint, NWUnit3D));
tempChangeSW = double(dot(gradientTAtPoint, SWUnit3D));

%2D
zz = 0.00125 * exp(-((0.5*yy.^2)+(xx-3).^2)).*(sin(2*xx)+2*sin(0.75*(0.5*yy-2).^2)).*(16*xx+64*xx.^2 + yy.^2);
ColorT = (-0.1*(zz.^2))+17*exp(-0.1*((0.1*xx-2)-(0.05*yy-1).^2-(zz-1).^2))-10;
figure(4)
surf(xx,yy,zz,ColorT);
legend=colorbar;
xlabel('East-West Distance (km)');
ylabel('North-South Distance (km)');
ylabel(legend, 'Temperature (Degrees Celsius)');
zlabel('Elevation (km)');
title('Terrain Plot with Temperature Colour Map');

%2E
figure(5)
surf(xx,yy,ColorT)
xlabel('East-West Distance (km)');
ylabel('North-South Distance (km)');
zlabel('Temperature (Degrees Celsius');
title('Temperature Plot');

%2F
C=0.00125.*exp(-((x-3).^2+0.5.*y.^2)).*(sin(2.*x)+2.*sin(0.75*(0.5*y-2).^2)).*(16.*x+64.*x.^2+y.^2)-zzz;
T(x,y) = 0.1.*(0.00125.*exp(-((x-3).^2+0.5.*y.^2)).*(sin(2.*x)+2.*sin(0.75*(0.5*y-2).^2)).*(16.*x+64.*x.^2+y.^2)).^2+17*exp(-0.1.*((0.1.*x-2)-(0.05.*y-1).^2-((0.00125.*exp(-((x-3).^2+0.5.*y.^2)).*(sin(2.*x)+2.*sin(0.75*(0.5*y-2).^2)).*(16.*x+64.*x.^2+y.^2))-1).^2))-10;
Lagrange = T - l*C;
dLdx = diff(Lagrange, x);
dLdy = diff(Lagrange, y);
dLdz = diff(Lagrange, zzz);
dLdl = diff(Lagrange, l);

lagVF = [dLdx, dLdy, dLdz, dLdl];
vars= [x y zzz l];
LagGuess = [3, -0.8, -1.2, 0];
lagVFmatlab = matlabFunction(lagVF,'Vars',{vars});
pointMaxTemp = fsolve(lagVFmatlab, LagGuess);


function [determination]= deriv2Test(A, B, C, cp, x, y)
    A1 = double(subs(A, [x y], cp));
    B1 = double(subs(B, [x y], cp));
    C1 = double(subs(C, [x y], cp));
    D = (B1).^2 - (A1).*(C1);
    if D<0
        if A1<0
            determination = 'Local Max';
        elseif A1>0
            determination = 'Local Min';
        end
    elseif D>0
        determination = 'Saddlepoint';
    else
        determination = 'No conclusion';
    end
end