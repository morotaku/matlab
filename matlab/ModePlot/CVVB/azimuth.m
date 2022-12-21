function [Ex,Ey]=azimuth(E,l,phi,thita)
ex=-sin(l.*phi).*E;
ey=cos(l.*phi).*E;
[Ex,Ey]=polarizer(thita,ex,ey);
end