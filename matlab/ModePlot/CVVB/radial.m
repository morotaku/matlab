function [Ex,Ey]=radial(E,l,phi,thita)
ex=cos(l.*phi).*E;
ey=sin(l.*phi).*E;
[Ex,Ey]=polarizer(thita,ex,ey);
end

