clear all
close all
clc
format shortEng
%%% Topological charge
l=0;
%%% Radial index
p=1;

%%% Parameters(µm)
w=10000; %Beam waist
lam=1.06; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size

f=20000; %Focus length
s=200000; %Distance from the input plane to lens
G0=1;
i=1j;

%ABCD matrix
A=-z/f;
B=-(z/f)*s+f+z;
C=-1/f;
D=1-s/f;

%%%%%%%%%%%%%%%%Intinsity plot%%%%%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N1=1001;
L=20; %Display lange
X1=linspace(-L,L,N1);
Y1=linspace(-L,L,N1);
[x1,y1]=meshgrid(X1,Y1);
[phi2,r1] = cart2pol(x1,y1);

%%%Mode plot
HGB=hgb(p,l,r1,z,w,k,A,B,C,D,G0);



%%%%%%%%%%%%%%%%Polarization plot%%%%%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=30;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);
[phi2,r2] = cart2pol(x2,y2);

HGB2=hgb(p,l,r2,z,w,k,A,B,C,D,G0);
max(HGB2(:))

%%% angular frequency ω
one_q=ones(size(phi2,1));
o_q=one_q.*(pi/20);

for t = 0:200
    imagesc(X1,Y1,HGB);
    hold on
    %%%HGB
    E1_q=exp(1j.*(l.*phi2+t*o_q)).*HGB2;  %ℓ=1
    %%%|LHC> ℓ=1
    ex1_q=real(E1_q);
    ey1_q=real(-1j.*E1_q);

    ex_q=ex1_q;
    ey_q=ey1_q;


    q=quiver(x2,y2,ex_q,ey_q);
    q.LineWidth=1;
    q.Color='red';
    colormap("gray")
    shading interp; lighting phong; view(2); axis equal; axis tight; axis off;
    hold off
    pause(0.05)
end


