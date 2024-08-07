for t = 0:0
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g1)+t.*o_q(g1))).*lg1_1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g1)+t.*o_q(g1))).*lg2_1; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_1=real(E1_1);
    ey1_1=real(-1j.*E1_1);
    %%%|RHC> ℓ=-1
    ex2_1=real(E2_1);
    ey2_1=real(1j.*E2_1);
    
    ex_1=ex1_1+ex2_1;
    ey_1=ey1_1+ey2_1;
    [Ex_1,Ey_1]=polarizer(angle,ex_1,ey_1);
    Ex1=Ex_1+x_1;
    Ey1=Ey_1+y_1;
    axis equal; axis off;
    
    scatter(Ex1,Ey1,'.',map)
    %hold off
    pause(0.005)
end
hold on
for t = 0:60
    %%%North pole of Poincare sphere
    E1_2=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(lg1)+t.*o_q(lg1))).*lg1_2;  %ℓ=1
    %%%South pole
    E2_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(lg1)+t.*o_q(lg1))).*lg2_2; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_2=real(E1_2);
    ey1_2=real(-1j.*E1_2);
    %%%|RHC> ℓ=-1
    ex2_2=real(E2_2);
    ey2_2=real(1j.*E2_2);
    
    ex_2=ex1_2+ex2_2;
    ey_2=ey1_2+ey2_2;
    [Ex_2,Ey_2]=polarizer(angle,ex_2,ey_2);
    Ex2=Ex_2+x_2;
    Ey2=Ey_2+y_2;
    axis equal; axis off;
    
    scatter(Ex2,Ey2,'.','blue')
    %hold off
    pause(0.005)
end