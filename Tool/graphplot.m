clear all
close all
clc

x1=[5.5,5.3,5,4.8,4.6,4.4,4.2,4,3.8,3.6,3.4];
y1=[98,87,79,67,61,53,45,36,31,26,20];

x2=[5.5,5.3,5,4.8,4.6,4.4,4.2,4,3.7];
y2=[59,54,48,43,37,31,26,22,18];

size=25;



%%
p1=polyfit(x1,y1,1);
x1_2=1:1:10;
y1_2=polyval(p1,x1_2);

%%
p2=polyfit(x2,y2,1);
x2_2=1:1:10;
y2_2=polyval(p2,x2_2);
scatter(x1,y1,size,"filled",'black')
hold on
plot(x1_2,y1_2,'black')
scatter(x2,y2,25,"filled",'red')
hold on
plot(x2_2,y2_2,'red')
xlim([2.5 6]);
ylim([0 110]);
xlabel('Pump power (W)')
ylabel('Output power (mW)')
a1=legend('AP vortex','Theoretical fit','RP vortex','Theoretical fit');
%a1=legend('AP vortex','1.4% Slope efficiency','RP vortex','1.0% Slope efficiency');
%saveas(jpg,AP)