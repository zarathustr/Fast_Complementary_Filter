clear all;
clc;

addpath('quaternion_library');

conv=diag([-1 1 -1]);
data=load('flat_motion2.txt');
quaternion_FCF=load('output.txt');

euler_FCF=quatern2euler(quaternConj(quaternion_FCF))*180/pi;

time=1/500*(1:length(data(:,1)));

figure(1);
subplot(3,1,1);
plot(time,-data(:,1),time,euler_FCF(:,1));
legend('Reference','Complementary');
xlabel('Time (s)');
ylabel('Roll (deg)');

subplot(3,1,2);
plot(time,data(:,2),time,euler_FCF(:,2));
legend('Reference','Complementary');
xlabel('Time (s)');
ylabel('Pitch (deg)');

subplot(3,1,3);
plot(time,-data(:,3)+180,time,euler_FCF(:,3));
legend('Reference','Complementary');
xlabel('Time (s)');
ylabel('Yaw (deg)');