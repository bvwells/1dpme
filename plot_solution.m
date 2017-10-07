clear all;

load 'variables.m';

m = variables(1);
t_init = variables(4);
output_t = variables(5);
numberReports = variables(6);
dt = (output_t-t_init)/numberReports;

figure(1)
title(strcat('Porous Medium Equation with m=', num2str(m)))
legendArray = []

for report = 0:numberReports 
   
   filename = strcat('solution', sprintf('%03d',[report]), '.m');

   u = load(filename);
   
   time = t_init + dt*report
 
   plot(u(:,1),u(:,2));

   legendArray = [legendArray ; strcat('t=', num2str(time))];
   hold on;
 end

 axis square;
 axis tight;
 xlabel('x');
 ylabel('u');
 legend(legendArray);