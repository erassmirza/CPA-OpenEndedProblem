clear, clc, close

%%                              QUESTION 01
% EFI Engine given data
EFI_rpm = [1522 1691 1803 1898 2000];
EFI_bhp = [4.029 4.916 5.644 5.800 6.025];

% DIESEL Engine given data
DIESEL_rpm = [1517 1707 1922 2009 2224];
DIESEL_bhp = [1.023 1.207 1.437 1.809 1.744];

% PETROL Engine given data
PETROL_rpm = [994 1097 1200 1309 1395];
PETROL_bhp = [0.823 0.785 1.853 1.911 1.841];

% For curve fitting using polyfit function
p1 = polyfit(EFI_rpm,EFI_bhp,2);
p2 = polyfit(DIESEL_rpm,DIESEL_bhp,2);
p3 = polyfit(PETROL_rpm,PETROL_bhp,2);

% New bhp values after curve fitting
eb_p = polyval(p1,EFI_rpm);                 % EFI
db_p = polyval(p2,DIESEL_rpm);              % DIESEL
pb_p = polyval(p3,PETROL_rpm);              % PETROL

% Error estimation
errorEFI = (eb_p - EFI_bhp)*100             % EFI
errorDIESEL = (db_p - DIESEL_bhp)*100       % DIESEL
errorPETROL = (pb_p - PETROL_bhp)*100       % PETROL

% Prediction
efi_at_1500 = polyval(p1,1500)              % EFI
diesel_at_1500 = polyval(p2,1500)           % DIESEL
petrol_at_1500 = polyval(p3,1500)           % PETROL

% Graph
figure
subplot(3,2,1); hold on;
plot(EFI_rpm,EFI_bhp,'rx',EFI_rpm,eb_p);
plot(1500,efi_at_1500,'o','MarkerSize',6,'Linewidth',2); hold off;
xlabel('Speed(rpm)'); ylabel('Power (bhp)');
title('EFI Engine (Quardatic)'); grid;
legend('Given Data','Model','BHP at 1500');

subplot(3,2,3); hold on;
plot(DIESEL_rpm,DIESEL_bhp,'rx',DIESEL_rpm,db_p)
plot(1500,diesel_at_1500,'o','MarkerSize',6,'Linewidth',2); hold off;
xlabel('Speed(rpm)'); ylabel('Power (bhp)');
title('DIESEL Engine (Quardatic)'); grid;
legend('Given Data','Model','BHP at 1500');

subplot(3,2,5); hold on;
plot(PETROL_rpm,PETROL_bhp,'rx',PETROL_rpm,pb_p)
plot(1500,petrol_at_1500,'o','MarkerSize',6,'Linewidth',2); hold off;
xlabel('Speed(rpm)'); ylabel('Power (bhp)');
title('PETROL Engine (Quardatic)'); grid;
legend('Given Data','Model','BHP at 1500');

% Linear Relation
% For curve fitting using polyfit function
pl1 = polyfit(EFI_rpm,EFI_bhp,1);
pl2 = polyfit(DIESEL_rpm,DIESEL_bhp,1);
pl3 = polyfit(PETROL_rpm,PETROL_bhp,1);

% New bhp values after curve fitting
eb_lp = polyval(pl1,EFI_rpm);                 % EFI
db_lp = polyval(pl2,DIESEL_rpm);              % DIESEL
pb_lp = polyval(pl3,PETROL_rpm);              % PETROL

% Graph
subplot(3,2,2);
plot(EFI_rpm,EFI_bhp,'rx',EFI_rpm,eb_lp);
xlabel('Speed(rpm)'); ylabel('Power (bhp)');
title('EFI Engine (Linear)'); grid;
legend('Given Data','Model');

subplot(3,2,4);
plot(DIESEL_rpm,DIESEL_bhp,'rx',DIESEL_rpm,db_lp)
xlabel('Speed(rpm)'); ylabel('Power (bhp)');
title('DIESEL Engine (Linear)'); grid;
legend('Given Data','Model');

subplot(3,2,6);
plot(PETROL_rpm,PETROL_bhp,'rx',PETROL_rpm,pb_lp)
xlabel('Speed(rpm)'); ylabel('Power (bhp)');
title('PETROL Engine (Linear)'); grid;
legend('Given Data','Model');

%%                              QUESTION 02
% log(IT) = - beta*L + log(Io*(1-R)^2) 

% Intensity of incident beam (I_B) = 5 W/m2
I_B = 5;

% Length of a transparent solid
L = [0.5 1.2 1.7 2.2 4.5 6.0];

% Transmitted intensity
I_T = [4.2 4.0 3.8 3.6 2.9 2.5];

% For curve fitting using polyfit function
p = polyfit(L,log(I_T),1);

% Absorption coefficient
beta = - p(1);
fprintf('\tAbsoption coefficient = %f\n',beta);

% Fraction of light which is reflected at the interface
% p(2) = log(I_B*(1-R)^2)
R = 1 - sqrt(exp(p(2)) / I_B);
fprintf('\tFraction of light = %f\n',R);

% Index of refraction for the transparent solid
% R = sqrt((n-1)/(n+1))
n = 2/(1-sqrt(R))-1;
fprintf('\tIndex of refraction = %f\n',n);

% Graph
Lp = linspace(0.5,6,100);
F = @ (x) I_B*(1-R)^2*exp(-beta*x);
I_Tp = F(Lp);
figure
plot(L,I_T,'ro',Lp,I_Tp)
xlabel('L (cm)'); ylabel('I_T (W/m^2)')
legend('Given Data','Model'); grid;
title('Intensity of light through transparent solid')

% Absorption coefficient & Index of refraction Values writen on graph
str = {'n = 1.6242','\beta = 0.0956'};
text(5.3,4.1,str)