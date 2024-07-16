%% Making Plots with Trendlines
load('Hcsf_Age.txt')
load('Hcsf_Age_2.txt')
load('Hcsf_Age_3.txt')
load('DGCA3_DR.txt')
load('CA1CA2_DR.txt')
load('SUB_DR.txt')
load('ERC_DR.txt')

x1 = Hcsf_Age;
x2 = Hcsf_Age_2;
x3 = Hcsf_Age_3;
y1 = DGCA3_DR;
y2 = CA1CA2_DR;
y3 = SUB_DR;
y4 = ERC_DR;

figure;
p = polyfit(x,y1,1); 
f = polyval(p,x); 
plot(x,y1,'o',x,f,'-') 
legend('data','linear fit') 
hold on;

%figure;
p = polyfit(x,y2,1); 
f = polyval(p,x); 
plot(x,y2,'o',x,f,'-') 
legend('data','linear fit') 
hold on;

%figure;
p = polyfit(x,y3,1); 
f = polyval(p,x); 
plot(x,y3,'o',x,f,'-') 
legend('data','linear fit') 
hold on;

%figure;
p = polyfit(x,y4,1); 
f = polyval(p,x); 
plot(x,y4,'o',x,f,'-') 
legend('data','linear fit') 
hold off;

%% Quadratic Fit
figure;
p = polyfit(x1,y1,2); 
f = polyval(p,x1); 
plot(x1,y1,'o',x1,f,'-') 
legend('data','linear fit') 
hold on;

%figure;
p = polyfit(x2,y2,2); 
f = polyval(p,x2); 
plot(x2,y2,'o',x2,f,'-') 
legend('data','linear fit') 
hold on;

%figure;
p = polyfit(x1,y3,2); 
f = polyval(p,x1); 
plot(x1,y3,'o',x1,f,'-') 
legend('data','linear fit') 
hold on;

%figure;
p = polyfit(x3,y4,2); 
f = polyval(p,x3); 
plot(x3,y4,'o',x3,f,'-') 
legend('data','linear fit') 
hold off;