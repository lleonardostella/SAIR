x=[1,2,4,6,8]';
y=[100,140,160,170,175].';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(1,8,50);
plot(x,y,'o',xx,f0(xx),'r-');