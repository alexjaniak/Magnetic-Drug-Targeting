%{
Script to generate 3d graphs used in presentation
%}

% Generate graph of random points in a cylinder
% Data
n = 2000;
radius = 1e-4;
yc = 0;
zc = 0;
% Engine
theta = rand(1,n)*(2*pi);
r = sqrt(rand(1,n))*radius;
y = yc + r.*cos(theta);
z = zc + r.*sin(theta);
x = (0.01+0.04)*rand(1,n)-0.04;
% Check
figure(1);
scatter3(x,y,z,'.');
title('Initial Positions ( X_0  Y_0  Z_0 )');
xlabel('X_0');
ylabel('Y_0');
zlabel('Z_0');
savefig('graphs/sample-1.fig');

%Generate graph for data with a colormap
T = readtable('data/xyzvsP01.dat');
T = removevars(T,{'Var2','Var4','Var6'});
T.Properties.VariableNames = {'x0','y0','z0','P'};

figure(2);
scatter3(T.x0,T.y0,T.z0,5,T.P);
title('Initial Positions ( X_0  Y_0  Z_0 )');
xlabel('X_0');
ylabel('Y_0');
zlabel('Z_0');
colormap(hot);
cb = colorbar;
title(cb,'Capture Probability');
savefig('graphs/xyz01-1.fig');

%Generate graph for predicted Capture Probability using the Gaussian RQ
%model
gpr = gaussianrq(T);
yfit = gpr.predictFcn(T);
figure(3);
scatter3(T.x0,T.y0,T.z0,5,yfit);
title('Initial Positions ( X_0  Y_0  Z_0 )');
xlabel('X_0');
ylabel('Y_0');
zlabel('Z_0');
colormap(hot);
cb = colorbar;
title(cb,'Predicted Capture Probability');
savefig('graphs/xyz01-2.fig');

%Generate graph for a more complete look at the model
% Data
n = 10000;
radius = 1e-4;
yc = 0;
zc = 0;
% Engine
theta = rand(1,n)*(2*pi);
r = sqrt(rand(1,n))*radius;
y = transpose(yc + r.*cos(theta));
z = transpose(zc + r.*sin(theta));
x = transpose((0.01+0.04)*rand(1,n)-0.04);
p = zeros(n,1);
T2 = table(x,y,z,p);
T2.Properties.VariableNames = {'x0','y0','z0','P'};

yfit2 = gpr.predictFcn(T2);
figure(4);
scatter3(T2.x0,T2.y0,T2.z0,5,yfit2);
title('Initial Positions ( X_0  Y_0  Z_0 )');
xlabel('X_0');
ylabel('Y_0');
zlabel('Z_0');
colormap(hot);
cb = colorbar;
title(cb,'Predicted Capture Probability For Random');
savefig('graphs/xyz01-3.fig');

%Generate a graph with only >= 0.5 capture probability
hp_filter = [];
lp_filter = [];
for i = 1:size(T2.x0)
    if yfit2(i) >= 0.5
        a = [T2.x0(i) T2.y0(i) T2.z0(i) yfit2(i)];
        hp_filter = [hp_filter;a];
    else 
        b = [T2.x0(i) T2.y0(i) T2.z0(i) yfit2(i)];
        lp_filter = [lp_filter;a];
    end 
end
hpx = hp_filter(:,1);
hpy = hp_filter(:,2);
hpz = hp_filter(:,3);
hpp = hp_filter(:,4);

lpx = lp_filter(:,1);
lpy = lp_filter(:,2);
lpz = lp_filter(:,3);
lpp = lp_filter(:,4);


figure(5);
scatter3(hpx,hpy,hpz,'.');
title({'Initial Positions ( X_0  Y_0  Z_0 )';'Capture Probability >= 0.5'});
xlabel('X_0');
ylabel('Y_0');
zlabel('Z_0');
savefig('graphs/xyz01-4.fig');

%Generate a slice of the graph
hp_slice = [];
for i = 1:size(hpx)
    if (hpx(i) >= -0.01) & (hpx(i)<= -0.005)
        a = [hpx(i) hpy(i) hpz(i) hpp(i)];
        hp_slice = [hp_slice;a];
    end 
end

shpx = hp_slice(:,1);
shpy = hp_slice(:,2);
shpz = hp_slice(:,3);
shpp = hp_slice(:,4);

figure(6);
scatter3(hpx,hpy,hpz,'.');
hold on
scatter3(shpx,shpy,shpz,'r','.');
hold off
title({'Initial Positions ( X_0  Y_0  Z_0 )';'Capture Probability >= 0.5';'X_0 = [-0.01,-0.005]'});
xlabel('X_0');
ylabel('Y_0');
zlabel('Z_0');
legend('Initial Positions','Initial positions in the slice [-0.01,-0.005]');
savefig('graphs/xyz01-5.fig');

%Generate a slice of the low prob graph
lp_slice = [];
for i = 1:size(lpx)
    if (lpx(i) >= -0.01) & (lpx(i)<= -0.005)
        a = [lpx(i) lpy(i) lpz(i) lpp(i)];
        lp_slice = [lp_slice;a];
    end 
end

slpx = lp_slice(:,1);
slpy = lp_slice(:,2);
slpz = lp_slice(:,3);
slpp = lp_slice(:,4);

figure(7);
scatter3(slpx,slpy,slpz,'.');
savefig('graphs/xyz01-6.fig');

for i = 1:size(shpp)
    if shpp(i) > 1
        shpp(i) = 1;
    elseif shpp(i) < 0
        shpp(i) = 0;
    end
end

for i = 1:size(slpp)
    if slpp(i) > 1
        slpp(i) = 1;
    elseif slpp(i) < 0
        slpp(i) = 0;
    end
end

% calculate average capture probability in the slice
av_p = (sum(shpp)+sum(slpp))/(size(shpp,1)+size(slpp,1));
hpratio = size(shpx)/(size(shpx)+size(slpx));
disp(av_p);