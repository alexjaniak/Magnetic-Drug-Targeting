%{
Script to model & graph x0 vs P
%}


T = readtable('data/x0vsprob02.dat');
T = removevars(T,{'Var2'});
T.Properties.VariableNames = {'x0','P'};

graph_x = linspace(-0.065,0.025,1000);
graph_y = fittedmodelpspline(graph_x);

figure(1);
hold on
scatter(T.x0,T.P,'.');
plot(graph_x,graph_y,'color',[0.9290 0.6940 0.1250],'linewidth',2);
hold off

title('Initial Position ( X_0 ) vs Capture Probability');
xlabel('X_0');
ylabel('P');
savefig('graphs/x0vsP01.fig');

disp(max(graph_y)) %max
disp(graph_x(find(graph_y == max(graph_y))));

%disp(fittedmodelpspline);
%disp(goodnesspspline);
%disp(outputpspline);