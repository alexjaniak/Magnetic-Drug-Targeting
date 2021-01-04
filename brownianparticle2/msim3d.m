
%{
3d (xyz) run file for creating brownian motion simulations
using the fbrownian.m function and varying initial conditions
%}

%Setup for .dat file
CREATE_TITLE = false; 
a_w = 'w'; % Either 'a' for append or 'w' for write

%Create data
num_epochs = 1000; %Number of simulations (intial positions)
a_prob = zeros(1,num_epochs); %Array to hold probabilities

%{
%For equal increments of Initial Position
if num_epochs == 1
    X = [-0.01];
else
    X = logspace(-11, -6, num_epochs);
end
%}

% Generating random points in a cylinder
radius = 1e-4;
yc = 0;
zc = 0;
% Engine
theta = rand(1,num_epochs)*(2*pi);
r = sqrt(rand(1,num_epochs))*radius;
y = yc + r.*cos(theta);
z = zc + r.*sin(theta);
x = (0.01+0.04)*rand(1,num_epochs)-0.04;
scatter3(x,y,z);


%Simulations
    for i = 1:num_epochs
            fprintf("Simulation %d\n",i);
            seed = randi(100);
            prob = fbrownian(x(i),y(i),z(i),seed);
            a_prob(i) = prob;
    end

%Append or write everything to a .dat file
if true
    fid=fopen('xyzvsprob01.dat',a_w);
    if CREATE_TITLE
        fprintf(fid,'x0,y0,z0,P\n');
    end
    
    if a_w == 'a'
        fprintf(fid,'\n');
    end
    
    for i = 1:num_epochs 
        fprintf(fid,'%12.5g %s %12.5g %s %12.5g %s %12.5g\n',z(i),', ',y(i),',',z(i),',',a_prob(i));
    end
    fclose(fid);
    
end
