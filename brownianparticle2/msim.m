
%{
Main run file for creating brownian motion simulations
using the fbrownian.m function and varying initial conditions
%}

%Setup for .dat file
CREATE_TITLE = false; 
a_w = 'a'; % Either 'a' for append or 'w' for write

%Create data
num_epochs = 10; %Number of simulations (intial positions)
a_prob = zeros(1,num_epochs); %Array to hold probabilities

%{
For equal increments of Initial Position
if num_epochs == 1
    X = [-0.01];
else
    X = logspace(-11, -6, num_epochs);
end
%}

x = (-0.06+0.065)*rand(1,num_epochs)-0.065;
%Simulations
    for i = 1:num_epochs
            fprintf("Simulation %d\n",i);
            seed = randi(100);
            prob = fbrownian(x(i),0,0,seed);
            a_prob(i) = prob;
    end

%Append or write everything to a .dat file
if true
    fid=fopen('x0vsprob02.dat',a_w);
    if CREATE_TITLE
        fprintf(fid,'x0,P\n');
    end
    
    if a_w == 'a'
        fprintf(fid,'\n');
    end
    
    for i = 1:num_epochs 
        fprintf(fid,'%12.5g %s %12.5g\n',x(i),',',a_prob(i));
    end
    fclose(fid);
    
end
