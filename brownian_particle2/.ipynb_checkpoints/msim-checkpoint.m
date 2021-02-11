
%{
2019-07-10
Main run file for creating brownian motion simulations
using the fbrownian.m function and varying initial position of
clusters (X)
%}

%Setup for .dat file
CREATE_TITLE = true; 
a_w = 'w'; % Either 'a' for append or 'w' for write

%Create data
num_epochs = 100; %Number of simulations (intial positions)
a_prob = zeros(1,num_epochs); %Array to hold probabilities
X = zeros(1, num_epochs);

%if num_epochs == 1;
%    X = [-0.01];
%else
%    X = linspace(0, 0.02, num_epochs);
%end

%Simulations
tic;
    for i = 1:num_epochs
            fprintf("Simulation %d\n",i);
            xinit = -0.01*rand;
            seed = randi(100);
            X(i) = xinit;
            prob = fbrownian(X(i),seed);
            a_prob(i) = prob;
    end
time = toc;
disp(time);

%Append everything to a .dat file
if true
    fid=fopen('x0vsprob01.dat',a_w);
    if CREATE_TITLE
        
        fprintf(fid,'x0,P\n');
    
    end
    for i = 1:num_epochs
        
        fprintf(fid,'%12.5g %s %12.5g\n',X(i),', ',a_prob(i));
        
    end
    fclose(fid);
    
end
