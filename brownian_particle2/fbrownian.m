function prob = fbrownian(xinit,yinit,zinit,seed)
        global m_dpl chi mu0 viscosity rad_pipe rad_particle xdipole vmax vosc omega;

        %
        mu0=4e-7*pi;            
        kb=1.3806503e-23;       %boltzmann constant
        temp=312;               %temperature (kelvin)
        %
        %mm=1e-3;                %magnitude of dipole moment
        %magnetsize=1.e-3;       %magnet size
        mm=1065.52;
        magnetsize=0.0337218;
        rad_particle=1.e-10;    %particle radius(m)
        %rad_particle = rp;
        rho_particle=1.e3;      %particle density(kg/m^3)
        rad_pipe=1e-4;          %pipe radius(m)
        %
        lx=100*rad_pipe;         %particles are released at x=-lx
        tumor_radius=1.e-2;     %tumor radius
        %x0=[-lx 0 -0.0*rad_pipe];           %initial position of particle
        x0=[xinit yinit zinit];   
        %
        viscosity=4.e-3;        %blood viscosity (pa s)
        chi=0.2;                %susceptibility
        vmax=2/pi*1.e-2;        %average axial velocity
        vosc=0;                 %amplitude of oscilaltion
        omega=18;           
        %
        output_trajectory=true;
        %
        t0=0;               %starting time
        tf=10;               %ending time
        n=50000;            %number of time steps to be simulated
        numsim=100;         %number of random simulations
        dt=(tf-t0)/n;       %time step
        %
        xdipole=[0 0 -rad_pipe-magnetsize];
        theta=0;
        %


        %tf=tf+dt;
        %n=n+1;

        % setup bigtheta
        bigtheta=zeros(1,100);  %100 is a big enough number to store all the parameters
        dim=3;                  %dimension of the brownian motion, no need to change
        ii=1;
        %
        bigtheta(ii:ii+dim-1)=x0;	% initial position of particle
        ii=ii+dim;
        bigtheta(ii:ii+dim-1)=xdipole;	% position of target
        ii=ii+dim;
        bigtheta(ii:ii+dim-1)=mm*[sind(theta) 0 cosd(theta)];	% target dipole moment
        ii=ii+dim;
        %
        bigtheta(ii)=rad_particle;	% particle radius
        ii=ii+1;
        bigtheta(ii)=rho_particle;	% particle density
        ii=ii+1;
        bigtheta(ii)=mu0;	% vacuum permeability
        ii=ii+1;
        bigtheta(ii)=kb;	% boltzmann constant
        ii=ii+1;
        bigtheta(ii)=temp;	% temperature (kelvin)
        ii=ii+1;
        bigtheta(ii)=chi;	% suspectibility of particle
        ii=ii+1;
        bigtheta(ii)=viscosity;	% blood viscosity
        ii=ii+1;
        bigtheta(ii)=rad_pipe;	% blood vessel radius
        ii=ii+1;
        bigtheta(ii)=vmax;	% blood velocity at vessel center
        ii=ii+1;
        bigtheta(ii)=vosc;	% amplitude of velocity oscillation
        ii=ii+1;
        bigtheta(ii)=omega;	% angular frequency of oscillation
        ii=ii+1;
        %

        problem='Brownian';
        owntime=[t0:dt:tf];
        numdepvars=dim;
        sdetype='ito';
        %seed = [];
        %seed=1;

        %the structure of xhat: 
        % xhat(i,(j-1)*dim+1:(j+1)*dim) stores the position of particle at t(i), in
        % the jth simulation.
        xhat=SDE_euler(bigtheta,problem,owntime,numdepvars,numsim,sdetype,seed);

        %get the probability of hitting target

        idx=ones(1,numsim)*(n+1);

        x=zeros(n+1,numsim);
        y=zeros(n+1,numsim);
        z=zeros(n+1,numsim);

        x(1:n+1,1:numsim)=xhat(1:n+1,1:dim:dim*numsim);
        y(1:n+1,1:numsim)=xhat(1:n+1,2:dim:dim*numsim);
        z(1:n+1,1:numsim)=xhat(1:n+1,3:dim:dim*numsim);
        %

        rr(1:n+1,1:numsim)=rad_pipe-sqrt(y(1:n+1,1:numsim).^2+ ...
            z(1:n+1,1:numsim).^2);

        hit=0;
        for j=1:numsim
            for i=1:n+1
                if(rr(i,j)<=0)  %rr(i-1,j)>0, r(i,j)<=0
                    idx(j)=i;
                    
                    c1=abs(rr(i,j))/(abs(rr(i-1,j))+abs(rr(i,j)));
                    c2=1-c1;
                    x(i,j)=c1*x(i-1,j)+c2*x(i,j);   %reset the particle position
                    y(i,j)=c1*y(i-1,j)+c2*y(i,j);   %to the pipe wall
                    z(i,j)=c1*z(i-1,j)+c2*z(i,j);
                    if(abs(x(i,j))<tumor_radius)
                        hit=hit+1;
                    end
                    break;
                end
            end
        end





        probability=hit/numsim;   %the number of hits divided total simulations

        fprintf('probability of hitting target = %g\n',probability); 
        %visualization

        %if(output_trajectory)
        %    fid=fopen('trajectory.dat','a');
        %    fprintf(fid,'title="particle trajectories" \n');
        %    fprintf(fid,'variables= "x" "y" "z" "t" \n');
        %    for j=1:numsim
        %        fprintf(fid,'zone t="%4.4d", i=%d, f=point\n', ...
        %                 j,idx(j));
        %        for i=1:idx(j);
        %            fprintf(fid,'%12.5g %12.5g %12.5g %12.5g\n',...
        %                x(i,j),y(i,j),z(i,j),dt*(i-1));
        %        end
        %    end
        %    fclose(fid);
        %end
           
        %for j=1:numsim
        %        for i=idx(j)+1:n+1
                % remove the portion of trajectories after hitting wall
        %        x(i,j)=x(idx(j),j);
        %        y(i,j)=y(idx(j),j);
        %        z(i,j)=z(idx(j),j);
        %        end
        %end
        %plot(xhat(:,dim*([1:numsim]-1)+1),xhat(:,dim*[1:numsim]));

        %figure(1);
        %plot(x(:,1:numsim),z(:,1:numsim));
        %title(['particle trajectories (numsim=',num2str(numsim,'%d'),')']);
        %axis([-lx*2, lx*6, -rad_pipe, -rad_pipe*0.8]);
        %xlabel('x(m)');
        %ylabel('z(m)');

        %figure(2);
        %plot(x(:,1:numsim),y(:,1:numsim));
        %title(['particle trajectories (numsim=',num2str(numsim,'%d'),')']);
        %axis([-lx*2, lx*6, -rad_pipe*0.1, rad_pipe*0.1]);
        %xlabel('x(m)');
        %ylabel('y(m)');

        %figure(3);
        %plot(y(:,1:numsim),z(:,1:numsim));
        %title(['particle trajectories (numsim=',num2str(numsim,'%d'),')']);
        %axis([-rad_pipe*0.1,rad_pipe*0.1, -rad_pipe, -rad_pipe*0.8]);
        %xlabel('y(m)');
        %ylabel('z(m)');
        prob = probability;
end
