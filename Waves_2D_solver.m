function [p,f,FB,vFB,FBdelx,v_array,vvv_array] = Waves_2D_solver(p0,f0,v0,vvv0,FB0,vFB0,FBdelx0,r,x,dx,dt,tfinal,p_small,p_cr,gr,w,va,ww,lam0)
% authors: Mariya Savinov, Alex Mogilner

% Solves the advection-reaction equation on 2D domain

% Inputs:
%     p0         :  initial density
%     f0         :  initial 'gas' density regrowing behind previous wave
%     v0         :  initial velocity
%     vvv0       :  initial velocity (connected density)
%     FB0        :  initial free boundary location
%     vFB0       :  initial free boundary velocity
%     r          :  radius of exclusion zone
%     x          :  space array from r
%     dx         :  spatial step
%     dt         :  time step
%     tfinal     :  final time
%     p_small    :  gas <--> connected threshold
%     p_cr       :  connected <--> contractile threshold
%     gr,w,ww    :  density-dependent assembly parameters
%     va         :  flux parameter
%     lam0       :  contractile strength


% Outputs:
%     p          :  density at all times
%     f          :  gas density at all times
%     FB         :  free boundary location at all times
%     vFB        :  free boundary velocity at all times
%     FBdelx     :  free boundary growth at final time
%     v_array    :  velocity at all times
%     vvv_array  :  velocity of connected network at all times



global makeanimation

% number of time steps
nsteps = round(tfinal / dt);
if abs(dt*nsteps - tfinal) > 1e-5
   % The last step won't go exactly to tfinal.
   disp(' ')
   disp(sprintf('WARNING *** dt does not divide tfinal, dt = %9.5e',dt))
   disp(' ')
end

%----------------------------Initialize--------------------------%

% number of x-points
Nx = length(x);

% initialize vectors
p = zeros(nsteps+1,Nx);             % to store density at all times
v_array = zeros(nsteps+1,Nx);       % to store velocity at all times
vvv_array = zeros(nsteps+1,Nx);     % to store velocity for connected portion at all times
f = zeros(nsteps+1,Nx);             % to store 'gas' density (unconnected portion) at all times
FB = zeros(nsteps+1,1);             % to store free boundary location
vFB = zeros(nsteps+1,1);            % to store free boundary velocity

% Initial Conditions
p(1,:) = p0;                        % initial density   % e.g. = (p_cr + 0.1)*ones(Nx,1)
v = v0;                             % initial velocity
v_array(1,:) = v;                   % store initial velocity
tn = 0;                             % start time

% gas network density
f(1,:) = f0;

% velocity corresponding to only the connected portion of the network
vvv = vvv0;                         % initial velocity of connected portion
vvv_array(1,:) = vvv;               % store initial vvv

% free boundary position and velocity
FB(1) = FB0;
vFB(1) = vFB0;


% free boundary growth per timestep
FBdx = va*dt;
% free boundary growth tracking
FBdelx = FBdelx0;

if makeanimation==1
    ss=VideoWriter('test_waves_2D.avi');open(ss);
end

% assembly dependent assembly + disassembly
% = a(p) - p = w + gr*(1-exp(-p/ww)) - p
alphabeta = @(xp) w+gr*(1-exp(-xp/ww))-xp;

    
%----------------------------MainLoop--------------------------%
for n = 1:nsteps 
    % current time
    tnp = tn + dt;     % = t_{n+1}
    
    % update the density
    p(n+1,1) = p(n,1) + lam0*dt/dx*(v(2)*p(n,2)) + ...  % first point
        + lam0*dt*(v(1)*p(n,1))/(r+x(1)) + ...
       dt * (sign(p(n,1)-p_small)+1)/2 * alphabeta(p(n,1));           
    for j = 2:Nx-1                                      % interior points
    p(n+1,j) = p(n,j) + lam0*dt/dx*(v(j+1)*p(n,j+1)-v(j)*p(n,j)) + ...
        + lam0*dt*(v(j)*p(n,j))/(r+x(j)) + ...
       dt * (sign(p(n,j)-p_small)+1)/2 * alphabeta(p(n,j));
    end
   p(n+1,Nx) = p(n,Nx) - lam0*dt/dx*(v(end)*p(n,Nx))+ ... % end point
       + lam0*dt*(v(end)*p(n,Nx))/(r+x(end)) + ...
       dt * (sign(p(n,Nx)-p_small)+1)/2 * alphabeta(p(n,Nx));
   
  
   % growth of 'gas' portion of network
   f(n+1,:) = f(n,:) + dt*alphabeta(f(n,:));
   f(n+1,:) = f(n+1,:).*(sign(-p(n+1,:)+p_small)+1)/2;
            % there is no gas density if the density is connected
   
   
   % find first point where density is not connected
   qq = find(p(n+1,:)<p_small,1); 
   if min(p(n+1,:))>=p_small        % condition to avoid the above line returning null
       qq = Nx;
   end
   
   % growth of free boundary
   FBdelx = FBdelx + FBdx;
   if FBdelx>=dx
       % free boundary extends 1 point
       qq = min(qq+1,Nx); % free boundary can never go past end
       for j=1:qq
          p(n+1,j) = max(p(n+1,j),p_small+0.01);
       end
       FBdelx = 0.0;      % reset free boundary growth term
   end
   % update free boundary arrays
   FB(n+1) = dx*qq;
   vFB(n+1) = lam0*v(qq);
   
   % global connection reset
   if  f(n+1,min(qq+1,Nx))>p_small
       % if gas density at this point has exceeded connected density
       % threshold, then all of the network is connected (and the gas
       % density is reset in preparation for a new wave)
       p(n+1,:)=p(n+1,:)+f(n+1,:); % all of the network is connected
       f(n+1,:)=0;                 % zero out gas density
   end   
   
   
    % update the velocity
    FBindx = find(x>=FB(n+1),1); 
    if FB(n+1)>1 
        FBindx = length(x);
    end
    
    v = 0*x;

    % find the number of regions we have before FB
    j = find(p(n+1,1:FBindx)<p_cr,1); % find first point where density is not contractile
    if min(p(n+1,1:FBindx))>=p_cr        % condition to avoid the above line returning null
        % solve for just one u(r) between r0 and FB 
        vcoeffs = [r 1/r;  1 -1/FB(n+1)^2]\[0;-1];
        v(:) = vcoeffs(1)*(x(:)+r) +vcoeffs(2)./(x(:)+r);
    else
    rtrans = [r];
    rtrans = [rtrans x(j)+r]; 
    if j==length(x)
        % solve for just one u(r) between r0 and FB 
        vcoeffs = [r 1/r;  1 -1/FB(n+1)^2]\[0;-1];
        v(:) = vcoeffs(1)*(x(:)+r) +vcoeffs(2)./(x(:)+r);
    else
    contractileBOOL = 0;
    findingRegionsBOOL = 1;
    while findingRegionsBOOL == 1
        if contractileBOOL == 0
            j2 = j + find(p(n+1,j+1:FBindx)>=p_cr,1);
            if max(p(n+1,j+1:FBindx))<p_cr       % condition to avoid the above line returning null
                j2 = j + find(p(n+1,j+1:FBindx)<p_small,1);
                if min(p(n+1,j+1:FBindx))>=p_small
                    j2 = length(x);
                end
                rtrans = [rtrans x(j2)+r];
                findingRegionsBOOL = 0;
            else
                rtrans = [rtrans x(j2)+r];
                contractileBOOL = 1;
                j = j2;
            end
        else
            j2 = j + find(p(n+1,j+1:FBindx)<p_cr,1);
            if min(p(n+1,j+1:FBindx))>=p_cr        % condition to avoid the above line returning null
               rtrans = [rtrans x(FBindx)+r];
               findingRegionsBOOL = 0;
            else
                rtrans = [rtrans x(j2)+r];
                j = j2;
                contractileBOOL = 0;
            end
        end
        if j==length(x)
            findingRegionsBOOL = 0;
        end
    end

    nreg = length(rtrans)-1; % number of regions -- 2*nreg unknowns to solve for
    A = zeros(2*nreg,2*nreg);
    % B.C. at r0
    A(1,1:2) = [r 1/r];
    % B.C. at free boundary
    A(end,end-1:end) = [1 -1/rtrans(end)^2];
    % continuity of velocity eqns
    for j=1:nreg-1
        A(j+1,2*(j-1)+1:2*(j+1)) = [rtrans(j+1) 1/rtrans(j+1) -rtrans(j+1) -1/rtrans(j+1)];
        A(j+nreg,2*(j-1)+1:2*(j+1)) = [1 -1/rtrans(j+1)^2 -1 +1/rtrans(j+1)^2];
    end

    y = [zeros(nreg,1); (-1).^(transpose(1:1:nreg-1)); -((-1)^(nreg-1)+1)/2];
    % last term = 0 if nreg is even, = -1 if nreg is odd

    vcoeffs = A\y;

    for j = 1:nreg
        xindx = find(rtrans(j)-r<= x & rtrans(j+1)-r >=x);
        if j==nreg
            xindx = find(rtrans(j)-r<= x);
        end
        v(xindx) = vcoeffs(2*j-1)*(x(xindx)+r) +vcoeffs(2*j)./(x(xindx)+r);
    end

    end
    end
    v = -v;
    
   % Separate out velocity corresponding to connected portion
   vvv = v.*(sign(p(n+1,:)-p_small)+1)/2;
   
   
   v_array(n+1,:) = v;
   vvv_array(n+1,:) = vvv;
   
   if makeanimation==1
       if ((n/5)-floor(n/5))==0
             plot(x,p(n,:)+f(n,:),'b',x,lam0*vvv_array(n,:),'r',...
               x,p_small*ones(size(x)),'g',x,p_cr*ones(size(x)),'y'),
           hold on; xline(FB(n),'k--'); hold off
       axis([0 1 0 3]),
       txt = ['time: ' num2str(tn) ' units']; text(0.5,2,txt)
       txt = ['FreeBound: ' num2str(FB(n)) ' units']; text(0.5,1.7,txt)
       F(n)=getframe; writeVideo(ss,F(n));
       end
   end
   
   % move timestep forward
   tn = tnp;
end

if makeanimation==1
    close(ss);
end
