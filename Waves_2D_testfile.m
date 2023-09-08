% authors: Mariya Savinov, Alex Mogilner

clear all
close all

global makeanimation % if you wish to make an animation when calling the solver
makeanimation = 0;

%---------------------IndepVar/Parameters-------------------%
% domain extends from r to r+b
r = 0.1; b = 1;
% final time
tfinal = 40.0;

% spatial step
dx = 0.005;
% timestep1
dt = 0.0025;

% density thresholds
% gas <-> connected
p_small = 0.1;
% connected <-> contractile
p_cr = 0.6;

% density-dependent assembly parameters
% a(y) = w + gr*(1-exp(-y/ww))
gr = 0.9;
w = 0.1;
ww = 0.5;

% flux parameter to mimic free boundary
va = 0.154;

% velocity contractile strength
lam0 = 0.55;
    
%-----------------------------Vectors-----------------------%
% space vector
x = (0:dx:b);
Nx = (b)/dx+1; % number of x-points
% time vector
t = (0:dt:tfinal);

%-----------------------------runSim-----------------------%
% make directory, if it does not exist, to store simulations
if ~exist('simulations_2D','dir')
    mkdir('simulations_2D')
    sim_vars_list = {};
    save('simulations_2D/sim_vars.mat','sim_vars_list')
end

% assume the only identifying variabls are p_cr, p_small, lam0, and va

% check if there is already a directory for this parameter set
myFolder = strcat('simulations_2D/pcr_',num2str(p_cr),...
    '__psmall_',num2str(p_small),'__lam0_',num2str(lam0),'__va_',num2str(va));
if ~exist(myFolder, 'dir')
    mkdir(myFolder)
end

%saving variables
myvar = strcat('pcr_',num2str(p_cr),'__psmall_',num2str(p_small),...
            '__lam0_',num2str(lam0),'__va_',num2str(va));
    
filename = strcat(myFolder,'/variables_',myvar,'.mat');

if exist(filename,'file')
    usr_input = input('Do you wish to continue from last timepoint of results?\n 1 for yes, 2 for quit, else overwrites\n');
%     usr_input = 6;
    if usr_input==1
        
        filename = strcat(myFolder,'/sim_p_f_FB_v_vvv_',myvar,'.mat');
        load(filename,'p','f','FB','vFB','FBdelx','v','vvv');
        p0 = p(end,:);
        f0 = f(end,:);
        v0 = v(end,:);
        vvv0 = vvv(end,:);
        FB0 = FB(end);
        vFB0 = vFB(end);
        FBdelx0 = FBdelx;
        
        saveoverwrite = 0;
    elseif usr_input==2
        
        msg = "You chose to quit this script";
        h = msgbox(msg);
        error(msg)
    else
        
        p0 =  (p_cr + 0.1)*ones(Nx,1);
        v0 = x;           % initial velocity profile
        vvv0 = v0;        
        f0 = zeros(Nx,1); % gas density (=0)
        FB0 = b;
        vFB0 = lam0*max(v0);
        FBdelx0 = 0;

        saveoverwrite = 1;
    end
else
    load('simulations_2D/sim_vars.mat','sim_vars_list');
    sim_vars_list{end+1} = {p_cr,p_small,lam0,va};
    save('simulations_2D/sim_vars.mat','sim_vars_list');
    
    p0 =  (p_cr + 0.1)*ones(Nx,1);
    v0 = x;           % initial velocity profile
    vvv0 = v0;        
    f0 = zeros(Nx,1); % gas density (=0)
    FB0 = b;
    vFB0 = lam0*max(v0);
    FBdelx0 = 0;
        
    saveoverwrite = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve PDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% call PDE function
[p,f,FB,vFB,FBdelx,v,vvv] = Waves_2D_solver(p0,f0,v0,vvv0,FB0,vFB0,FBdelx0,r,x,dx,dt,tfinal,p_small,p_cr,gr,w,va,ww,lam0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% measure period by tracking changes in the free boundary location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% only interested in waves in last 40 time units
timeindx = find(t>=60,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% measure period by tracking changes in the free boundary location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% only interested in waves after 40 time units
timeindx = find(t>=60,1);

% if model breaks
modelNA_BOOL = 0;

if ~isempty(timeindx)

    % estimate the maximum distance traveled by FB
    disc_est = 1 - min(FB(timeindx:end));

    if disc_est < 5*dx
        wavesBOOL = 0;
        period_avg = 0;
        period_sd_avg = 0;
    else
        wavesBOOL = 1;
        disc_times = [];
        for j = 2:length(t)
            if FB(j)-FB(j-1)>disc_est/2
                disc_times = [disc_times t(j)];
                if p(j,find(x>=0.95,1))>p_cr
                    modelNA_BOOL = 1;
                end
            end
        end

        waveindx = find(disc_times>=60,1);
        disc_times_post60 = disc_times(waveindx:end);
        period_avg = mean(diff(disc_times_post60));
        period_sd_avg = std(diff(disc_times_post60));
        
        if isempty(disc_times)
            wavesBOOL = 0;
            period_avg = 0;
            period_sd_avg = 0;
        end
    end
else
    
    wavesBOOL = NaN;
    period_avg = NaN;
    period_sd_avg = NaN;

end

%% saving


if saveoverwrite==1
    %saving variables
    filename = strcat(myFolder,'/variables_',myvar,'.mat');
    save(filename,'r','b','dt','dx','gr','lam0','Nx','p_cr','p_small','w','ww','va')


    %save t,x,rho,v
    filename = strcat(myFolder,'/sim_p_f_FB_v_vvv_',myvar,'.mat');
    save(filename,'t','x','p','f','FB','vFB','FBdelx','v','vvv','wavesBOOL','modelNA_BOOL','period_avg','period_sd_avg')
    
else
    %saving variables
    filename = strcat(myFolder,'/variables_',myvar,'.mat');
    save(filename,'r','b','dt','dx','gr','lam0','Nx','p_cr','p_small','w','ww','va')
    
    % combine old and new data
    filename = strcat(myFolder,'/sim_p_f_FB_v_vvv_',myvar,'.mat');
    
    t_new = t;
    p_new = p;
    FB_new = FB;
    vFB_new = vFB;
    v_new = v;
    vvv_new = vvv;
    f_new = f;
    
    load(filename,'t','p','f','FB','vFB','v','vvv');
    
    t_old = t(1:end-1); % initial conditions of new sim are last timepoint of old sim
    p_old = p(1:end-1,:);
    FB_old = FB(1:end-1);
    vFB_old = vFB(1:end-1);
    v_old = v(1:end-1,:);   
    vvv_old = vvv(1:end-1,:);
    f_old = f(1:end-1,:);
    
    t = [t_old t(end)*ones(size(t_new))+t_new];
    p = [p_old; p_new];
    FB = [FB_old; FB_new];
    vFB = [vFB_old; vFB_new];
    v = [v_old; v_new];
    vvv = [vvv_old; vvv_new];
    f = [f_old; f_new];
    
    %save t,x,rho,v
    filename = strcat(myFolder,'/sim_p_f_FB_v_vvv_',myvar,'.mat');
    save(filename,'t','x','p','f','FB','vFB','FBdelx','v','vvv','wavesBOOL','modelNA_BOOL','period_avg','period_sd_avg')
    
end

%% plotting


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kymographs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_spaced = t(1):0.01:t(end);

[T,X] = meshgrid(x,t);
[Tq,Xq] = meshgrid(x,t_spaced);

f_spaced = interp2(T,X,f,Tq,Xq);
p_spaced = interp2(T,X,p,Tq,Xq);
vvv_spaced = interp2(T,X,vvv,Tq,Xq);

% kymographs of  density and velocity
fig_v_and_colorplot = figure(1);
clf
fig_v_and_colorplot.Position = [100, 100, 1200, 1200];

ax1 = subplot(2,1,1);
imagesc(x,t_spaced,f_spaced+p_spaced);
set(gca, 'YDir','reverse')
ylabel('time $t$','Interpreter','latex','FontSize',17)
xlabel('position $x$','Interpreter','latex','FontSize',16)
colorbar
set(ax1,'XMinorTick','on','YMinorTick','on')
title('Gas density $f(x,t)$','interpreter','latex','FontSize',18)

ax2 = subplot(2,1,2);
imagesc(x,t_spaced,lam0*vvv_spaced);
set(gca, 'YDir','reverse')
ylabel('time $t$','Interpreter','latex','FontSize',17)
xlabel('position $x$','Interpreter','latex','FontSize',16)
colorbar
set(ax1,'XMinorTick','on','YMinorTick','on')
title('Velocity $v(x,t)$','interpreter','latex','FontSize',18)

figfilename = strcat(myFolder,'/kymographs.png');
saveas(fig_v_and_colorplot,figfilename)
