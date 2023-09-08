% authors: Alex Mogilner, Mariya Savinov

clear,
k=0.2; % contraction rate
pp=0.01; % prob of disassembly at each step
L=1;N=1*L; % size and number of new nodes assembled at each step
delta=0.02; del=0.04; % critical distance between nodes
n=N/pp; % expected average node density in steady state
dt=0.1; tend=80; T=tend/dt; tt=dt*(1:T);% vector of time
z=(0:0.05:1);zz=length(z); % spatial coordinates for computing density
dens=zeros(size(z)); dd=tt; Nclust=tt; % set up density vector
x=rand(1,1); x=sort(x); % initializing nodes
%v=VideoWriter('zz.avi');open(v); % setting up movie recording
for i=1:T %k=0.1+i/T;
    %k=0.2+0.15*0.5*(sign(i/T-0.3)+1)+0.65*0.5*(sign(i/T-0.5)+1);
    p=rand(1,length(x)); ss=find(p<pp); x(ss)=[];  
    % line above is disassembly; p is rand assigned to each node;
    % if p<pp, delete that node
    y=L*rand(1,N); %yy=-log(rand(1,2*N))/6;%yy=0.4*rand(1,2*N); 
    yy=[]; %yy=x(floor(length(x)*rand)+1)+del*(rand-0.5);
    x=[x y yy]; nn=length(x); %x=x+0.003*randn(1,nn);
    % adding N new nodes in random places and in a biased way
    x=sort(x); % sorting new array of nodes
    d=[x(1), diff(x)]; jj=find(d>del); % distances between nodes;
    % jj are distances greater than critical
    if length(jj)<1 x(1)=x(1)*(1-k*dt*0.5*(sign(delta-x(1))+1));
            for j=1:nn-1 x(j+1)=x(j)+d(j+1)*...
                    (1-k*dt*0.5*(sign(delta-d(j+1))+1)); zz=[]; end, else...
        % if ALL distances are short, the whole array contracts
kk=jj(1)-1; x(1)=x(1)*(1-k*dt*0.5*(sign(delta-x(1))+1));
    for j=1:kk-1 x(j+1)=x(j)+d(j+1)*...
            (1-k*dt*0.5*(sign(delta-d(j+1))+1)); end
        % else, first connected sub-array contracts to aggregate
for s=1:length(jj)-1 c=(x(jj(s))+x(jj(s+1)-1))/2;
    % find centers of all other connected sub-arrays
for j=jj(s):jj(s+1)-2 
    x(j+1)=x(j)+d(j+1)*(1-k*dt*0.5*(sign(delta-d(j+1))+1)); end
    cc=(x(jj(s))+x(jj(s+1)-1))/2; 
    x(jj(s):jj(s+1)-1)=x(jj(s):jj(s+1)-1)+c-cc; end
        % contract all connected sub-arrays to their centers
        c=(x(jj(end))+x(nn))/2;
for j=jj(end):nn-1  
    x(j+1)=x(j)+d(j+1)*(1-k*dt*0.5*(sign(delta-d(j+1))+1)); end
    cc=(x(jj(end))+x(nn))/2; 
    x(jj(end):nn)=x(jj(end):nn)+c-cc; end
        % contract the most distal (from aggregate) connected sub-array
plot(x,i*dt*ones(size(x)),'.k'), axis([0 L 0 tend]), hold on
        % plot all nodes at this time step
for b=1:21 dens(b)=length(find(x>0.05*b & x<0.05*(b+1))); end,
dd(i)=dens(16)+dens(17)+dens(18);
%xs=x(find(x<0.8));ds=[xs(1), diff(xs)]; jjj=find(ds>del);
%Nclust(i)=length(jjj); 
        % compute density by counting node # in each bin
%stairs(z,dens), axis([0 1 0 200]), F(i)=getframe; writeVideo(v,F(i));
        % plot density and capture the plot as movie frame
end, hold off
%NN=0.5*(1-sign(Nclust-0.5)); Y = fft(NN);
%P2 = abs(Y/T);P1 = P2(1:T/2+1); P1(2:end-1)=2*P1(2:end-1);
%fn = (1/dt)*(0:(T/2))/T;plot(fn,P1),axis([0 0.2 0 0.05])
%close(v); % make movie
%plot(tt,dd)