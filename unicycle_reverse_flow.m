%% Forward control system: 
% clear all
clc
syms xi1 xi2 xi3 xi4 xi5 dt

F	 = [xi4*cos(xi3);
        xi4*sin(xi3);
        xi5;
        0;
        0];
p = 5;

tic
[Phi,Psi_p,JPhi] = compute_Phi_and_JPhi(p,F,[xi1 xi2 xi3 xi4 xi5],dt);
toc


%%
dt = 0.01;
T = 5;
iter_max = ceil(T/dt);
x0 = [2;2;0];
% x00 = [0.5;2.5;-pi/2];
x_target = [1;-1;0];
rng('shuffle');
u = 0.01*ones(2*iter_max,1); 

while true
    x = x0;
    x_traj = [];

    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(dt,x(1),x(2),x(3),u(2*iter-1),u(2*iter));        
        R_store{iter} = R_big(1:3,1:3);       % A in Ax+Bu 
        B_store{iter} = R_big(1:3,4:5);       % B in Ax+Bu

        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u(2*iter-1);u(2*iter)]); % x[k]-->x[k+1] might internally be broken into smaller steps to have better accuracy
        x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
        x = x(1:3);
        x_traj = [x_traj x];      
    end
    
    error = norm(x-x_target)
    if error <= 0.001
        break
    end
    % Calculate H
    H = B_store{1};
    for iter = 2:iter_max
        H = [R_store{iter}*H, B_store{iter}];
    end
    u = u - (H'*H+0.0001*eye(2*iter_max))\(H'*(x-x_target));
end 
u_forward = u;
traj_forward = x_traj;
figure
hold on
optimal_control = plot(x_traj(1,:),x_traj(2,:),'g:','LineWidth',2.5);
start = plot(x0(1),x0(2),'ko','MarkerSize',8,'LineWidth',2.5);
% new_start = plot(x00(1),x00(2),'ro','MarkerSize',8,'LineWidth',1.5);
target = plot(x_target(1),x_target(2),'rx','MarkerSize',8,'LineWidth',1.5);
% legend([start, target],{'start','target'})

%% Backward control system:
syms x1 x2 x3 x4 x5 dt

F_rev=-[x4*cos(x3);
        x4*sin(x3);
        x5;
        0;
        0];
p = 5;

tic
[Phi_rev,Psi_p_rev,JPhi_rev] = compute_Phi_and_JPhi(p,F_rev,[x1 x2 x3 x4 x5],dt);
toc

%%
dt = 0.01;
T = 5;
iter_max = ceil(T/dt);

u = flip(blkdiag(kron(eye(500),[0 1; 1 0]))*u_forward);

x0 = x_target;
x_target = [0.5;2.5;-pi/2];

x = x0;
x_traj = [];

for iter = 1:iter_max
    % The actual trajectory using adaptive step size mechanism
    [~, x_trajJ_fine] = adaptive_taylor(p,Phi_rev,Psi_p_rev,[0 dt],[x;u(2*iter-1);u(2*iter)]); % x[k]-->x[k+1] might internally be broken into smaller steps to have better accuracy
    x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
    x = x(1:3);
    x_traj = [x_traj x];      
end
traj_backward = x_traj;
%
while true
    x = x0;
    x_traj = [];

    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi_rev(dt,x(1),x(2),x(3),u(2*iter-1),u(2*iter));        
        R_store{iter} = R_big(1:3,1:3);       % A in Ax+Bu 
        B_store{iter} = R_big(1:3,4:5);       % B in Ax+Bu

        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi_rev,Psi_p_rev,[0 dt],[x;u(2*iter-1);u(2*iter)]); % x[k]-->x[k+1] might internally be broken into smaller steps to have better accuracy
        x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
        x = x(1:3);
        x_traj = [x_traj x];      
    end
    
    error = norm(x-x_target)
    if error <= 0.001
        break
    end
    % Calculate H
    H = B_store{1};
    for iter = 2:iter_max
        H = [R_store{iter}*H, B_store{iter}];
    end
    u = u - (H'*H+0.0001*eye(2*iter_max))\(H'*(x-x_target));
end 

hold on
new_start = plot(x_target(1),x_target(2),'ro','MarkerSize',8,'LineWidth',2.5);
optimal_control = plot(x_traj(1,:),x_traj(2,:),'r:','LineWidth',2.5);    


