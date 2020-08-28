%% Forward control system: 
close all
clear all
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

% Boundary points: 
boundary_pts = circle(2,0,1);
% scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');
% [distance, closest_element, index] = distance_to_a_set([1;1],boundary_pts);

%%
L = linspace(0,2*pi,6);
xv = cos(L)';
yv = sin(L)';
scatter(xv,yv,12,'r','filled')



%%
[parameter,linear_data] = linear_regression([x,y]);
scatter(linear_data(:,1),linear_data(:,2),12,'r','filled');

%%
x0 = [0;0;0];
x_target = [4;0;0];
u = 0.01*ones(2*200,1);

x_min=-1;
x_max=5;
y_min=-2;
y_max=2;
theta_min=-pi/2;
theta_max=pi/2;

[traj_steer,u_steer] = steer(Phi,Psi_p,JPhi,p,x0,x_target,u,boundary_pts);

if norm(traj_steer(:,end)-x_target)<1e-3
    display('done!')
end



X_data = traj_steer(:,end);
U_data = u_steer;

while true     
for rep = 1:5
    
% 1/ Random goal
x_goal = [x_min+rand(1)*(x_max-x_min); y_min+rand(1)*(y_max-y_min); theta_min+rand(1)*(theta_max-theta_min)];

% 2/ Obtain U_best from X_data, U_data
[~, ~, closest_index] = distance_to_a_set(x_goal,X_data);
U_best = U_data(:,closest_index);

% 3/ Steer 
[traj_steer,u_steer] = steer(Phi,Psi_p,JPhi,p,x0,x_goal,U_best,boundary_pts);

if norm(traj_steer(:,end)-x_target)<1e-3
    display('done!')
    break
end

% 4/ Store Xn, u_steer
X_data(:,end+1) = traj_steer(:,end);
U_data(:,end+1) = u_steer;
end

% Introducing bias by having a random goal as x_target
% 2/ Obtain U_best from X_data, U_data
[~, ~, closest_index] = distance_to_a_set(x_target,X_data);
U_best = U_data(:,closest_index);
% 3/ Steer
[traj_steer,u_steer] = steer(Phi,Psi_p,JPhi,p,x0,x_target,U_best,boundary_pts);
if norm(traj_steer(:,end)-x_target)<1e-3
    display('done!')
    break
end

end






% 
function [traj_steer,u_steer] = steer(Phi,Psi_p,JPhi,p,x0,x_target,u,boundary_pts)
error_old = 1e6;
dt = 0.02;
T = 4;
iter_max = ceil(T/dt);
% x0 = [0;0;0];
% x_target = [4;0;0];
% rng('shuffle');
% u = 0.01*ones(2*iter_max,1); 
n_of_pts = 10;

figure(1)
% clf
hold on
plot_traj = plot(x0(1),x0(2),'ko','MarkerSize',8,'LineWidth',2.5);
scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');
grid on
% tic
while true
% tic    
    x = x0;
    x_traj = [];

    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(dt,x(1),x(2),x(3),u(2*iter-1),u(2*iter));        
        A_store{iter} = R_big(1:3,1:3);       % A in Ax+Bu 
        B_store{iter} = R_big(1:3,4:5);       % B in Ax+Bu

        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u(2*iter-1);u(2*iter)]); % x[k]-->x[k+1] might internally be broken into smaller steps to have better accuracy
        x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
        x = x(1:3);
        x_traj = [x_traj x];      
    end
    
    error = norm(x-x_target)
    improve = error_old - error;
    
    if error <= 1e-4 || improve < 1e-4
        break
    end
    error_old = error;
    
    % Calculate HH
    H = B_store{1};
    HH_store{1} = [H,zeros(3,2*(iter_max-1))];
    for iter = 2:iter_max
        H = [A_store{iter}*H, B_store{iter}];   
        HH_store{iter} = [H, zeros(3,2*(iter_max-iter))];
    end
    
    % Some plots
    figure(1)
    delete(plot_traj)
    plot_traj = plot(x_traj(1,:),x_traj(2,:),'k','LineWidth',2.5);
    start = plot(x0(1),x0(2),'ko','MarkerSize',8,'LineWidth',2.5);
    target = plot(x_target(1),x_target(2),'rx','MarkerSize',8,'LineWidth',1.5);
    axis([-5 5 -3 3]);
    drawnow
    
    
    % Establish constraint if needed 
    A = [];
    B = [];
    current_row = 1;
    for iter = 1:iter_max
        dist_to_boundary = vecnorm(x_traj(1:2,iter)-boundary_pts);
        [sorted_dist,indices] = sort(dist_to_boundary);
        
        if sorted_dist(1) < 0.1
            data = boundary_pts(:,indices(1:n_of_pts))'; % n-nearest boundary points 
            param = linear_regression(data);             % interpolating hyperplane            
            normal_vec = [param(1:end-1);-1];            % normal vector of the hyperplane
            
            % Positive if state and normal vector are on the same side of the hyperplane, else negative
            correct_sign = sign([normal_vec;param(end)]'*[x_traj(1:2,iter);1]);             
            A_k = [correct_sign*(normal_vec(1:2)/norm(normal_vec))', 0];
            
            A(current_row,:) = A_k*HH_store{iter};
            B(current_row,:) = 0.9*(0.05-sorted_dist(1));
            current_row = current_row + 1;
        end
        
    end
    
    
    if isempty(A)
        du = - (H'*H+1*error^2*eye(2*iter_max))\(H'*(x-x_target));
    else
        options = optimset('display','off');
        du = quadprog(H'*H + 0.1*error^2*eye(2*iter_max) , (x-x_target)'*H , -A, -B, [],[],[],[],[],options);
    end
    u = u + du;
% toc    
end 
% toc
u_steer = u;
traj_steer = x_traj;

end

function [parameter] = linear_regression(data) 
% data & linear_data: pxn matrix
% p: number of data points
% n: number of dimensions

[p, n] = size(data);
parameter = [data(:,1:end-1), ones(p,1)]\data(:,end);

% % Return a data using linear model if needed
% xn_linear = [data(:,1:end-1), ones(p,1)]*parameter;
% linear_data = [data(:,1:end-1),xn_linear];
end


function points = circle(x,y,r)
th = 0:pi/100:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
% hold on
% h = scatter(xunit, yunit,10,'b','filled');
% hold off
points = [xunit;yunit];
end

function [distance, closest_element, index] = distance_to_a_set(x,set)
% Arg = {x: a point; set: a set}
% Return = {distance to the set, closest element of the set, the closest element's index}
norm_ = vecnorm(x-set);
[distance,index] = min(norm_);
closest_element = set(:,index);
end