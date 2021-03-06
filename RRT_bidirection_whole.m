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


tic
[Phi_rev,Psi_p_rev,JPhi_rev] = compute_Phi_and_JPhi(p,-F,[xi1 xi2 xi3 xi4 xi5],dt);
toc


% Boundary points: 
% boundary_pts = circle(2,0,1);
% scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');
% [distance, closest_element, index] = distance_to_a_set([1;1],boundary_pts);

%% Another boundary points: 
r = 1.5;
left_angle=-pi/2:0.005:pi/2; 
left_xp=r*cos(left_angle);
left_yp=r*sin(left_angle);

x_hori = 0:0.01:4;
y_hori_up = 2*ones(1,length(x_hori));
y_hori_down = -2*ones(1,length(x_hori));

y_verti = -2:0.01:2;
x_verti = 4*ones(1,length(y_verti));

boundary_pts=[0+left_xp x_hori x_hori x_verti; 0+left_yp  y_hori_up y_hori_down y_verti];

pts1 = pts_between([0;-1.5],[0;-2],30);
pts2 = pts_between([0;1.5],[0;2],30);
boundary_pts = [boundary_pts,pts1(2:end-1,:)',pts2(2:end-1,:)'];

scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');

%% Draft zone
figure(1)
hold on
plot_traj = plot(x0(1),x0(2),'ko','MarkerSize',10,'LineWidth',2.5);
plot_target = plot(x_target(1),x_target(2),'rx','MarkerSize',10,'LineWidth',2.5);
scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');
grid on
axis([-3 8 -5 5]);
% [parameter,linear_data] = linear_regression([x,y]);
% scatter(linear_data(:,1),linear_data(:,2),12,'r','filled');
% [traj_steer,u_steer] = steer(Phi,Psi_p,JPhi,p,x_target,[3;-3;0],zeros(2*200,1),boundary_pts);
% [traj_steer,u_steer] = steer(Phi,Psi_p,JPhi,p,[4;-4;0],[6;-2;0],0.01*zeros(2*200,1),boundary_pts);
% [traj_steer,u_steer,reach] = steer(Phi_rev,Psi_p_rev,JPhi_rev,p,[3;-3;0],x_target,zeros(2*200,1),boundary_pts);


%% Bidirectional sampling + iterative steering (whole trajectory)
tic;

x0 = [0;0;0];
x_target = [5;0;0];

x_min=-3;
x_max=8;
y_min=-5;
y_max=5;
theta_min=-pi/2;
theta_max=pi/2;

figure(1)
clf
hold on
plot_traj = plot(x0(1),x0(2),'ko','MarkerSize',10,'LineWidth',2.5);
plot_target = plot(x_target(1),x_target(2),'rx','MarkerSize',10,'LineWidth',2.5);
scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');
grid on
axis([-3 8 -5 5]);


% Create X_data and U_data
X_data_start = [];
U_data_start = [];

X_data_target = [];
U_data_target = [];

while true     
% 1/ Randomize a goal, and extend x0-tree
x_goal = [x_min+rand(1)*(x_max-x_min); y_min+rand(1)*(y_max-y_min); theta_min+rand(1)*(theta_max-theta_min)];
plot_x_goal = plot(x_goal(1),x_goal(2),'kx','MarkerSize',7,'LineWidth',1);

[traj_steer,u_steer_1,X_data_start,U_data_start,reach] = steer(Phi,Psi_p,JPhi,p,x0,x_goal,X_data_start,U_data_start,boundary_pts);
delete(plot_x_goal)
plot(traj_steer(1,end),traj_steer(2,end),'ko','MarkerSize',5,'LineWidth',2);

% pause()

% 2/ Force random goal = Xn, and extend x_target-tree
x_goal = traj_steer(:,end);
plot_x_goal = plot(x_goal(1),x_goal(2),'rx','MarkerSize',7,'LineWidth',1);

[traj_steer,u_steer_2,X_data_target,U_data_target,reach] = steer(Phi_rev,Psi_p_rev,JPhi_rev,p,x_target,x_goal,X_data_target,U_data_target,boundary_pts);
delete(plot_x_goal)
plot(traj_steer(1,end),traj_steer(2,end),'r.','MarkerSize',10,'LineWidth',5);
if reach == 1
    part1_control = [u_steer_1;flip(kron(eye(200),[0 1; 1 0])*u_steer_2)];
    display('up-done!')
    break
end

% pause()


% 3/ Randomize a goal, and extend x_target-tree
x_goal = [x_min+rand(1)*(x_max-x_min); y_min+rand(1)*(y_max-y_min); theta_min+rand(1)*(theta_max-theta_min)];
plot_x_goal = plot(x_goal(1),x_goal(2),'rx','MarkerSize',7,'LineWidth',1);

[traj_steer,u_steer_3,X_data_target,U_data_target,reach] = steer(Phi_rev,Psi_p_rev,JPhi_rev,p,x_target,x_goal,X_data_target,U_data_target,boundary_pts);
delete(plot_x_goal)
plot(traj_steer(1,end),traj_steer(2,end),'r.','MarkerSize',10,'LineWidth',5);

% pause()


% 4/ Force random goal = Xn, and extend x0-tree
x_goal = traj_steer(:,end);
plot_x_goal = plot(x_goal(1),x_goal(2),'kx','MarkerSize',7,'LineWidth',1);

[traj_steer,u_steer_4,X_data_start,U_data_start,reach] = steer(Phi,Psi_p,JPhi,p,x0,x_goal,X_data_start,U_data_start,boundary_pts);
delete(plot_x_goal)
plot(traj_steer(1,end),traj_steer(2,end),'ko','MarkerSize',5,'LineWidth',2);

if reach == 1
    part1_control = [u_steer_4;flip(kron(eye(200),[0 1; 1 0])*u_steer_3)];
    display('down-done!')
    break
end


end

toc

%% part2: 
dt = 0.02;
T = 8;
iter_max = ceil(T/dt);

% u = part1_control;
n_of_pts = 10;
mu = 10;

% tic
for repeat=1:500   
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
    
    error = norm(x-x_target);
    
    [error, norm(u)]
    
    % Calculate HH
    H = B_store{1};
    HH_store{1} = [H,zeros(3,2*(iter_max-1))];
    for iter = 2:iter_max
        H = [A_store{iter}*H, B_store{iter}];   
        HH_store{iter} = [H, zeros(3,2*(iter_max-iter))];
    end
    
%     % Some plots
    figure(1)
    hold on
    delete(plot_traj)
    plot_traj = plot(x_traj(1,:),x_traj(2,:),'m','LineWidth',2.5);
%     start = plot(x0(1),x0(2),'ko','MarkerSize',8,'LineWidth',2.5);
%     target = plot(x_target(1),x_target(2),'rx','MarkerSize',8,'LineWidth',1.5);
%     axis([-3 8 -5 5]);
    drawnow
%     
    
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
    
    options = optimset('display','off');
    du = quadprog((1+mu)*eye(2*iter_max) , u , -A, -B, H, x_target - x,[],[],[],options);
        
    u = u + du;
    
end 
% toc


% 
function [traj_steer,u_steer,X_data,U_data,reach] = steer(Phi,Psi_p,JPhi,p,x0,x_target,X_data,U_data,boundary_pts)
clc
error_old = 1e6;
dt = 0.02;
T = 4;
iter_max = ceil(T/dt);
% x0 = [0;0;0];
% x_target = [4;0;0];
% rng('shuffle');
% u = 0.01*ones(2*iter_max,1); 
if isempty(X_data)
    u = zeros(2*iter_max,1);
else
    [~, ~, closest_index] = distance_to_a_set(x_target,X_data);
    u = U_data(:,closest_index);
end


n_of_pts = 10;

figure(1)
% % clf
% hold on
plot_traj = plot(x0(1),x0(2));
% scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');
% grid on
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
    
    if error <= 1e-3 
        reach = 1;
        break
    end
    if improve < 1e-3
        reach = 0;
%         x_traj(:,end)
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
    
%     % Some plots
    figure(1)
    hold on
    delete(plot_traj)
    plot_traj = plot(x_traj(1,:),x_traj(2,:),'m','LineWidth',2.5);
%     start = plot(x0(1),x0(2),'ko','MarkerSize',8,'LineWidth',2.5);
%     target = plot(x_target(1),x_target(2),'rx','MarkerSize',8,'LineWidth',1.5);
%     axis([-3 8 -5 5]);
%     drawnow
%     
    
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
        du = - (H'*H+1*error^1.3*eye(2*iter_max))\(H'*(x-x_target));
    else
%         size(A)
%         B = zeros(size(A,1),1);
        options = optimset('display','off');
        du = quadprog(H'*H + 5*error^1.3*eye(2*iter_max) , (x-x_target)'*H , -A, -B, [],[],[],[],[],options);
    end
    u = u + du;
% toc    
end 
% toc
delete(plot_traj)
u_steer = u;
traj_steer = x_traj;
X_data(:,end+1) = traj_steer(:,end);
U_data(:,end+1) = u_steer;

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

function pts = pts_between(x,y,n_inner_points)
xvals = linspace(x(1), y(1), n_inner_points+2);
yvals = linspace(x(2), y(2), n_inner_points+2);
pts = [xvals(:), yvals(:)];
end