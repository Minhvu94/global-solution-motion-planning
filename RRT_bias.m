clear all; close all; clc;

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

%% Boundary points: 
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

% scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');

%% Main RRT:
K=2000;
xMin=-3; xMax=8;
yMin=-5; yMax=5;
theta_Min=-pi/2;
theta_Max=pi/2;

% xInit=0; yInit=0; %initial point for planner
% xGoal=10; yGoal=10; %goal for planner
% thresh=1;           %acceptable threshold for solution

x0 = [0;0;0];
x_target = [5;0;0];

tree.vertex(1).x = x0(1);
tree.vertex(1).y = x0(2);
tree.vertex(1).theta = x0(3);
tree.vertex(1).xPrev = x0(1);
tree.vertex(1).yPrev = x0(2);
tree.vertex(1).thetaPrev = x0(3);
tree.vertex(1).u = zeros(100,1);


figure(1); hold on; grid on;
plot(x0(1), x0(2), 'ko', 'MarkerSize',10, 'MarkerFaceColor','k');
plot(x_target(1), x_target(2), 'go', 'MarkerSize',10, 'MarkerFaceColor','g');
scatter(boundary_pts(1,:), boundary_pts(2,:),10,'b','filled');
axis([-3 8 -5 5]);

for iter = 1:200
    % Randomize
    if mod(iter,10)==0
        X_Rand = x_target;
    else
        X_Rand = [xMin+rand*(xMax-xMin); yMin+rand*(yMax-yMin); theta_Min+rand*(theta_Max-theta_Min)];
    end
    
    % Find a nearest node 
    dist = zeros(1,length(tree.vertex));
    for j = 1:length(tree.vertex)
        dist(j) = sqrt( (X_Rand(1)-tree.vertex(j).x)^2 + (X_Rand(2)-tree.vertex(j).y)^2 + (X_Rand(3)-tree.vertex(j).theta)^2 );
    end    
    [val, ind] = min(dist);    
    X_near = [tree.vertex(ind).x; tree.vertex(ind).y; tree.vertex(ind).theta];
    
    % Steering
    [x_traj, u, reach] = steer(Phi,Psi_p,JPhi,p,X_near,X_Rand,boundary_pts);
           
    % Store data 
    tree.vertex(iter).x = x_traj(1,end); 
    tree.vertex(iter).y = x_traj(2,end);
    tree.vertex(iter).theta = x_traj(3,end);
    tree.vertex(iter).xPrev = X_near(1);
    tree.vertex(iter).yPrev = X_near(2);
    tree.vertex(iter).thetaPrev = X_near(3);
    tree.vertex(iter).u = u;
    tree.vertex(iter).traj = x_traj;

    scatter([tree.vertex(iter).x; tree.vertex(iter).xPrev],[tree.vertex(iter).y; tree.vertex(iter).yPrev], 12, 'r', 'filled');
    plot([X_near(1),x_traj(1,:)],[X_near(2),x_traj(2,:)],'b')
    drawnow
    if norm(x_traj(1:2,end)-x_target(1:2)) <= 1e-3
        break
    end
    
    pause(0);
    iter
end
CAC



%%
% function [xN, u] = steer(x0,x_target)
% xN = x_target;
% u = ones(100,1);
% end

% 
function [traj_steer,u_steer,reach] = steer(Phi,Psi_p,JPhi,p,x0,x_target,boundary_pts)
clc
error_old = 1e6;
dt = 0.01;
T = 0.2;
iter_max = ceil(T/dt);

u = zeros(2*iter_max,1);

n_of_pts = 10;

figure(1)
% % clf
% hold on
plot_traj = plot(x0(1),x0(2));
plot_target = plot(x_target(1),x_target(2),'rx','MarkerSize',8,'LineWidth',1.5);
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
    
    if isempty(A)
        du = - (H'*H+0.1*error^1.3*eye(2*iter_max))\(H'*(x-x_target));
    else
%         size(A)
%         B = zeros(size(A,1),1);
        options = optimset('display','off');
        du = quadprog(H'*H + 0.1*error^1.3*eye(2*iter_max) , (x-x_target)'*H , -A, -B, [],[],[],[],[],options);
    end
    u = u + du;
% toc    
end 
% toc
delete(plot_traj)
delete(plot_target)

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

function pts = pts_between(x,y,n_inner_points)
xvals = linspace(x(1), y(1), n_inner_points+2);
yvals = linspace(x(2), y(2), n_inner_points+2);
pts = [xvals(:), yvals(:)];
end