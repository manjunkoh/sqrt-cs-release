wall_y_pos = 1.5;
obstacle_radius = 1;

pos_idx = 1:2;

num_nodes = 50;
t_f = 30;
dt = t_f / num_nodes;
nx = 4;
nu = 2;

mu_0 = [0; 0; 0; 0];
mu_f = [10; 0; 0; 0];

Sigma_0 = diag([0.1, 0.1, 0.01, 0.01]);
Sigma_f = diag([0.04, 0.04, 0.01, 0.01]);

% Q = 0.001 * eye(nx);
Q = zeros(nx);
R = eye(nu);

A = [eye(2), dt*eye(2); zeros(2, 2), eye(2)];
B = [0.5*dt^2*eye(2); dt*eye(2)];
q = 0.005;
G = sqrt(q * dt) * eye(4);

A_sys = repmat(A, [1,1,num_nodes]);
B_sys = repmat(B, [1,1,num_nodes]);
G_sys = repmat(G, [1,1,num_nodes]);

control_risk = 0.005;
u_max = 0.15;
state_risk = 0.005;