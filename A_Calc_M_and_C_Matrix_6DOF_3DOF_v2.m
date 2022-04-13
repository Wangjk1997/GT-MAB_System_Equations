clear; clc; close all;
addpath('.\GNC')
syms rz_gb_b      real  % vertical position of CG in {b}
syms Igx Igy Igz  real  % inertial tensor at CG
syms m_RB         real  % rigid-body mass
syms m_Ax m_Ay m_Az I_Ax I_Ay I_Az real %added mass and inertia at CB
syms g real % 9.81 m/s/s

syms Dcbx Dcby Dcbz Dcbwx Dcbwy Dcbwz real%damping coefficients at CB

%body-fixed velocities of CG and CB in {n}, expressed in {b}
syms vx_bn_b vy_bn_b vz_bn_b real
syms vx_gn_b vy_gn_b vz_gn_b real
syms vx_bn_b_dot vy_bn_b_dot vz_bn_b_dot real
syms vx_gn_b_dot vy_gn_b_dot vz_gn_b_dot real
syms omega_x_bn_b omega_y_bn_b omega_z_bn_b  real    
syms omega_x_gn_b omega_y_gn_b omega_z_gn_b  real  
syms omega_x_bn_b_dot omega_y_bn_b_dot omega_z_bn_b_dot  real    
syms omega_x_gn_b_dot omega_y_gn_b_dot omega_z_gn_b_dot  real    

%Euler angles and their rate of change
syms phi phi_dot theta theta_dot psi psi_dot real

%body-fixed velocities of CG in {n}, expressed in {n}
syms vx_gn_n vy_gn_n vz_gn_n real


%% Equations at CB in {b}
r_g = [0 0 rz_gb_b]';        % location of the CG in {b}
r_b = [0;0;0];               % location of the CB in {b}
% nu_bn_b      = [vx_bn_b, vy_bn_b, vz_bn_b, omega_x_bn_b, omega_y_bn_b, omega_z_bn_b]'; % velocity vector, {b} in {n} seen in {b}

nu_gn_b      = [vx_gn_b, vy_gn_b, vz_gn_b, omega_x_gn_b, omega_y_gn_b, omega_z_gn_b]'; % velocity vector, {g} in {n} seen in {b}
nu_bn_b = inv(Hmtrx(r_g))*nu_gn_b % Eq 11.216, Fossen 2011
nu_bn_b_dot  = [vx_bn_b_dot, vy_bn_b_dot, vz_bn_b_dot, omega_x_bn_b_dot, omega_y_bn_b_dot, omega_z_bn_b_dot]'; 
nu_gn_b_dot  = [vx_gn_b_dot, vy_gn_b_dot, vz_gn_b_dot, omega_x_gn_b_dot, omega_y_gn_b_dot, omega_z_gn_b_dot]'; 

% 3DOF motion assumptions
% % Motion of CB in {n}, expressed in {b}
assume(phi==0); 
assume(psi==0); 
assume(vy_bn_b==0); 
assume(vy_bn_b_dot==0); 
assume(omega_x_bn_b==0);
assume(omega_x_bn_b_dot==0);
assume(omega_z_bn_b==0);
assume(omega_z_bn_b_dot==0);

% 3DOF motion assumptions
% % Motion of CG in {n}, expressed in {b}
assume(phi==0); 
assume(psi==0); 
assume(vy_gn_b==0); 
assume(vy_gn_b_dot==0); 
% assume(vz_gn_b==0); 
% assume(vz_gn_b_dot==0); 
assume(omega_x_gn_b==0);
assume(omega_x_gn_b_dot==0);
assume(omega_z_gn_b==0);
assume(omega_z_gn_b_dot==0);

% Rigid-body system inertia matrix at CG
MRB_CG=diag([m_RB,m_RB,m_RB,Igx,Igy,Igz]);
% Added mass and inertia at CB
MA_CB = diag([m_Ax, m_Ay, m_Az, I_Ax, I_Ay, I_Az]); 

% Total mass at CB
MRB_CB=Hmtrx(r_g)'*MRB_CG*Hmtrx(r_g);
M_CB=MA_CB+MRB_CB;

% Damping at CB
D_CB=diag([Dcbx Dcby Dcbz Dcbwx Dcbwy Dcbwz]);

% Coriolis and centripetal matrix at CB
CRB_CB = simplify(m2c(MRB_CB,nu_bn_b));
CA_CB  = simplify(m2c(MA_CB, nu_bn_b));
C_CB = simplify(CRB_CB+CA_CB)

% Restoration torques at CB
W=m_RB*g;
B=m_RB*g;
g_CB = gvect(W,B,theta,phi,r_g,r_b);

% Overall Equations at CB
simplify(M_CB*nu_bn_b_dot + C_CB*nu_bn_b + D_CB*nu_bn_b + g_CB)


%% Transfer to CG (refer to Sec 7.5.4 of Fossen 2011)
% Total mass at CG
MA_CG =inv(Hmtrx(r_g)')*MA_CB*inv(Hmtrx(r_g));
MRB_CG=MRB_CG;
M_CG=MA_CG+MRB_CG

% Coriolis and centripetal matrix at CG 
CA_CG =inv(Hmtrx(r_g)')*CA_CB*inv(Hmtrx(r_g));
CRB_CG=inv(Hmtrx(r_g)')*CRB_CB*inv(Hmtrx(r_g));
C_CG=CA_CG+CRB_CG

% Damping at CG
D_CG=inv(Hmtrx(r_g)')*D_CB*inv(Hmtrx(r_g))

% Restoration torques at CG
g_CG = inv(Hmtrx(r_g)')*g_CB

% Overall Equations at CG
simplify(M_CG*nu_gn_b_dot + C_CG*nu_gn_b + D_CG*nu_gn_b + g_CG)




