% Project: 'Optimal control in water distribution'

% Given parameters of a network, the code first creates a model and then
% runs simulation test on the model based on demand curve and the input
% curve.
% The network is constructed of only pipes being edges, the pumps(or sources) and the
% consumers are modelled as independent flows in and out the network at
% different nodes.
% Tanks are also modelled as independent flows at node.

% Required parameters of the network for creating the model
% 1. Incidence Matrix
% 2. Nodes at which sources and consumers are connected
% 3. Nodes at which tanks are connected
% 4. Geodesic level of each nodes
% 5. Pipe parameter
%   5a. Length of pipe section
%   5b. Diameter of pipe section
%   5c. Roughness Height of pipe section
% 6. Reynolds number
% 7. Density of fluid
% 8. Gravitational acceleration


% Author: Saruch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear all;
close all;

%% Parameters of the network

% 1. Incidence Matrix

H=[ 1 0 1 0 0 0 0;         % Specify the incidence matrix
    -1 1 0 0 0 0 0;
    0 0 -1 1 1 0 0;
    0 -1 0 -1 0 1 0;
    0 0 0 0 -1 0 -1;
    0 0 0 0 0 -1 1];

% 2. Nodes at which sources and consumers are connected

dp_index=[1 6];         % [] Nodes connected to the pump
dc_index=[2 5];     % [] Nodes connected to the consumers

% d_index=sort([df_pump_index, df_consumer_index]);          % [] Nodes with non-zero demand


% 3. Nodes at which tanks are connected

d_tau_index=4;    % Nodes at which tanks are connected

% 4. Geodesic level of each nodes

height_node=[0; 0.9; 0; 3; 0.9; 0];  % [m] Geodesic level of each node

% 5. Pipe parameter

length_pipe=[10; 20; 20; 15; 10; 10; 25]; % [m] Length of pipes
diameter_pipe=[25; 20; 25; 15; 25; 20; 25]*1e-03; % [m] Diameter of pipes
epsilon=[0.05; 0.05; 0.05; 0.05; 0.05; 0.05; 0.05]*1e-03; % [m] Height of roughness inside the pipe

% 6. Reynolds number

R=4000;    % [] Reynolds number

% 7. Density of fluid

rho_fluid=997;  % [kg/m^3] Density of water


% 8. Gravitational acceleration

g=9.81;      % [m/s^2] Gravitational constant

%9. Tank parameters

A_er= 0.283;      % [m^2] Cross sectional area of the tank


%% Creating Model of the network
disp('Creating Model of the network...')

network_model

% Display the static and dynamic equations of the model
disp(' ');
fprintf(2, 'Static equations of the model are:\n')
disp(vpa(eq1_sym,3))
disp(vpa(eq2_sym,3))
disp(vpa(eq3_sym,3))
fprintf(2, 'Dynamic equations of the model are:\n')
disp(vpa(eq4_sym,3))
disp(' ');


%% Selection to run NMPC control simulation

prompt = 'Do you want to run nonlinear MPC control simulation? Press "Y" for yes, otherwise press enter to continue     ';
NMPC_order = input(prompt,'s');   % Ask for running open loop simulation
disp(' ');

if NMPC_order=='Y'
    
    %% NMPC Control
    
    controller_NMPC
    
end

%% Selection to run Kalman filter

prompt = 'Do you want to test Kalman filter predictor? Press "Y" for yes, otherwise press enter to continue     ';
Predictor_order = input(prompt,'s');   % Ask for running open loop simulation
disp(' ');

if Predictor_order=='Y'
    
    %% Kalman filter predictor
    
    predictor_kalman_filter

%% Selection to run NMPC with predicted consumer demand

prompt = 'Do you want to run NMPC with predicted consumer demand? Press "Y" for yes, otherwise press enter to continue     ';
NMPC_predictor_order = input(prompt,'s');   % Ask for running open loop simulation
disp(' ');

if NMPC_predictor_order=='Y'

%% NMPC with predictor
    
    controller_predictor

end

    
end

prompt = 'Do you want to run On-Off control simulation? Press "Y" for yes, otherwise press enter to continue     ';
open_loop = input(prompt,'s');   % Ask for running open loop simulation
disp(' ');

if open_loop=='Y'
    
    %% On-Off Control
    
on_off_control
    
end


disp('Thank you')
