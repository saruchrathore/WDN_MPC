% Project: 'Optimal control in water distribution'

% Given network parameters this function gives static and dynamic model equations as output.

% Incidence matrix of the network needs to be provided as an input
% It constructs reduced incidence matrix, chords, spanning tree, loop matrix
% of the network.
% These matrix are further used for creating model of the network


% Author: Saruch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Operation on the matrices

[n_H,m_H]=size(H);              % Calculate size of incidence matrix

% Check if input provided are correct

if m_H<(n_H-1)
    error('Number of edges cannot be less than (number of edge-1). Disconnected Graph')
elseif length(epsilon)~=m_H
    error('Number of elements in epsilon does not match with number of edge')
elseif length(diameter_pipe)~=m_H
    error('Number of elements in diameter_pipe does not match with number of edge')
elseif length(length_pipe)~=m_H
    error('Number of elements in length_pipe does not match with number of edge')
elseif length(height_node)~=n_H
    error('Number of elements in heigth_node does not match with number of nodes')
end

% Matrix operations. Finding reference node, chords, loop matrix

disp(['Total number of edges: ',num2str(m_H)]);
disp(['Total number of nodes: ',num2str(n_H)]);
disp(' ');

prompt = ['Please provide a reference node. Suggestion:', num2str(d_tau_index) ,'th node   '];
ref_node = input(prompt);   % Ask for the reference node
ref_node(isempty(ref_node))=4;    % If empty, select 4th node as reference node
non_ref_nodes=setdiff(1:n_H,ref_node);  % Non-reference nodes
disp(' ');



H_bar=H(non_ref_nodes,:);       % Generate reduced incidence matrix.(H_bar)

[n_H_bar,m_H_bar]=size(H_bar);              % Calculate size of reduced incidence matrix

n_chords= m_H_bar-n_H_bar;                  % Calculate number of chords

prompt = 'To select chords manually type "Y", otherwise press enter to automatically select chords';
selection = input(prompt,'s');   % Ask for whether to select chords manually
disp(' ');

if selection=='Y'
    prompt = ['Select ', num2str(n_chords), ' edges as chords'] ;
    set_chord=input(prompt);    % Select chords manually
    disp(' ');
    
    edge_tree=setdiff(1:m_H_bar,set_chord); % Remaining edges are edges of tree
    H_T_bar=H_bar(:,edge_tree);             % Generate reduced incidence matrix for tree
    det_H_T_bar=det(H_T_bar);
    
    if det_H_T_bar==0                       % Check if generated tree is valid. If yes continue
        error('Selected chords are not valid');  %Otherwise generate error message
    end
    
else
    chords = combnk(1:m_H_bar,n_chords);        % Possible combination for chords
    
    for i=1:size(chords,1)
        
        set_chord=chords(i,:);                  % Select a set of possible chords
        edge_tree=setdiff(1:m_H_bar,set_chord); % Remaining edges are edges of tree
        H_T_bar=H_bar(:,edge_tree);             % Generate reduced incidence matrix for tree
        
        det_H_T_bar=det(H_T_bar);
        
        if det_H_T_bar~=0                       % Check if generated tree is valid. If yes break
            break                               % break the loop.
        end
    end
end

edge_tree=sort(edge_tree);                  % Tree edges are in ascending order
set_chord=sort(set_chord);                  % Chord edges are in ascending order

disp (['Selected refernce node is  [',num2str(ref_node),']'])    % Display the selected reference node
disp (['Selected chords are  [',num2str(set_chord),']'])    % Display the selected chords


H_T=H(:,edge_tree);                         % Incidence matrix for the tree. H_T
H_T_bar=H_bar(:,edge_tree);                 % Reduce incidence matrix for the tree. H_T_bar
H_C=H(:,set_chord);                         % Incidence matrix for the chords. H_C
H_C_bar=H_bar(:,set_chord);                 % Reduce incidence matrix for the chords. H_C_bar

B_unsort=[eye(n_chords) -H_C_bar'*inv(H_T_bar')];   % Generate loop matrix with edge index as [chords tree].

idx=[set_chord edge_tree];                      % set index for unsorted B matrix

for i =1:size(B_unsort,1)
    B(i,idx)=B_unsort(i,:);
end


%% Finding M_p, M_c and M_tau matrix

% Finding M_p matrix, such that d=M_p*dp

% Start by finding M_p matrix such that d=M_p*dp. Then remove the row for the
% reference node to get M_p_bar

m_M_p=length(dp_index);             % Columns of M_p matrix
n_M_p=n_H;                         % Rows of M_p matrix

M_p=zeros(n_M_p,m_M_p);           % Preallocating for M_p matrix

for i=1:m_M_p
    M_p(dp_index(i),i)=1;     % Creating M_p matrix
end

M_p_bar=M_p(non_ref_nodes,:);   % Removing row corresponding to reference node

% Finding M_c matrix, such that d=M_c*dc

% Start by finding M_c matrix such that d=M_c*dc. Then remove the row for the
% reference node to get M_c_bar

m_M_c=length(dc_index);             % Columns of M_c matrix
n_M_c=n_H;                         % Rows of M_c matrix

M_c=zeros(n_M_c,m_M_c);           % Preallocating for M_c matrix

for i=1:m_M_c
    M_c(dc_index(i),i)=1;     % Creating M_c matrix
end

M_c_bar=M_c(non_ref_nodes,:);   % Removing row corresponding to reference node

% Finding M_tau matrix, such that d=M_tau*d_tau

% Start by finding M_tau matrix such that d=M_tau*d_tau. Then remove the row for the
% reference node to get M_tau_bar

m_M_tau=length(d_tau_index);             % Columns of M_tau matrix
n_M_tau=n_H;                         % Rows of M_tau matrix

M_tau=zeros(n_M_tau,m_M_tau);           % Preallocating for M_tau matrix

for i=1:m_M_tau
    M_tau(d_tau_index(i),i)=1;     % Creating M_tau matrix
end

M_tau_bar=M_tau(non_ref_nodes,:);   % Removing row corresponding to reference node


%% Finding pressure loss due to friction in elements(pipes in this case)

syms e_sym D_sym L_sym f_sym kf_sym       % Symbolic (epsilon, diameter, length, friction coeff, foam loss coeff.)

% Defining the friciton coeff. eq symbolically
friction_coeff_func=1.325*(log(e_sym/(3.7*D_sym)+5.74/R^0.9))^(-2);
friction_coeff_func=matlabFunction(friction_coeff_func);

% Calculating friction coeff. for each edge
friction_coeff=zeros(m_H,1);
for i=1:m_H
    friction_coeff(i,1)=friction_coeff_func(diameter_pipe(i),epsilon(i));
end

% Defining the friction loss eq symbolically
friction_loss_func=f_sym*8*L_sym*rho_fluid/(pi^2*D_sym^5*1e05*(36e02)^2);
friction_loss_func=matlabFunction(friction_loss_func);

% Calculating friction loss for each edge
friction_loss=zeros(m_H,1);
for i=1:m_H
    friction_loss(i,1)=friction_loss_func(diameter_pipe(i),length_pipe(i),friction_coeff(i));
end


% Calculating overall loss for each edge
lambda=2*friction_loss;          % Considering foam loss is equal to the friction loss

q_sym=sym('q_sym', [m_H 1]);             % Symbolic variable for flow in all the edges

% Calculating overall loss for each edge in [bar]
lambda_q_sym=lambda.*abs(q_sym).*q_sym;

%% Defining flow, demand and pressure as symbolic variables

d_sym=sym('d_sym',[n_H 1]);             % Symbolic variable for demand at each nodes

qc_sym=q_sym(set_chord);              % Symbolic variable for flow in the chords
qt_sym=q_sym(edge_tree);              % Symbolic variable for flow in the edges of tree
dp_sym=d_sym(dp_index);
dc_sym=d_sym(dc_index);
d_tau_sym=d_sym(d_tau_index);

%% Static equation for flows in the system

% Eq. q_T=-inv(H_T_bar)*H_C_bar*qc_sym+inv(H_T_bar)*F'*df
eq1_sym=qt_sym==-inv(H_T_bar)*H_C_bar*qc_sym+inv(H_T_bar)*M_p_bar*dp_sym+inv(H_T_bar)*M_c_bar*dc_sym+inv(H_T_bar)*M_tau_bar*d_tau_sym;

% Eq. B*lamda(q)=0
eq2_sym_qt_qc=B*lambda_q_sym==0;

rhs_eq1_sym=rhs(eq1_sym);

eq2_sym=subs(eq2_sym_qt_qc,qt_sym,rhs_eq1_sym);

%% Static equation for pressure at nodes

p_sym=sym('p_sym',[n_H 1]);     % Symbolic variable for pressure at each nodes

height_node_bar=height_node*rho_fluid*g/1e5;    % [bar] Pressure due to geodesic level at nodes
h_bar=height_node_bar(non_ref_nodes);           % [bar] Pressure due to geodesic level at non-reference nodes
h_0=height_node_bar(ref_node);                  % [bar] Pressure due to geodesic level at reference nodes
p_sym_bar=p_sym(non_ref_nodes);                 % Symbolic variable for pressure at non-reference nodes
p_0_sym=p_sym(ref_node);

% Eq. p_bar=inv(H_T_bar')*lambda(q_t)+(h_bar-1*h_0)+1*p_0
eq3_sym_qt=p_sym_bar==inv(H_T_bar')*lambda_q_sym(edge_tree)-(h_bar-ones(n_H-1,1)*h_0)+ones(n_H-1,1)*p_0_sym;

eq3_sym=subs(eq3_sym_qt,qt_sym,rhs_eq1_sym);

rhs_eq3_sym=rhs(eq3_sym);

%% Dynamic equations for the system

p_dot_sym=sym('p_dot_sym', [n_H 1]);          % Symbolic variable for pressure_dot at nodes

p_dot_tau_sym=p_dot_sym(d_tau_index);   % Symbolic variable for pressure_dot at nodes connected to tank

tau=rho_fluid*g*1./A_er/(1e5*36e2);            % Tank constant

% Eq. p_dot_tau=-T*d_tau
eq4_sym=p_dot_tau_sym==-tau*d_tau_sym;