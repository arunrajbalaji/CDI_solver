clear;clc;
%V16
for initial_condition_solver = [1 0]
%% Input Numerics
% Mesh Points
N1 = 10;
N2 = 10;

N = N1+N2;
optimal_mesh = 1;
uniform = 1;
steady_state_solver = 0;
driving_force = 1;% 1 - Potential Driven; 2 - Current Driven
current_collect_index = floor(N1/2);
t_restart = 0;

if initial_condition_solver == 1
    initial_guess = 1;
    steady_state_solver = 1;
    save_all_name = 'CV_initial_ALL.mat';
    save_name = 'CV_initial.mat';
else
    initial_guess = 0;
    load_name = 'CV_initial.mat';
    save_all_name = 'CV_ALL.mat';
    save_name = 'CV.mat';
end

% Species
s_mo = 2;       % Mobile species Num
eqn_ma = s_mo+1;   % Macropore Equation number for each mesh unit
s_ch = 0;       % Immobile chemical species Num
eqn_mi = eqn_ma+s_ch;   % Macropore Equation number for each mesh unit
iteration_num = 6;
Bulk_Reaction_switch = 0;
Faradaic_switch = 0;

% Time & Plot
dt_initial = 1e-1;
dt = dt_initial;
if initial_condition_solver == 1
    dt_interval = 1e1;
else
    dt_interval = 1e-1;
end
plot_dt_check = @(t) mod(t,10);

dt_count = 0;
if steady_state_solver == 1
    t_final = 1e1;        % final time
    dt_max = dt_interval/2;
    dt_min = 0.001;
    t_plan = linspace(0,t_final,ceil(t_final/dt_interval)+1);   % desired plot times
else
    t_final = 1e4;        % final time
    dt_max = dt_interval;
    dt_min = 0.001;
%     t_plan = linspace(0,t_final,ceil(t_final/dt_interval)+1);   % desired plot times
    t_plan = linspace(t_restart,t_final,ceil((t_final-t_restart)/dt_interval)+1);   % desired plot times
end

pn = length(t_plan);          % number of desired plots
pc = 1;               % plot counter
dt_stored_switch = 0;

% Saving Data
save_path_format = 'CV_curve//';
save_path = sprintf(save_path_format);
mkdir(save_path)
t_saved = nan(t_final*2,1);
t_plotting_array = nan(t_final*2,1);
c_eff_array = nan(t_final*2,1);
current_array = nan(t_final*2,1);
V_plotting_array = nan(t_final*2,1);
ratio_phi_D_array = nan(N,t_final*2);
dt_plotting_array = nan(t_final*2,1);
t_plotting_index = 1;

v_saved = nan(N*eqn_ma,t_final*2);
v2_saved = nan(N*eqn_mi,t_final*2);
phi_el_saved = nan(t_final*2,1);
saved_data_index = 1;

%% Initialization
% Input Physics
if initial_condition_solver == 1
    l_e = 1e-4/1e2; %m
else
    l_e = 6e-4; %m
end
l_sep = 0;%65e-6+5/2*l_e;%65e-6; %m
l_res = 1*l_e;
L = 1*l_res+1*l_e+l_sep;

R = 8.3145;    %J/(mol·K)
T = 298;    %K
F = 96485;   %C/mol
V_T = R*T/F;    %V

A_c = 25*1e-6;  %m^2
Q_vol = 0;%1e-6/60;    %m^3/s
u = Q_vol/A_c;

C_S = 200e6;  %F/m^3, Stern Capacitance,200e6,200e-7
C_external = 29.4;% F/m^2, External Capacitance
m_e = 0.109*1e-3; %kg
EER = 0; %Ohm, 2.57

% Operation Parameters
scan_rate = 10e-3;%V/s
V_magnitude = 0.6;
FCT = V_magnitude/scan_rate*4;%s
Switch_Time = FCT/2;
Smoothing_Zone = 0.1;%s

if driving_force == 1
    V_discharge = 0;%V 
else
    if initial_condition_solver == 1
        I_app = 0;% A/m^2 
    else
        I_app = 0.3;% A/m^2
    end
end
p_ma_value = 0.7;
p_mi_value = 0.172;%0.172;
p_res_value = 1;

% Mobile group
index_R = 10;
index_Cl = 2;
index_O = 30;
index_K = 1;
index_H = 4;
index_OH = 5;

% Immobile group
index_COOH = eqn_ma+1;
index_COO = eqn_ma+2;

z_R = 1;
z_Cl = -1;
z_O = 2;
z_K = 1;
z_H = 1;
z_OH = -1;

z_COOH = 0;
z_COO = -1;
z_mo_input = [z_K;z_Cl;z_O;z_K;z_H;z_OH];
% z_mo_input = [z_R;z_Cl;z_O];
z_chem = [z_COOH;z_COO];

D_R = 2.57e-9;  %m^2/s
D_Cl = 2.03e-9;%2.03e-9;  %m^2/s
D_H = 9.13e-9;
D_OH = 5.16e-9;
D_O = 2.57e-9;
D_K = 1.8e-9;
D_input = [D_K;D_Cl;D_O;D_K;D_H;D_OH];

mu_R = 0;
mu_Cl = 0;
mu_O = 0;
mu_H = 0;
mu_OH = 0;
mu_K = 0;
mu_att = [mu_R;mu_Cl;mu_O;mu_K;mu_H;mu_OH];

pH_feed = 7;

c_f_Cl = 1000;    % mM or mol/m^3
c_f_H = 10^(3-pH_feed);
c_f_OH = 10^(-11+pH_feed);
c_f_O = 10;
c_f_R = 0;
c_f_K = 1000;
fixed_chem_charge = -0;

c_i_HX = 0;
c_i_X = 0;
c_i_YH = c_i_X;
c_i_Y = c_i_HX;
c_i_COOH = 100;
c_i_COO = 0;
c_f = [c_f_K;c_f_Cl;c_f_O;c_f_K;c_f_H;c_f_OH];  %mol/m^3, feed concentration
c_i = [c_f_K;c_f_Cl;c_f_O;c_f_K;c_f_H;c_f_OH;0];
c_i_ch = [c_i_COOH;c_i_COO];
c_ref = 1;  %mol/m^3

k_ma = [1e0,1e8];
R_bulk = [0, index_H;
          0, index_OH];
P_bulk = [index_H,  0;
          index_OH, 0];

k_mi = [1e0,1e8];

R_bulk_mi = [0, index_H;
             0, index_OH];
P_bulk_mi = [index_H, 0;
             index_OH, 0];

k_chem = [1e4,1e8]/1e2;

R_chem = [index_COOH, index_H;
          0,          index_COO];
P_chem = [index_H,    index_COOH;
          index_COO,  0];


L_reservoir_domain_check = @(x) x<0;
anode_domain_check = @(x) 0;
seperator_domain_check = @(x) 0;
cathode_domain_check = @(x) (0<=x).*(x<=l_e);
R_reservoir_domain_check = @(x) 0;

% Initialization - Mesh
% x_f = linspace(-l_res,l_e,N+1)';
if optimal_mesh > 0
    x_f = optimal_mesh_gen(1,N1,-l_res,0,3);
    mesh = optimal_mesh_gen(1,N2,0,l_e,1);
    x_f = [x_f; mesh(2:end)];
else
    if uniform > 0
        x_f = linspace(-l_res,0,N1+1)';
        x_f = [x_f; linspace(l_e/N2,l_e,N2)'];
        x_f = [x_f; linspace(l_e+l_sep/N3,l_e+l_sep,N3)'];
        x_f = [x_f; linspace(l_e+l_sep+l_e/N4,2*l_e+l_sep,N4)'];
        x_f = [x_f; linspace(2*l_e+l_sep+l_res/N5,L-l_res,N5)'];
    else
        x_f(1,1) = -l_res;
        x_f = [x_f; nonuniform_generator(l_res,0,N1,x_f)];
        x_f = [x_f; nonuniform_generator(l_e,l_e,N2,x_f)];
        x_f = [x_f; nonuniform_generator(l_sep,l_e+l_sep,N3,x_f)];
        x_f = [x_f; nonuniform_generator(l_e,2*l_e+l_sep,N4,x_f)];
        x_f = [x_f; nonuniform_generator(l_res,2*l_e+l_sep+l_res,N5,x_f)];
        N = length(x_f)-1;
    end
end
% Strecth the first mesh to make the first mesh center overlap with the physical one
x_f(1) = x_f(1)-(x_f(2)-x_f(1));

x_c = 0.5*(x_f(1:end-1)+x_f(2:end));

delta_x_c = (x_f(2:end)-x_f(1:end-1));

delta_x_f = (x_c(2:end)-x_c(1:end-1));



figure(10)
plot(x_c,1,'-o');
xline([-l_res 0 l_e])

% Initialization - Memory
v_plot = nan(N*eqn_ma,length(t_saved));
v = nan(N*eqn_ma,1);
v2 = nan(N*eqn_mi,1);

% Initialization - Physics Property Vector
z_mo = nan(N*eqn_ma,1);
for i = 1:s_mo
    z_mo(i:eqn_ma:end) = z_mo_input(i);
end

p = nan(N*eqn_ma,1);
for i = 1:s_mo+1
    p(i:eqn_ma:end) = L_reservoir_domain_check(x_c)*p_res_value+cathode_domain_check(x_c).*p_ma_value;
end

tau = p.^(-1/2);

p_mi = nan(N*eqn_ma,1);
for i = 1:s_mo+1
    p_mi(i:eqn_ma:end) = L_reservoir_domain_check(x_c)*0+cathode_domain_check(x_c).*p_mi_value;
end

D = nan(N*eqn_ma,1);
for i = 1:s_mo
    D(i:eqn_ma:end) = D_input(i);
end
D_eff = D.*p./tau;
D_eff_f = (D_eff(1:(N-1)*eqn_ma)+ D_eff(1+eqn_ma:eqn_ma*N))/2;


%% Faradaic Initialization

E_A = 0;
E_C = 0;

Fara_Rea_num_A = 1;

index1_J2_Fara_A = nan((N*eqn_mi)^2,1);
index2_J2_Fara_A = index1_J2_Fara_A;
element_J2_Fara_A = index1_J2_Fara_A;

index1_J_Fara_A = nan(Fara_Rea_num_A*(N*eqn_ma)^2,1);
index2_J_Fara_A = index1_J_Fara_A;
element_J_Fara_A = index1_J_Fara_A;

%%Cathode
Fara_Rea_num_C = 1;
Oxi_C = index_O;
Red_C = index_R;
s_coeff_Oxi_C = 1;
s_coeff_Red_C = 1;

n_ele_C = [1];

J_C = 1e-9; %A/m^2
K_R = 1e-12;
K_O = K_R/exp(-0.45/V_T);
beta_C = sqrt(K_R*K_O);
OCV_rel_dimensional = log(K_R/K_O)*V_T;%E_C - E_A;

index1_J2_Fara_C = nan((N*eqn_mi)^2,1);
index2_J2_Fara_C = index1_J2_Fara_C;
element_J2_Fara_C = index1_J2_Fara_C;

index1_J_Fara_C = nan(Fara_Rea_num_C*(N*eqn_ma)^2,1);
index2_J_Fara_C = index1_J_Fara_C;
element_J_Fara_C = index1_J_Fara_C;


%% Initial Condition
if initial_guess == 1
    for i = 1:s_mo
        v(i:eqn_ma:end) = c_i(i);   %c_s
    end
    v(s_mo+1:eqn_ma:end) = linspace(0,0,N);   %phi
    
    for i = 1:s_mo
        v2(i:eqn_mi:end) = L_reservoir_domain_check(x_c)*0+cathode_domain_check(x_c)*(c_i(i)*exp(mu_att(i)));   %cs_mi
    end
    v2(s_mo+1:eqn_mi:end) = OCV_rel_dimensional/F*C_S;  %Sigma_ele
    if s_ch ~= 0
        for i = 1:s_ch
            v2(eqn_ma+i:eqn_mi:end) = L_reservoir_domain_check(x_c)*0+cathode_domain_check(x_c)*c_i_ch(i);   %cs_mi
        end
    end
    t = 0;

    % Calculate Initial Current
    current_mag = current_calculator(current_collect_index,eqn_ma,s_mo,z_mo_input,D_eff_f,v,delta_x_f,F,A_c);

    if driving_force == 1
        phi_el_value = 0;
    else
        phi_el_value = -OCV_rel_dimensional/V_T;
    end
    phi_el = (-phi_el_value).*cathode_domain_check(x_c);
else
    %%Inherit intial condition
    load([save_path,load_name])
    t = 0;
    % Calculate Initial Current
    current_mag = current_calculator(current_collect_index,eqn_ma,s_mo,z_mo_input,D_eff_f,v,delta_x_f,F,A_c);

    phi_el = (-phi_el_value).*cathode_domain_check(x_c);
end

%% Memory Allocation
index1_J = nan((N*eqn_ma)^2,1);
index2_J = index1_J;
element_J = index1_J;

b = nan(N*eqn_ma,1);

index1_J2 = nan((N*eqn_mi)^2,1);
index2_J2 = index1_J2;
element_J2 = index1_J2;

index1_J3 = nan((N*eqn_mi)^2,1);
index2_J3 = index1_J3;
element_J3 = index1_J3;

b2 = zeros(N*eqn_mi,1); % Zero helps to prevent extra nan in b2

index1_J2_chem = nan((N*eqn_mi)^2,1);
index2_J2_chem = index1_J2_chem;
element_J2_chem = index1_J2_chem;

index1_J2_chem_Rea = index1_J2_chem;
index2_J2_chem_Rea = index1_J2_chem;
element_J2_chem_Rea = index1_J2_chem;

index1_J_Bulk_Rea = nan(2*(length(R_bulk)^2)*(N*eqn_ma)^2,1);
index2_J_Bulk_Rea = index1_J_Bulk_Rea;
element_J_Bulk_Rea = index1_J_Bulk_Rea;

index1_J_Bulk_Rea_mi = nan(2*(length(R_bulk_mi)^2)*(N*eqn_ma)^2,1);
index2_J_Bulk_Rea_mi = index1_J_Bulk_Rea_mi;
element_J_Bulk_Rea_mi = index1_J_Bulk_Rea_mi;

index1_J_chem_Rea = nan(2*(length(R_chem)^2)*(N*eqn_ma)^2,1);
index2_J_chem_Rea = index1_J_chem_Rea;
element_J_chem_Rea = index1_J_chem_Rea;



%% Processor
if t_plan(saved_data_index) == t
    t_saved(saved_data_index) = t;   % for plot y vs time
    v_saved(:,saved_data_index) = v;
    v2_saved(:,saved_data_index) = v2;
    phi_el_saved(saved_data_index) = -phi_el_value;
    saved_data_index = saved_data_index+1; % for plot y vs time
end
if initial_condition_solver ~= 1
    save([save_path,save_all_name]);
    save_function(0,initial_condition_solver,v,v2,t,phi_el_value,save_name,save_all_name,save_path);
end
while 1
    
    % Applied voltage to the whole circuit
    Charging_smoothing_ss = 1 - exp(-(t+dt)/Smoothing_Zone);
    Charging_smoothing = 1 - exp(-mod((t+dt),FCT)/Smoothing_Zone);
    Discharging_smoothing = 1 - exp(-(mod(t+dt-Switch_Time,FCT))/Smoothing_Zone);

    if steady_state_solver == 1
        V_app = V_discharge*Charging_smoothing_ss;
    else
        V_app = scan_rate*t.*(0<=t).*(t<=FCT/4)+...
                 (V_magnitude-scan_rate*mod(t-FCT/4,FCT)).*(0<mod(t-FCT/4,FCT)).*(mod(t-FCT/4,FCT)<Switch_Time)+...
                 (-V_magnitude+scan_rate*mod(t-FCT*3/4,FCT)).*(0<=mod(t-FCT*3/4,FCT)).*(mod(t-FCT*3/4,FCT)<=Switch_Time).*(t>FCT/4);
        V_app = -V_app;
    end
    
    %%Plot while Solving
    t_plotting_array(t_plotting_index) = t;

    % Effluent Concentration Calculator
%     c_index = @(N) (N-1)*eqn_ma;
%     c_eff_array(t_plotting_index) = delta_x_c(N-N5)*(v(1+c_index(N-N5))-v(1+c_index(N-N5-1)))/delta_x_f(N-N5-1) + (v(1+c_index(N-N5))+v(1+c_index(N-N5-1)))/2;
    % Current Magnitude Saving    
    current_array(t_plotting_index) = current_mag;
    V_plotting_array(t_plotting_index) = V_app;
%     ratio_phi_D_array(:,t_plotting_index) = (phi_el-v(eqn_ma:eqn_ma:end)-v2(eqn_ma:eqn_mi:end)*F/V_T/C_S)/(-phi_el_value);
    dt_plotting_array(t_plotting_index) = dt;
    t_plotting_index = t_plotting_index+1;


    if mod(t_plotting_index,100) == 0

%         figure(1)
%         clf
%         hold on
%         for i = [1 2 3]%[1 3]
%             plot(x_c,v(i:eqn_ma:N*eqn_ma),'-o')
%         end
%         hold off
%         text(0,1,sprintf('t=%f',t))
% %         xlim([0.75*l_e l_e+l_sep+0.25*l_e])
%         xlabel('x')
%         ylabel('c')
%         legend('R','O')
%         drawnow
% %             
%         figure(2)
%         clf
%         hold on
%         for i = [1 2 3]%[1 3]
%             plot(x_c,v2(i:eqn_mi:N*eqn_mi),'-o')
%         end    
% %         for i = eqn_ma+1:eqn_mi
% %             plot(x_c,v2(i:eqn_mi:N*eqn_mi),'-o')
% %         end    
% %         plot(x_c,v2(index_HX:eqn_mi:N*eqn_mi)+v2(index_X:eqn_mi:N*eqn_mi),'-o')
% %         plot(x_c,v2(index_HY:eqn_mi:N*eqn_mi)+v2(index_Y:eqn_mi:N*eqn_mi),'-o')
%         hold off
%         xlabel('x')
%         ylabel('c_{mi}')
%         legend('R','Cl','O')
%         drawnow
% %     
%         figure(3)
%         clf
%         hold on
%         plot(x_c,v(eqn_ma:eqn_ma:N*eqn_ma),'-o')
%         plot(x_c,v2(eqn_ma:eqn_mi:N*eqn_mi)*F/V_T/C_S,'-o')
%         plot(x_c,phi_el-v(eqn_ma:eqn_ma:N*eqn_ma)-v2(eqn_ma:eqn_mi:N*eqn_mi)*F/V_T/C_S,'-o')
%         plot(x_c,phi_el,'-o')
%         hold off
%         xlabel('x')
%         ylabel('\phi')
%         legend('\phi','\phi_{St}','\phi_{Donnan}','\phi_{electrode}')
%         drawnow
% %     
%         figure(4)
%         clf
%         hold on
%         plot(x_c,v2(s_mo+1:eqn_mi:N*eqn_mi),'-o')
%         hold off
%         xlabel('x')
%         ylabel('\sigma_{ele}')
%         drawnow
    % 
%         figure(5)
%         clf
%         hold on
%         plot(x_f(2:end-1)',flux_1_E,'-o')
%         plot(x_f(2:end-1)',flux_1_D,'-o')
%         hold off
%         legend('ELE','Diff')
%         xlabel('x')
%         ylabel('Flux')
%         drawnow
    %     
%         figure(6)
%         clf
%         hold on
%         for i = index_H
%             plot(x_c,-log10(v(i:eqn_ma:N*eqn_ma)/1e3),'-o')
%         end
%         plot(x_c,-log10(v(index_H:eqn_ma:N*eqn_ma)/1e3) -log10(v(index_OH:eqn_ma:N*eqn_ma)/1e3),'-o')
%         hold off
%         text(3e-3,14,sprintf('t=%.2f',t))
%     %     xlim([0 L-2*l_res])
%     %     ylim([12 16])
%         xlabel('x')
%         ylabel('pH_{mA},(pH+pOH)_{mA}')
%         drawnow
%     
%         figure(7)
%         clf
%         hold on
%         for i = index_H
%             plot(x_c,-log10(v2(i:eqn_mi:N*eqn_mi)/1e3),'-o')
%         end
%         plot(x_c,-log10(v2(index_H:eqn_mi:N*eqn_mi)/1e3) -log10(v2(index_OH:eqn_mi:N*eqn_mi)/1e3),'-o')
% %         plot(x_c,-log10(v2(index_H:eqn_mi:N*eqn_mi).*v2(index_X:eqn_mi:N*eqn_mi)./v2(index_HX:eqn_mi:N*eqn_mi)/1e3),'-o')
% %         plot(x_c,-log10(v2(index_H:eqn_mi:N*eqn_mi).*v2(index_Y:eqn_mi:N*eqn_mi)./v2(index_HY:eqn_mi:N*eqn_mi)/1e3),'-o')
%         plot(x_c,-log10(v2(index_H:eqn_mi:N*eqn_mi).*v2(index_COO:eqn_mi:N*eqn_mi)./v2(index_COOH:eqn_mi:N*eqn_mi)/1e3),'-o')
% 
%         hold off
%         text(3e-3,14,sprintf('t=%.2f',t))
%         legend('pH','pH+pOH','pKa5','pKa9','pKa9.5')
%     %     ylim([12 16])
%         xlabel('x')
%         ylabel('pH & pKa')
%         drawnow

%         figure(8)
%         clf
%         hold on
%         plot(t_plotting_array,c_eff_array,'g')
%         xlabel('t')
%         ylabel('c_eff')
% 
%         figure(9)
%         clf
%         hold on
% 
%         plot(x_c,v2(index_COO:eqn_mi:end),'-o')
%         plot(x_c,v2(index_COO:eqn_mi:end)+v2(index_COOH:eqn_mi:end),'-o') % Total
% 
%         hold off
%         text(3e-3,14,sprintf('t=%.2f',t))
%         xlim([0 2*l_e+l_sep])
%         legend('COO','COOH+COO')
% %         ylim([12 16])
%         xlabel('x')
%         ylabel('Concentration[mM]')
%         drawnow

        figure(11)
        clf
        hold on
%         c_index = @(N) (N-1)*eqn_ma;
        plot(t_plotting_array,current_array,'g')
        xlabel('t')
        ylabel('current')
%         ylim([0.18 0.82])

        figure(12)
        clf
        hold on
%         c_index = @(N) (N-1)*eqn_ma;
        plot(V_plotting_array,current_array,'o')
        xlabel('V')
        ylabel('current')

        figure(13)
        clf
        hold on
%         c_index = @(N) (N-1)*eqn_ma;
        plot(t_plotting_array,dt_plotting_array,'g')
        xlabel('t')
        ylabel('dt')
    end

    if ( t == t_saved(pc) )
        v_plot(:,pc) = v;  % v_plot is the saved variable
        pc = pc + 1;
        if (pc > pn)
            break
        end
    end

    v_n = v;        % The last time step solution before the next step solution. Initial Guess.
    v2_n = v2;
    

    if t >= t_final
        break
    end

    for iterations = 1:iteration_num
    
        % Calculate electrode potential
        % Calculate CDI Current
        current_mag = current_calculator(current_collect_index,eqn_ma,s_mo,z_mo_input,D_eff_f,v,delta_x_f,F,A_c);

        if driving_force == 1
            phi_el_value = (V_app-current_mag*EER)/V_T;
        else
            phi_el_value = phi_el_value+(I_app-current_mag)/C_external/V_T*dt;
        end
        phi_el = (-phi_el_value).*cathode_domain_check(x_c);

        % Updating exp_array
        exp_array = zeros(N*eqn_ma,1);
        
        for i = 1:s_mo
            exp_array(i:eqn_ma:end) = exp( -z_mo_input(i)*(phi_el-v(eqn_ma:eqn_ma:end)-v2(eqn_ma:eqn_mi:end)*F/V_T/C_S)+mu_att(i) );
        end

        %% 1.1. Constructing J - Diagonal sub matrix
        current_index_J = 1;
        for n = 2:N-1
            
            i = eqn_ma*(n-1)+1; % First Index (i,i) for the unit mesh, the top left index in the unit mesh; 
            % i - first element; eqn - s+1, equation num, next unit
            % i - c1,n ; i+s-1 - last specices c,n ; i+s - phi 
            % Constrcuting Diagonal, for transport eq.
            index1 = i:i+s_mo-1;
            index2 = index1;
            element =   p(i:i+s_mo-1)/dt...
                        + D_eff_f(i:i+s_mo-1)./delta_x_f(n)./delta_x_c(n).*...
                        ( 1-z_mo(i:i+s_mo-1).* (v(i+s_mo+eqn_ma)-v(i+s_mo))/2 )...
                        + D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma)./delta_x_f(n-1)./delta_x_c(n).*...
                        ( 1+z_mo(i:i+s_mo-1).* (v(i+s_mo)-v(i+s_mo-eqn_ma))/2 );
            num = length(element);

            index1_J(current_index_J:current_index_J+num-1,1) = index1;
            index2_J(current_index_J:current_index_J+num-1,1) = index2;
            element_J(current_index_J:current_index_J+num-1,1) = element;
            current_index_J = current_index_J+num;

            % Constructing the last column, for transport eq.
            index1 = i:i+s_mo-1;
            index2 = i+eqn_ma-1;
            element =   D_eff_f(i:i+s_mo-1)./delta_x_f(n)./delta_x_c(n).*z_mo(i:i+s_mo-1).*( (v(i+eqn_ma:i+s_mo-1+eqn_ma)+v(i:i+s_mo-1))/2 )+...
                        D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma)./delta_x_f(n-1)./delta_x_c(n).*z_mo(i:i+s_mo-1).*( (v(i-eqn_ma:i+s_mo-1-eqn_ma)+v(i:i+s_mo-1) )/2 );
            num = length(element);

            index1_J(current_index_J:current_index_J+num-1,1) = index1;
            index2_J(current_index_J:current_index_J+num-1,1) = index2;
            element_J(current_index_J:current_index_J+num-1,1) = element;
            current_index_J = current_index_J+num;

            % Constructing the last row, for electroneutrality eq.
            index1 = i+eqn_ma-1;
            index2 = i:i+s_mo-1;
            element = z_mo(i:i+s_mo-1);
            num = length(element);

            index1_J(current_index_J:current_index_J+num-1,1) = index1;
            index2_J(current_index_J:current_index_J+num-1,1) = index2;
            element_J(current_index_J:current_index_J+num-1,1) = element;
            current_index_J = current_index_J+num;
        end
       

        % 1.2. Constructing J - Lower Diagonal sub matrix
        for n = 2:N-1
            
            i = eqn_ma*(n-1)+1; % First Index (i,i) for the unit mesh, the top left index in the unit mesh; 
            % i - first element; eqn - s+1, equation num, next unit
            % i - c1,n ; i+s-1 - last specices c,n ; i+s - phi 
            % Constrcuting Diagonal, for transport eq.
            index1 = i:i+s_mo-1;
            index2 = i-eqn_ma:i+s_mo-1-eqn_ma;
            element =   D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma)./delta_x_f(n-1)/delta_x_c(n).*( -1+z_mo(i:i+s_mo-1)/2.*(v(i+s_mo)-v(i+s_mo-eqn_ma)) ) ...
                        -u/2/delta_x_c(n);
            num = length(element);

            index1_J(current_index_J:current_index_J+num-1,1) = index1;
            index2_J(current_index_J:current_index_J+num-1,1) = index2;
            element_J(current_index_J:current_index_J+num-1,1) = element;
            current_index_J = current_index_J+num;

            % Constructing the last column, for transport eq.
            index1 = i:i+s_mo-1;
            index2 = i-1;
            element = -D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma)./delta_x_f(n-1)./delta_x_c(n).*z_mo(i:i+s_mo-1).*( (v(i-eqn_ma:i+s_mo-1-eqn_ma)+v(i:i+s_mo-1))/2 );
            num = length(element);

            index1_J(current_index_J:current_index_J+num-1,1) = index1;
            index2_J(current_index_J:current_index_J+num-1,1) = index2;
            element_J(current_index_J:current_index_J+num-1,1) = element;
            current_index_J = current_index_J+num;

        end

        % 1.3. Constructing J - Upper Diagonal sub matrix
        for n = 2:N-1
            
            i = eqn_ma*(n-1)+1; % First Index (i,i) for the unit mesh, the top left index in the unit mesh; 
            % i - first element; eqn - s+1, equation num, next unit
            % i - c1,n ; i+s-1 - last specices c,n ; i+s - phi 
            % Constrcuting Diagonal, for transport eq.
            index1 = i:i+s_mo-1;
            index2 = i+eqn_ma:i+s_mo-1+eqn_ma;
            element =   -D_eff_f(i:i+s_mo-1)./delta_x_f(n)./delta_x_c(n).*( 1+z_mo(i:i+s_mo-1)/2.*(v(i+s_mo+eqn_ma)-v(i+s_mo)) ) ...
                        +u/2/delta_x_c(n);
            num = length(element);

            index1_J(current_index_J:current_index_J+num-1,1) = index1;
            index2_J(current_index_J:current_index_J+num-1,1) = index2;
            element_J(current_index_J:current_index_J+num-1,1) = element;
            current_index_J = current_index_J+num;

            % Constructing the last column, for transport eq.
            index1 = i:i+s_mo-1;
            index2 = i+2*eqn_ma-1;
            element = -D_eff_f(i:i+s_mo-1)./delta_x_f(n)./delta_x_c(n).*z_mo(i:i+s_mo-1).*( (v(i+eqn_ma:i+s_mo-1+eqn_ma)+v(i:i+s_mo-1))/2 );
            num = length(element);

            index1_J(current_index_J:current_index_J+num-1,1) = index1;
            index2_J(current_index_J:current_index_J+num-1,1) = index2;
            element_J(current_index_J:current_index_J+num-1,1) = element;
            current_index_J = current_index_J+num;

        end


        % 1.4. Constructing J - Fist and Last Rows Boundary Conidtions
        % Fist Rows - Constrcuting Diagonal resepct to the first unit, for fixed concentration & zero current
        n = 1; 
        i = eqn_ma*(n-1)+1; % First Index (i,i) for the unit mesh, the top left index in the unit mesh; 
        % i - first element; eqn - s+1, equation num, next unit
        % i - c1,n ; i+s-1 - last specices c,n ; i+s - phi 

        index1 = i:i+s_mo;
        index2 = i:i+s_mo;
        element = ones(eqn_ma,1);
        num = length(element);

        index1_J(current_index_J:current_index_J+num-1,1) = index1;
        index2_J(current_index_J:current_index_J+num-1,1) = index2;
        element_J(current_index_J:current_index_J+num-1,1) = element;
        current_index_J = current_index_J+num;
        
        % Last Rows - Constrcuting the last cell, for advective flux only
        n = N; 
        i = eqn_ma*(n-1)+1; % First Index (i,i) for the unit mesh, the top left index in the unit mesh; 
        % i - first element; eqn - s+1, equation num, next unit
        % i - c1,n ; i+s-1 - last specices c,n ; i+s - phi 

        % Constrcuting Diagonal, for transport eq.
        index1 = i:i+s_mo-1;
        index2 = index1;
        element =   p(i:i+s_mo-1)/dt...
                    + u/delta_x_f(n-1)+...
                    + D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma)./delta_x_f(n-1)./delta_x_c(n).*...
                    ( 1+z_mo(i:i+s_mo-1).* (v(i+s_mo)-v(i+s_mo-eqn_ma))/2 );
        num = length(element);
        
        index1_J(current_index_J:current_index_J+num-1,1) = index1;
        index2_J(current_index_J:current_index_J+num-1,1) = index2;
        element_J(current_index_J:current_index_J+num-1,1) = element;
        current_index_J = current_index_J+num;

        % Constructing the last column, for transport eq.
        index1 = i:i+s_mo-1;
        index2 = i+eqn_ma-1;
        element = D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma)./delta_x_f(n-1)./delta_x_c(n).*z_mo(i:i+s_mo-1).*( (v(i-eqn_ma:i+s_mo-1-eqn_ma)+v(i:i+s_mo-1) )/2 );
        num = length(element);

        index1_J(current_index_J:current_index_J+num-1,1) = index1;
        index2_J(current_index_J:current_index_J+num-1,1) = index2;
        element_J(current_index_J:current_index_J+num-1,1) = element;
        current_index_J = current_index_J+num;

        % Constructing the last row, for zero current.
        index1 = i+eqn_ma-1;
        index2 = i:i+s_mo-1;
        element = z_mo(i:i+s_mo-1);
        num = length(element);
        
        index1_J(current_index_J:current_index_J+num-1,1) = index1;
        index2_J(current_index_J:current_index_J+num-1,1) = index2;
        element_J(current_index_J:current_index_J+num-1,1) = element;
        current_index_J = current_index_J+num;

        % Last Rows - Constrcuting the second last unit, for advective flux only
        % Constrcuting Diagonal, for transport eq.
        index1 = i:i+s_mo-1;
        index2 = i-eqn_ma:i+s_mo-1-eqn_ma;
        element = D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma)./delta_x_f(n-1)/delta_x_c(n).*( -1+z_mo(i:i+s_mo-1)/2.*(v(i+s_mo)-v(i+s_mo-eqn_ma)) ) ...
                  - u/delta_x_f(n-1);
        num = length(index1);

        index1_J(current_index_J:current_index_J+num-1,1) = index1;
        index2_J(current_index_J:current_index_J+num-1,1) = index2;
        element_J(current_index_J:current_index_J+num-1,1) = element;
        current_index_J = current_index_J+num;


        % Constructing the last column, for transport eq.
        index1 = i:i+s_mo-1;
        index2 = i-1;
        element = -D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma)./delta_x_f(n-1)./delta_x_c(n).*z_mo(i:i+s_mo-1).*( (v(i-eqn_ma:i+s_mo-1-eqn_ma)+v(i:i+s_mo-1))/2 );
        num = length(element);
    
        index1_J(current_index_J:current_index_J+num-1,1) = index1;
        index2_J(current_index_J:current_index_J+num-1,1) = index2;
        element_J(current_index_J:current_index_J+num-1,1) = element;
        current_index_J = current_index_J+num;

        % Construct J
        J = sparse(index1_J(1:current_index_J-1),index2_J(1:current_index_J-1),element_J(1:current_index_J-1),N*eqn_ma,N*eqn_ma);

        %% 2. Constructing b
        b(1:s_mo,1) = -v(1:s_mo,1)+c_f(1:s_mo);
        b(s_mo+1,1) = -v(s_mo+1,1);

        b(end-eqn_ma+1:end-1,1) = -p(end-eqn_ma+1:end-1).*(v(end-eqn_ma+1:end-1)-v_n(end-eqn_ma+1:end-1))/dt...
                               + u* (v(end-2*eqn_ma+1:end-2*eqn_ma+s_mo,1)-v(end-eqn_ma+1:end-eqn_ma+s_mo,1)) /delta_x_f(end)...
                               -( D_eff_f(end-eqn_ma+1:end-1) .*( (v(end-eqn_ma+1:end-1,1)-v(end-2*eqn_ma+1:end-eqn_ma-1,1))./delta_x_f(end) + ...
                               + z_mo(end-eqn_ma+1:end-1,1).*(v(end-eqn_ma+1:end-1,1)+v(end-2*eqn_ma+1:end-eqn_ma-1,1))/2.*(v(end)-v(end-eqn_ma,1))./(delta_x_f(end)) )...
                               )/delta_x_c(end);


        n = N;
        i = eqn_ma*(n-1)+1;
        b(end,1) = -z_mo(i:i+s_mo-1)' *v(i:i+s_mo-1);
        for n = 2:N-1
         
            i = eqn_ma*(n-1)+1; % First Index for the unit mesh; 
            % i - first element; s+1 - eqn, equation num
            % i - c1,n ; i+s - phi 

            % Constructing the first s terms, transport eq.
            b(i:i+s_mo-1,1) = -p(i:i+s_mo-1).*(v(i:i+s_mo-1,1)-v_n(i:i+s_mo-1,1))/dt ...
                           + ...
                           ( D_eff_f(i:i+s_mo-1) .*( (v(i+eqn_ma:i+s_mo-1+eqn_ma,1)-v(i:i+s_mo-1,1))./delta_x_f(n) + ...
                           + z_mo(i:i+s_mo-1,1).*(v(i+eqn_ma:i+s_mo-1+eqn_ma,1)+v(i:i+s_mo-1,1))/2.*(v(i+s_mo+eqn_ma)-v(i+s_mo,1))./(delta_x_f(n)) )...
                           - u*(v(i+eqn_ma:i+s_mo-1+eqn_ma,1)+v(i:i+s_mo-1,1))/2 ...
                           )/delta_x_c(n) ...
                           - ...
                           ( D_eff_f(i-eqn_ma:i+s_mo-1-eqn_ma) .*( (v(i:i+s_mo-1,1)-v(i-eqn_ma:i+s_mo-1-eqn_ma,1))./delta_x_f(n-1) + ...
                           + z_mo(i:i+s_mo-1,1).*(v(i:i+s_mo-1,1)+v(i-eqn_ma:i+s_mo-1-eqn_ma,1))/2.*(v(i+s_mo)-v(i+s_mo-eqn_ma,1))./(delta_x_f(n-1)) )...
                           - u*(v(i:i+s_mo-1,1)+v(i-eqn_ma:i+s_mo-1-eqn_ma,1))/2 ...
                           )/delta_x_c(n);

            % Constructing the last term, electroneutrality eq.
            b(i+s_mo) = -z_mo(i:i+s_mo-1)' *v(i:i+s_mo-1);    
        end       
        %% 3.1. Constructing J2 - Diagonal sub matrix only
        current_index_J2 = 1;
        for n = 1:N
            
            i = eqn_mi*(n-1)+1; % First Index (i,i) for the unit mesh, the top left index in the unit mesh; 
            % i - first element; eqn - s+1, equation num, next unit
            % i - c_mu_1,n ; i+s-1 - last specices c_mu_s,n ; i+s - sigma_el
            

            % Constructing J2 matrix
            if ~(anode_domain_check(x_c(n)) || cathode_domain_check(x_c(n)))     % Constructing J2 matrix for reservoir & seperator domain
                % Constrcuting Diagonal
                index1 = i:i+s_mo;
                index2 = index1;
                element = ones(s_mo+1,1);
                num = length(element);

                index1_J2(current_index_J2:current_index_J2+num-1,1) = index1;
                index2_J2(current_index_J2:current_index_J2+num-1,1) = index2;
                element_J2(current_index_J2:current_index_J2+num-1,1) = element;
                current_index_J2 = current_index_J2+num;
            else         % Constructing J2 matrix for electrode domain
                % Constructing Diagnoal, for equilibrium eq.
                index1 = i:i+s_mo-1;
                index2 = index1;
                element = ones(s_mo,1);
                num = length(element);

                index1_J2(current_index_J2:current_index_J2+num-1,1) = index1;
                index2_J2(current_index_J2:current_index_J2+num-1,1) = index2;
                element_J2(current_index_J2:current_index_J2+num-1,1) = element;
                current_index_J2 = current_index_J2+num;

                % Constructing the sigma_ele column, for equilibrium eq.
                index1 = i:i+s_mo-1;
                index2 = i+eqn_ma-1;
                diff = (eqn_mi-eqn_ma)*(n-1);
                element =   -z_mo(i-diff:i+s_mo-1-diff).*F/V_T/C_S.*v(i-diff:i+s_mo-1-diff).*exp_array(i-diff:i+s_mo-1-diff);
                num = length(element);
    
                index1_J2(current_index_J2:current_index_J2+num-1,1) = index1;
                index2_J2(current_index_J2:current_index_J2+num-1,1) = index2;
                element_J2(current_index_J2:current_index_J2+num-1,1) = element;
                current_index_J2 = current_index_J2+num;
    
                % Constructing the row for micropore charge balance eq.
                index1 = i+eqn_ma-1;
                index2 = i:i+s_mo;
                element = [z_mo(i-diff:i+s_mo-1-diff);1];
                num = length(element);
    
                index1_J2(current_index_J2:current_index_J2+num-1,1) = index1;
                index2_J2(current_index_J2:current_index_J2+num-1,1) = index2;
                element_J2(current_index_J2:current_index_J2+num-1,1) = element;
                current_index_J2 = current_index_J2+num;
            end
        end
        % Construct J2
        J2 = sparse(index1_J2(1:current_index_J2-1),index2_J2(1:current_index_J2-1),element_J2(1:current_index_J2-1),N*eqn_mi,N*eqn_mi);

        %% 3.2. Constructing J3 - Diagonal sub matrix only
        current_index_J3 = 1;
        for n = 1:N
            
            i = eqn_mi*(n-1)+1; % First Index (i,i) for the unit mesh, the top left index in the unit mesh; 
            % i - first element; eqn - s+1, equation num, next unit
            % i - c_mu_1,n ; i+s-1 - last specices c_mu_s,n ; i+s - sigma_el
            
            % Constructing J3 matrix
            if  anode_domain_check(x_c(n)) || cathode_domain_check(x_c(n))   % Constructing J3 matrix for electrode domain                                                                                   
                
                % Constructing Diagnoal (locally), for equilibrium eq.
                diff = (eqn_mi-eqn_ma)*(n-1);
                index1 = i:i+s_mo-1;
                index2 = index1-diff;
                element = exp_array(i-diff:i+s_mo-1-diff);
                num = length(element);

                index1_J3(current_index_J3:current_index_J3+num-1,1) = index1;
                index2_J3(current_index_J3:current_index_J3+num-1,1) = index2;
                element_J3(current_index_J3:current_index_J3+num-1,1) = element;
                current_index_J3 = current_index_J3+num;

                % Constructing the last column, for equilibrium eq.
                index1 = i:i+s_mo-1;
                index2 = i+eqn_ma-1-diff;
                element = v(i-diff:i+s_mo-1-diff).*exp_array(i-diff:i+s_mo-1-diff).*z_mo(i-diff:i+s_mo-1-diff);
                num = length(element);
    
                index1_J3(current_index_J3:current_index_J3+num-1,1) = index1;
                index2_J3(current_index_J3:current_index_J3+num-1,1) = index2;
                element_J3(current_index_J3:current_index_J3+num-1,1) = element;
                current_index_J3 = current_index_J3+num;
            end
        end
        % Construct J3
        J3 = sparse(index1_J3(1:current_index_J3-1),index2_J3(1:current_index_J3-1),element_J3(1:current_index_J3-1),N*eqn_mi,N*eqn_ma);



        %% 3.3. Constructing b2 
        for n = 1:N
         
            i = eqn_mi*(n-1)+1; % First Index for the unit mesh; 
            % i - first element; eqn - s+1, equation num, next unit
            % i - c_mu_1,n ; i+s-1 - last specices c_mu_s,n ; i+s - sigma_el
            
            % Constructing b2 array
            diff = (eqn_mi-eqn_ma)*(n-1);
            if  ~(anode_domain_check(x_c(n)) || cathode_domain_check(x_c(n)))  % For reservoir/separator domain      
                b2(i:i+s_mo) = -v2(i:i+s_mo);
            else
                % Constructing the first s terms, equilibrium eq.
                b2(i:i+s_mo-1,1) = v(i-diff:i+s_mo-1-diff).*exp_array(i-diff:i+s_mo-1-diff)-v2(i:i+s_mo-1);
                
                % Constructing the term for micropore charge balance eq.
                b2(i+s_mo) = -v2(i+s_mo) -z_mo(i-diff:i+s_mo-1-diff)' *v2(i:i+s_mo-1) -fixed_chem_charge;    
            end
        end
        
        %% 4. Constructing b3, micropore time evolution
        b3=FullMatrix_Extract(v2-v2_n,N,eqn_ma,eqn_mi);


        %% 5. Constructing J2_chem_sum and b2_chem_sum, for surface group
        current_index_J2_chem = 1;

        current_index_J2_chem_Rea = 1;

        b2_chem = zeros(N*eqn_mi,1);
        b2_chem_Rea = zeros(N*eqn_mi,1);
        
        J2_chem_sum = 0;
        b2_chem_sum = 0;

        if s_ch ~= 0
            for n = 1:N
                
                i = eqn_mi*(n-1)+1; % First Index (i,i) for the unit mesh, the top left index in the unit mesh; 
                % i - first element; eqn_ma - s_mo+1, equation num
                % i - c_mu_1,n ; i+s_mo-1 - last mobile specices c_mu_s,n ; i+s_mo - sigma_el;
                % i+s_mo+1 - first chemical group; i+s_mo+s_ch - last chemical group
                
                if  ~(anode_domain_check(x_c(n)) || cathode_domain_check(x_c(n))) % Constructing J2 matrix for reservoir & seperator domain
                    % Constructing J2_chem matrix
                    % Constrcuting Diagonal
                    index1 = i+s_mo+1:i+s_mo+s_ch;
                    index2 = index1;
                    element = ones(s_ch,1);
                    num = length(element);
            
                    index1_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = index1;
                    index2_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = index2;
                    element_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = element;
                    current_index_J2_chem = current_index_J2_chem+num;
            
                    % Constructing b2_chem matrix
                    b2_chem(i+s_mo+1:i+s_mo+s_ch,1) = -v2(i+s_mo+1:i+s_mo+s_ch)+b2_chem(i+s_mo+1:i+s_mo+s_ch,1);
            
                else            % Constructing J2_chem & b2_chem matrix for electrode domain
                    % Constructing J2_chem Diagnoal, for time evolution part
                    diff = (eqn_mi-eqn_ma)*(n-1);
                    index1 = i+s_mo+1:i+s_mo+s_ch;
                    index2 = index1;
                    num = length(index1);
                    element = p_mi(i-diff)/dt*ones(num,1);

            
                    index1_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = index1;
                    index2_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = index2;
                    element_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = element;
                    current_index_J2_chem = current_index_J2_chem+num;
            
                    % Constructing the row for micropore charge balance eq.
                    index1 = i+s_mo;
                    index2 = i+s_mo+1:i+s_mo+s_ch;
                    element = z_chem(1:s_ch);
                    num = length(element);
            
                    index1_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = index1;
                    index2_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = index2;
                    element_J2_chem(current_index_J2_chem:current_index_J2_chem+num-1,1) = element;
                    current_index_J2_chem = current_index_J2_chem+num;

                  
                    % Constructing b2_chem matrix, for charge balance part
                    b2_chem(i+s_mo,1) = -z_chem(1:s_ch)' *v2(i+s_mo+1:i+s_mo+s_ch)+b2_chem(i+s_mo,1);

                    % Constructing b2_chem matrix, for time evolution part
                    b2_chem(i+s_mo+1:i+s_mo+s_ch,1) = -p_mi(i-diff)/dt* (v2(i+s_mo+1:i+s_mo+s_ch)-v2_n(i+s_mo+1:i+s_mo+s_ch))+b2_chem(i+s_mo+1:i+s_mo+s_ch,1);
                    
                    % Reaction Rate for J2_chem_Rea & b2_chem_Rea
                    ii = (n-1)*eqn_mi; % zero index for the present unit
                    for m = 1:min([s_ch,length(k_chem)]) % start to loop in all kinds of surface chemical group reactions
                        
                        Rea_array_L = zeros(1,eqn_mi);
                        for j = R_chem(:,m)'
                            if j ~= 0 
                                temp = 1;
                                for jj = R_chem(:,m)'
                                    if jj ~= j && jj ~= 0
                                        temp = v2(ii+jj)*temp;   
                                    end
                                end
                                Rea_array_L(j) = temp*k_chem(m)+Rea_array_L(j);
                            end
                        end
                    
                        temp = 1;
                        for j = R_chem(:,m)'
                            if j ~= 0 
                                temp = v2(ii+j)*temp;   
                            end
                        end
                        Rea_array_R = temp*k_chem(m);% for b, just an element

                        for j = R_chem(:,m)' % Construct J and b only in a unit cell
                            if j ~= 0 && j ~= index_H
                                index1 = (ii+j)*ones(eqn_mi,1);
                                index2 = ii+(1:eqn_mi)';
                                element = p_mi(ii+1-diff)*Rea_array_L';
                                num = length(element);

                                index1_J2_chem_Rea(current_index_J2_chem_Rea:current_index_J2_chem_Rea+num-1,1) = index1;
                                index2_J2_chem_Rea(current_index_J2_chem_Rea:current_index_J2_chem_Rea+num-1,1) = index2;
                                element_J2_chem_Rea(current_index_J2_chem_Rea:current_index_J2_chem_Rea+num-1,1) = element;
                                current_index_J2_chem_Rea = current_index_J2_chem_Rea+num;%J(i,:)

                            end
                            if j ~= 0 && j ~= index_H
                                b2_chem_Rea(ii+j,1) = -p_mi(ii+1-diff)*Rea_array_R+b2_chem_Rea(ii+j,1);
                            end
                        end
                
                        for j = P_chem(:,m)'
                            if j ~= 0 && j ~= index_H
                                index1 = (ii+j)*ones(eqn_mi,1);
                                index2 = ii+(1:eqn_mi)';
                                element = -p_mi(ii+1-diff)*Rea_array_L';
                                num = length(element);

                                index1_J2_chem_Rea(current_index_J2_chem_Rea:current_index_J2_chem_Rea+num-1,1) = index1;
                                index2_J2_chem_Rea(current_index_J2_chem_Rea:current_index_J2_chem_Rea+num-1,1) = index2;
                                element_J2_chem_Rea(current_index_J2_chem_Rea:current_index_J2_chem_Rea+num-1,1) = element;
                                current_index_J2_chem_Rea = current_index_J2_chem_Rea+num;
   
                            end
                            if j ~= 0 && j ~= index_H
                                b2_chem_Rea(ii+j,1) = p_mi(ii+1-diff)*Rea_array_R+b2_chem_Rea(ii+j,1);
                            end
                        end
                
                    end
                end
            end
           
            % Construct J2_chem
            J2_chem = sparse(index1_J2_chem(1:current_index_J2_chem-1),index2_J2_chem(1:current_index_J2_chem-1),element_J2_chem(1:current_index_J2_chem-1),N*eqn_mi,N*eqn_mi);

            if ~isnan(index1_J2_chem_Rea(1))
                J2_chem_Rea = sparse(index1_J2_chem_Rea(1:current_index_J2_chem_Rea-1),index2_J2_chem_Rea(1:current_index_J2_chem_Rea-1),element_J2_chem_Rea(1:current_index_J2_chem_Rea-1),N*eqn_mi,N*eqn_mi);
            else
                J2_chem_Rea = 0;
                b2_chem_Rea = 0;
            end
            J2_chem_sum = J2_chem+J2_chem_Rea;
            b2_chem_sum = b2_chem+b2_chem_Rea;
        end
        
        %% 6. Faradaic micropore
            %%Anode
            J2_Fara_A = 0;
            b2_Fara_A = 0;
            
            if Faradaic_switch == 1
                current_index_J2_Fara_A = 1;
                b2_Fara_A = zeros(N*eqn_mi,1);
                    
                
                for n = 1:N
                    if  anode_domain_check(x_c(n))  % For anode domain    
                        ii = (n-1)*eqn_mi; % zero index for the present unit
                        diff = (eqn_mi-eqn_ma)*(n-1);
                        for m = 1:Fara_Rea_num_A % start to loop in all kinds of surface reactions
            
                            Rea_array_L_Red = zeros(1,eqn_mi);
                            Rea_array_L_Oxi = zeros(1,eqn_mi);
                            Rea_array_L = zeros(1,eqn_mi);
            
                            temp = 1;
                            Rea_array_R = 0;
            
                            Pi_R = 1;
                            for j = Red_A(:,m)'
                                if j ~= 0
                                    Pi_R = Pi_R*(v2(ii+j)/c_ref);
                                end
                            end
                            Pi_O = 1;
                            for j = Oxi_A(:,m)'
                                if j ~= 0
                                    Pi_O = Pi_O*(v2(ii+j)/c_ref);
                                end
                            end
            
                            sinh_term = sinh(n_ele_A(m)/2* ((F/V_T/C_S)*v2(ii+eqn_ma)-0/V_T+1/n_ele_A(m)*log(Pi_R/Pi_O)) );
                            cosh_term = cosh(n_ele_A(m)/2* ((F/V_T/C_S)*v2(ii+eqn_ma)-0/V_T+1/n_ele_A(m)*log(Pi_R/Pi_O)) );
                            
                            for k = 1:length(Oxi_A(:,m))  % Construct the reaction horizontal vetcor for J & b later
                                j = Oxi_A(k,m); % Reaction species index in local J2 matrix domain
                                if j ~= 0
                                    factor = s_coeff_Oxi_A(k)*(v2(ii+j)/c_ref)^(s_coeff_Oxi_A(k)-1);
                                    for kk = 1:length(Oxi_A(:,m))
                                        jj = Oxi_A(kk,m);
                                        if (jj ~= j) && (jj ~= 0)
                                            factor = factor*(v2(ii+jj)/c_ref)^s_coeff_Oxi_A(kk);
                                        end
                                    end
                                    Rea_array_L_Oxi(j) = factor + Rea_array_L_Oxi(j);
                                end
                            end
                            for k = 1:length(Red_A(:,m))  % Construct the reaction horizontal vetcor for J & b later
                                j = Red_A(k,m); % Reaction species index in local J2 matrix domain
                                if j ~= 0
                                    factor = s_coeff_Red_A(k)*(v2(ii+j)/c_ref)^(s_coeff_Red_A(k)-1);
                                    for kk = 1:length(Red_A(:,m))
                                        jj = Red_A(kk,m);
                                        if (jj ~= j) && (jj ~= 0)
                                            factor = factor*(v2(ii+jj)/c_ref)^s_coeff_Red_A(kk);
                                        end
                                    end
                                    Rea_array_L_Red(j) = factor + Rea_array_L_Red(j);
                                end
                            end
                            Rea_array_L_Red = Rea_array_L_Red * beta_A(m)/2*sqrt(Pi_O/Pi_R)*(sinh_term+cosh_term);
                            Rea_array_L_Oxi = Rea_array_L_Oxi * beta_A(m)/2*sqrt(Pi_R/Pi_O)*(sinh_term-cosh_term);
                            Rea_array_L(eqn_ma) = Rea_array_L(eqn_ma) + beta_A(m)*sqrt(Pi_R*Pi_O)*cosh_term*n_ele_A(m)/2*F/V_T/C_S;   %Sigma_ele contribution
            
                            Rea_array_L = Rea_array_L+Rea_array_L_Red+Rea_array_L_Oxi;
                            Rea_array_R = beta_A(m)*sqrt(Pi_R*Pi_O)*sinh_term + Rea_array_R;
            
                            for k = 1:length(Oxi_A(:,m)) % Construct J and b only in a unit cell
                                j = Oxi_A(k,m);
                                if j > eqn_ma
                                    index1 = (ii+j)*ones(eqn_mi,1);
                                    index2 = ii+(1:eqn_mi)';
                                    element = -p_mi(ii+1-diff)*s_coeff_Oxi_A(k)*Rea_array_L';
                                    num = length(element);
                        
                                    index1_J2_Fara_A(current_index_J2_Fara_A:current_index_J2_Fara_A+num-1,1) = index1;
                                    index2_J2_Fara_A(current_index_J2_Fara_A:current_index_J2_Fara_A+num-1,1) = index2;
                                    element_J2_Fara_A(current_index_J2_Fara_A:current_index_J2_Fara_A+num-1,1) = element;
                                    current_index_J2_Fara_A = current_index_J2_Fara_A+num; 
            
                                    b2_Fara_A(ii+j) = p_mi(ii+1-diff)*s_coeff_Oxi_A(k)*Rea_array_R + b2_Fara_A(ii+j);
                                end
                            end
                            
                            for k = 1:length(Red_A(:,m))
                                j = Red_A(k,m);
                                if j > eqn_ma
                                    index1 = (ii+j)*ones(eqn_mi,1);
                                    index2 = ii+(1:eqn_mi)';
                                    element = p_mi(ii+1-diff)*s_coeff_Red_A(k)*Rea_array_L';
                                    num = length(element);
                    
                                    index1_J2_Fara_A(current_index_J2_Fara_A:current_index_J2_Fara_A+num-1,1) = index1;
                                    index2_J2_Fara_A(current_index_J2_Fara_A:current_index_J2_Fara_A+num-1,1) = index2;
                                    element_J2_Fara_A(current_index_J2_Fara_A:current_index_J2_Fara_A+num-1,1) = element;
                                    current_index_J2_Fara_A = current_index_J2_Fara_A+num; 
            
                                    b2_Fara_A(ii+j) = -p_mi(ii+1-diff)*s_coeff_Red_A(k)*Rea_array_R + b2_Fara_A(ii+j);
                                end
                            end
                    
                        end
                    end
                end
                if ~isnan(index1_J2_Fara_A(1))
                    J2_Fara_A = sparse(index1_J2_Fara_A(1:current_index_J2_Fara_A-1),index2_J2_Fara_A(1:current_index_J2_Fara_A-1),element_J2_Fara_A(1:current_index_J2_Fara_A-1),N*eqn_mi,N*eqn_mi);
                else
                    J2_Fara_A = 0;
                    b2_Fara_A = 0;
                end
            end
            %%Cathode
            J2_Fara_C = 0;
            b2_Fara_C = 0;
            
            if Faradaic_switch == 1
                current_index_J2_Fara_C = 1;
                b2_Fara_C = zeros(N*eqn_mi,1);
                    
                
                for n = 1:N
                    if  cathode_domain_check(x_c(n))  % For cathode domain    
                        ii = (n-1)*eqn_mi; % zero index for the present unit
                        diff = (eqn_mi-eqn_ma)*(n-1);
                        for m = 1:Fara_Rea_num_C % start to loop in all kinds of surface reactions
            
                            Rea_array_L_Red = zeros(1,eqn_mi);
                            Rea_array_L_Oxi = zeros(1,eqn_mi);
                            Rea_array_L = zeros(1,eqn_mi);
            
                            temp = 1;
                            Rea_array_R = 0;
            
                            Pi_R = 1;
                            for j = Red_C(:,m)'
                                if j ~= 0
                                    Pi_R = Pi_R*(v2(ii+j)/c_ref);
                                end
                            end
                            Pi_O = 1;
                            for j = Oxi_C(:,m)'
                                if j ~= 0
                                    Pi_O = Pi_O*(v2(ii+j)/c_ref);
                                end
                            end
            
                            sinh_term = sinh(n_ele_C(m)/2* ((F/V_T/C_S)*v2(ii+eqn_ma)-OCV_rel_dimensional/V_T+1/n_ele_C(m)*log(Pi_R/Pi_O)) );
                            cosh_term = cosh(n_ele_C(m)/2* ((F/V_T/C_S)*v2(ii+eqn_ma)-OCV_rel_dimensional/V_T+1/n_ele_C(m)*log(Pi_R/Pi_O)) );
                        
                            for k = 1:length(Oxi_C(:,m))  % Construct the reaction horizontal vetcor for J & b later
                                j = Oxi_C(k,m); % Reaction species index in local J2 matrix domain
                                if j ~= 0
                                    factor = s_coeff_Oxi_C(k)*(v2(ii+j)/c_ref)^(s_coeff_Oxi_C(k)-1);
                                    for kk = 1:length(Oxi_C(:,m))
                                        jj = Oxi_C(kk,m);
                                        if (jj ~= j) && (jj ~= 0)
                                            factor = factor*(v2(ii+jj)/c_ref)^s_coeff_Oxi_C(kk);
                                        end
                                    end
                                    Rea_array_L_Oxi(j) = factor + Rea_array_L_Oxi(j);
                                end
                            end
                            for k = 1:length(Red_C(:,m))  % Construct the reaction horizontal vetcor for J & b later
                                j = Red_C(k,m); % Reaction species index in local J2 matrix domain
                                if j ~= 0
                                    factor = s_coeff_Red_C(k)*(v2(ii+j)/c_ref)^(s_coeff_Red_C(k)-1);
                                    for kk = 1:length(Red_C(:,m))
                                        jj = Red_C(kk,m);
                                        if (jj ~= j) && (jj ~= 0)
                                            factor = factor*(v2(ii+jj)/c_ref)^s_coeff_Red_C(kk);
                                        end
                                    end
                                    Rea_array_L_Red(j) = factor + Rea_array_L_Red(j);
                                end
                            end
                            Rea_array_L_Red = Rea_array_L_Red * beta_C(m)/2*sqrt(Pi_O/Pi_R)*(sinh_term+cosh_term);
                            Rea_array_L_Oxi = Rea_array_L_Oxi * beta_C(m)/2*sqrt(Pi_R/Pi_O)*(sinh_term-cosh_term);
                            Rea_array_L(eqn_ma) = Rea_array_L(eqn_ma) + beta_C(m)*sqrt(Pi_R*Pi_O)*cosh_term*n_ele_C(m)/2*F/V_T/C_S;   %Sigma_ele contribution
            
                            Rea_array_L = Rea_array_L+Rea_array_L_Red+Rea_array_L_Oxi;
                            Rea_array_R = beta_C(m)*sqrt(Pi_R*Pi_O)*sinh_term + Rea_array_R;
            
                            for k = 1:length(Oxi_C(:,m))
                                j = Oxi_C(k,m); % Construct J and b only in a unit cell
                                if j > eqn_ma
                                    index1 = (ii+j)*ones(eqn_mi,1);
                                    index2 = ii+(1:eqn_mi)';
                                    element = -p_mi(ii+1-diff)*s_coeff_Oxi_C(k)*Rea_array_L';
                                    num = length(element);
                        
                                    index1_J2_Fara_C(current_index_J2_Fara_C:current_index_J2_Fara_C+num-1,1) = index1;
                                    index2_J2_Fara_C(current_index_J2_Fara_C:current_index_J2_Fara_C+num-1,1) = index2;
                                    element_J2_Fara_C(current_index_J2_Fara_C:current_index_J2_Fara_C+num-1,1) = element;
                                    current_index_J2_Fara_C = current_index_J2_Fara_C+num; 
            
                                    b2_Fara_C(ii+j) = p_mi(ii+1-diff)*s_coeff_Oxi_C(k)*Rea_array_R + b2_Fara_C(ii+j);
                                end
                            end
                            
                            for k = 1:length(Red_C(:,m))
                                j = Red_C(k,m);
                                if j > eqn_ma
                                    index1 = (ii+j)*ones(eqn_mi,1);
                                    index2 = ii+(1:eqn_mi)';
                                    element = p_mi(ii+1-diff)*s_coeff_Red_C(k)*Rea_array_L';
                                    num = length(element);
                    
                                    index1_J2_Fara_C(current_index_J2_Fara_C:current_index_J2_Fara_C+num-1,1) = index1;
                                    index2_J2_Fara_C(current_index_J2_Fara_C:current_index_J2_Fara_C+num-1,1) = index2;
                                    element_J2_Fara_C(current_index_J2_Fara_C:current_index_J2_Fara_C+num-1,1) = element;
                                    current_index_J2_Fara_C = current_index_J2_Fara_C+num; 
            
                                    b2_Fara_C(ii+j) = -p_mi(ii+1-diff)*s_coeff_Red_C(k)*Rea_array_R + b2_Fara_C(ii+j);
                                end
                            end
                    
                        end
                    end
                end
                if ~isnan(index1_J2_Fara_C(1))
                    J2_Fara_C = sparse(index1_J2_Fara_C(1:current_index_J2_Fara_C-1),index2_J2_Fara_C(1:current_index_J2_Fara_C-1),element_J2_Fara_C(1:current_index_J2_Fara_C-1),N*eqn_mi,N*eqn_mi);
                else
                    J2_Fara_C = 0;
                    b2_Fara_C = 0;
                end
            end

        %% 7. Preparing mapping from local micro matrix to the macro matrix
        J2_sum = J2+J2_chem_sum+J2_Fara_A+J2_Fara_C;
        Matrix_temp1_full = J2_sum\J3;
        Matrix_temp1 = Matrix_temp1_full;
        Matrix_temp1(s_mo+1:eqn_mi:end,:) = 0;        % Only adding the entries for transport equations     
        Matrix_temp1=FullMatrix_Extract(Matrix_temp1,N,eqn_ma,eqn_mi);  % Will be used for mapping
 
        b2_sum = b2+b2_chem_sum+b2_Fara_A+b2_Fara_C;
        Vector_temp1_full = J2_sum\b2_sum;
        Vector_temp1 = Vector_temp1_full;
        Vector_temp1(s_mo+1:eqn_mi:end,:) = 0;
        Vector_temp1=FullMatrix_Extract(Vector_temp1,N,eqn_ma,eqn_mi);  % Will be used for mapping
        
        b3(s_mo+1:eqn_ma:end,:) = 0;  % Will be used for mapping
        
        %% 8. Construct Bulk Reaction J and b
        % Initialization 
%         index1_J_Bulk_Rea = nan(2*(length(R_bulk)^2)*(N*eqn_ma)^2,1);
%         index2_J_Bulk_Rea = index1_J_Bulk_Rea;
%         element_J_Bulk_Rea = index1_J_Bulk_Rea;
        current_index_J_Bulk_Rea = 1;
        b_Bulk_Rea = zeros(N*eqn_ma,1);
        J_Bulk_Rea = 0;
        if Bulk_Reaction_switch == 1
            for n = 2:N
                ii = (n-1)*eqn_ma; % zero index for the present unit
                for m = 1:length(k_ma) % start to loop in all kinds of reactions
                    Rea_array_L = zeros(1,eqn_ma);

                    for j = R_bulk(:,m)'
                        if j ~= 0 
                            temp = 1;
                            for jj = R_bulk(:,m)'
                                if jj ~= j && jj ~= 0
                                    temp = v(ii+jj)*temp;   
                                end
                            end
                            Rea_array_L(j) = temp*k_ma(m)+Rea_array_L(j);
                        end
                    end
                
                    temp = 1;
                    for j = R_bulk(:,m)'
                        if j ~= 0 
                            temp = v(ii+j)*temp;   
                        end
                    end
                    Rea_array_R = temp*k_ma(m);% for b, just an element
            
                    for j = R_bulk(:,m)' % Construct J and b only in a unit cell
                        if j ~= 0
                            index1 = (ii+j)*ones(eqn_ma,1);
                            index2 = ii+(1:eqn_ma)';
                            element = p(ii+j)*Rea_array_L';
                            num = length(element);
    
                            index1_J_Bulk_Rea(current_index_J_Bulk_Rea:current_index_J_Bulk_Rea+num-1,1) = index1;
                            index2_J_Bulk_Rea(current_index_J_Bulk_Rea:current_index_J_Bulk_Rea+num-1,1) = index2;
                            element_J_Bulk_Rea(current_index_J_Bulk_Rea:current_index_J_Bulk_Rea+num-1,1) = element;%J(i,:)
                            current_index_J_Bulk_Rea = current_index_J_Bulk_Rea+num;
                        end
                        if j ~= 0
                            b_Bulk_Rea(ii+j) = -p(ii+j)*Rea_array_R+b_Bulk_Rea(ii+j);
                        end
                    end
            
                    for j = P_bulk(:,m)'
                        if j ~= 0
    
                            index1 = (ii+j)*ones(eqn_ma,1);
                            index2 = ii+(1:eqn_ma)';
                            element = -p(ii+j)*Rea_array_L';
                            num = length(element);
    
                            index1_J_Bulk_Rea(current_index_J_Bulk_Rea:current_index_J_Bulk_Rea+num-1,1) = index1;
                            index2_J_Bulk_Rea(current_index_J_Bulk_Rea:current_index_J_Bulk_Rea+num-1,1) = index2;
                            element_J_Bulk_Rea(current_index_J_Bulk_Rea:current_index_J_Bulk_Rea+num-1,1) = element;%J(i,:)
                            current_index_J_Bulk_Rea = current_index_J_Bulk_Rea+num;       
                        end
                        if j ~= 0
                            b_Bulk_Rea(ii+j) = p(ii+j)*Rea_array_R+b_Bulk_Rea(ii+j);
                        end
                    end
    
                end
            end
            if ~isnan(index1_J_Bulk_Rea(1))
                J_Bulk_Rea = sparse(index1_J_Bulk_Rea(1:current_index_J_Bulk_Rea-1),index2_J_Bulk_Rea(1:current_index_J_Bulk_Rea-1),element_J_Bulk_Rea(1:current_index_J_Bulk_Rea-1),N*eqn_ma,N*eqn_ma);
            else
                J_Bulk_Rea = 0;
                b_Bulk_Rea = 0;
            end
        end

        %% 9. Building Micropore Bulk Reaction Matrix for J
        current_index_J_Bulk_Rea_mi = 1;
        J_Bulk_Rea_mi = 0;
        b_Bulk_Rea_mi = zeros(N*eqn_ma,1);

        if Bulk_Reaction_switch == 1
            for n = 1:N
                diff = (eqn_mi-eqn_ma)*(n-1);
                if  anode_domain_check(x_c(n)) || cathode_domain_check(x_c(n))   % For electrode domain    
                    ii = (n-1)*eqn_ma; % zero index for the present unit
                    for m = 1:length(k_mi) % start to loop in all kinds of reactions
                        Rea_array_L = zeros(1,eqn_ma*N);
                        Rea_array_R = 0;

                        for j = R_bulk_mi(:,m)'
                            if j ~= 0 
                                temp = 1;
                                for jj = R_bulk_mi(:,m)'
                                    if jj ~= j && jj ~= 0
                                        temp = v2(ii+jj+diff)*temp;   
                                    end
                                end
                                Rea_array_L(1:eqn_ma*N) = temp*k_mi(m).*Matrix_temp1_full(ii+j+diff,:)+Rea_array_L(1:eqn_ma*N);
                                Rea_array_R = Rea_array_R+temp*k_mi(m).*Vector_temp1_full(ii+j+diff);
                            end
                        end
                    
                        temp = 1;
                        for j = R_bulk_mi(:,m)'
                            if j ~= 0 
                                temp = v2(ii+j+diff)*temp;   
                            end
                        end
                        Rea_array_R = Rea_array_R+temp*k_mi(m);% for b, just an element


                        for j = R_bulk_mi(:,m)' % Construct J and b only in a unit cell
                            if j ~= 0
                                index1 = (ii+j)*ones(eqn_ma*N,1);
                                index2 = (1:eqn_ma*N)';
                                element = p_mi(ii+j)*Rea_array_L';
                                num = length(element);
                
                                index1_J_Bulk_Rea_mi(current_index_J_Bulk_Rea_mi:current_index_J_Bulk_Rea_mi+num-1,1) = index1;
                                index2_J_Bulk_Rea_mi(current_index_J_Bulk_Rea_mi:current_index_J_Bulk_Rea_mi+num-1,1) = index2;
                                element_J_Bulk_Rea_mi(current_index_J_Bulk_Rea_mi:current_index_J_Bulk_Rea_mi+num-1,1) = element;
                                current_index_J_Bulk_Rea_mi = current_index_J_Bulk_Rea_mi+num;
                            end
                            if j ~= 0
                                b_Bulk_Rea_mi(ii+j) = -p_mi(ii+j)*Rea_array_R+b_Bulk_Rea_mi(ii+j);
                            end
                        end
                
                        for j = P_bulk_mi(:,m)'
                            if j ~= 0
                                index1 = (ii+j)*ones(eqn_ma*N,1);
                                index2 = (1:eqn_ma*N)';
                                element = -p_mi(ii+j)*Rea_array_L';
                                num = length(element);
                            
                                index1_J_Bulk_Rea_mi(current_index_J_Bulk_Rea_mi:current_index_J_Bulk_Rea_mi+num-1,1) = index1;
                                index2_J_Bulk_Rea_mi(current_index_J_Bulk_Rea_mi:current_index_J_Bulk_Rea_mi+num-1,1) = index2;
                                element_J_Bulk_Rea_mi(current_index_J_Bulk_Rea_mi:current_index_J_Bulk_Rea_mi+num-1,1) = element;
                                current_index_J_Bulk_Rea_mi = current_index_J_Bulk_Rea_mi+num;      
                            end
                            if j ~= 0
                                b_Bulk_Rea_mi(ii+j) = p_mi(ii+j)*Rea_array_R+b_Bulk_Rea_mi(ii+j);
                            end
                        end
                    end
                end
            end
            if ~isnan(index1_J_Bulk_Rea_mi(1))
                J_Bulk_Rea_mi = sparse(index1_J_Bulk_Rea_mi(1:current_index_J_Bulk_Rea_mi-1),index2_J_Bulk_Rea_mi(1:current_index_J_Bulk_Rea_mi-1),element_J_Bulk_Rea_mi(1:current_index_J_Bulk_Rea_mi-1),N*eqn_ma,N*eqn_ma);
            else
                J_Bulk_Rea_mi = 0;
                b_Bulk_Rea_mi = 0;
            end
        end
       
        %% 10.Building Micropore chemical group Reaction for J & b, gives Matrix J_chem, b_chem
        J_chem_Rea = 0;
        b_chem_Rea = 0;
        if s_ch ~= 0
            current_index_J_chem_Rea = 1;
            b_chem_Rea = zeros(N*eqn_ma,1);
                
            for n = 1:N
                if  anode_domain_check(x_c(n)) || cathode_domain_check(x_c(n))   % For electrode domain    
                    ii = (n-1)*eqn_mi; % zero index for the present unit
                    diff = (eqn_mi-eqn_ma)*(n-1);
                    for m = 1:min([s_ch,length(k_chem)]) % start to loop in all kinds of surface reactions

                        Rea_array_L = zeros(1,eqn_ma*N);
                        Rea_array_R = 0;

                        for j = R_chem(:,m)'
                            if j ~= 0 
                                temp = 1;
                                for jj = R_chem(:,m)'
                                    if jj ~= j && jj ~= 0
                                        temp = v2(ii+jj)*temp;   
                                    end
                                end
                                Rea_array_L(1:eqn_ma*N) = temp*k_chem(m).*Matrix_temp1_full(ii+j,:)+Rea_array_L(1:eqn_ma*N);
                                Rea_array_R = Rea_array_R+temp*k_chem(m).*Vector_temp1_full(ii+j);
                            end
                        end
                    
                        temp = 1;
                        for j = R_chem(:,m)'
                            if j ~= 0 
                                temp = v2(ii+j)*temp;   
                            end
                        end
                        Rea_array_R = Rea_array_R+temp*k_chem(m);% for b, just an element

                        for j = R_chem(:,m)' % Construct J and b only in a unit cell
                            if j == index_H
                                index1 = (ii+j)*ones(eqn_ma*N,1)-diff;
                                index2 = (1:eqn_ma*N)';
                                element = p_mi(ii+1-diff)*Rea_array_L';
                                num = length(element);
                    
                                index1_J_chem_Rea(current_index_J_chem_Rea:current_index_J_chem_Rea+num-1,1) = index1;
                                index2_J_chem_Rea(current_index_J_chem_Rea:current_index_J_chem_Rea+num-1,1) = index2;
                                element_J_chem_Rea(current_index_J_chem_Rea:current_index_J_chem_Rea+num-1,1) = element;
                                current_index_J_chem_Rea = current_index_J_chem_Rea+num; 
                            end
                            if j == index_H
                                b_chem_Rea((ii+j)-diff) = -p_mi(ii+1-diff)*Rea_array_R + b_chem_Rea((ii+j)-diff);
                            end
                        end
                
                        for j = P_chem(:,m)'
                            if j == index_H
                                index1 = (ii+j)*ones(eqn_ma*N,1)-diff;
                                index2 = (1:eqn_ma*N)';
                                element = -p_mi(ii+1-diff)*Rea_array_L';
                                num = length(element);
                
                                index1_J_chem_Rea(current_index_J_chem_Rea:current_index_J_chem_Rea+num-1,1) = index1;
                                index2_J_chem_Rea(current_index_J_chem_Rea:current_index_J_chem_Rea+num-1,1) = index2;
                                element_J_chem_Rea(current_index_J_chem_Rea:current_index_J_chem_Rea+num-1,1) = element;
                                current_index_J_chem_Rea = current_index_J_chem_Rea+num; 
                            end
                            if j == index_H
                                b_chem_Rea((ii+j)-diff) = p_mi(ii+1-diff)*Rea_array_R + b_chem_Rea((ii+j)-diff);
                            end
                        end
                
                    end
                end
            end
            if ~isnan(index1_J_chem_Rea(1))
                J_chem_Rea = sparse(index1_J_chem_Rea(1:current_index_J_chem_Rea-1),index2_J_chem_Rea(1:current_index_J_chem_Rea-1),element_J_chem_Rea(1:current_index_J_chem_Rea-1),N*eqn_ma,N*eqn_ma);
            else
                J_chem_Rea = 0;
                b_chem_Rea = 0;
            end
        end
        %% 11. Faradaic macropore
        J_Fara_A = 0;
        b_Fara_A = 0;

        if Faradaic_switch ~= 0
            current_index_J_Fara_A = 1;
            b_Fara_A = zeros(N*eqn_ma,1);
                
            for n = 1:N
                if  anode_domain_check(x_c(n))  % For anode domain    
                    ii = (n-1)*eqn_mi; % zero index for the present unit
                    diff = (eqn_mi-eqn_ma)*(n-1);
                    for m = 1:Fara_Rea_num_A % start to loop in all kinds of surface reactions
        
                        Rea_array_L_Red = zeros(1,eqn_ma*N);
                        Rea_array_L_Oxi = zeros(1,eqn_ma*N);
                        Rea_array_L = zeros(1,eqn_ma*N);
        
                        temp = 1;
                        Rea_array_R = 0;
        
                        Pi_R = 1;
                        for j = Red_A(:,m)'
                            if j ~= 0
                                Pi_R = Pi_R*(v2(ii+j)/c_ref);
                            end
                        end
                        Pi_O = 1;
                        for j = Oxi_A(:,m)'
                            if j ~= 0
                                Pi_O = Pi_O*(v2(ii+j)/c_ref);
                            end
                        end
        
                        sinh_term = sinh(n_ele_A(m)/2* ((F/V_T/C_S)*v2(ii+eqn_ma)-0/V_T+1/n_ele_A(m)*log(Pi_R/Pi_O)) );
                        cosh_term = cosh(n_ele_A(m)/2* ((F/V_T/C_S)*v2(ii+eqn_ma)-0/V_T+1/n_ele_A(m)*log(Pi_R/Pi_O)) );
                       
                        for k = 1:length(Oxi_A(:,m))  % Construct the reaction horizontal vetcor for J & b later
                            j = Oxi_A(k,m); % Reaction species index in local J2 matrix domain
                            if j ~= 0
                                factor = s_coeff_Oxi_A(k)*(v2(ii+j)/c_ref)^(s_coeff_Oxi_A(k)-1);
                                for kk = 1:length(Oxi_A(:,m))
                                    jj = Oxi_A(kk,m);
                                    if (jj ~= j) && (jj ~= 0)
                                        factor = factor*(v2(ii+jj)/c_ref)^s_coeff_Oxi_A(kk);
                                    end
                                end
                                Rea_array_L_Oxi = factor.*Matrix_temp1_full(ii+j,:) + Rea_array_L_Oxi;
                            end
                        end
                        for k = 1:length(Red_A(:,m))  % Construct the reaction horizontal vetcor for J & b later
                            j = Red_A(k,m); % Reaction species index in local J2 matrix domain
                            if j ~= 0
                                factor = s_coeff_Red_A(k)*(v2(ii+j)/c_ref)^(s_coeff_Red_A(k)-1);
                                for kk = 1:length(Red_A(:,m))
                                    jj = Red_A(kk,m);
                                    if (jj ~= j) && (jj ~= 0)
                                        factor = factor*(v2(ii+jj)/c_ref)^s_coeff_Red_A(kk);
                                    end
                                end
                                Rea_array_L_Red = factor.*Matrix_temp1_full(ii+j,:) + Rea_array_L_Red;
                            end
                        end
                        Rea_array_L_Red = Rea_array_L_Red * beta_A(m)/2*sqrt(Pi_O/Pi_R)*(sinh_term+cosh_term);
                        Rea_array_L_Oxi = Rea_array_L_Oxi * beta_A(m)/2*sqrt(Pi_R/Pi_O)*(sinh_term-cosh_term);
                        Rea_array_L = Rea_array_L + beta_A(m)*sqrt(Pi_R*Pi_O)*cosh_term*n_ele_A(m)/2*F/V_T/C_S.*Matrix_temp1_full(ii+eqn_ma,:);   %Sigma_ele contribution
        
                        Rea_array_L = Rea_array_L+Rea_array_L_Red+Rea_array_L_Oxi;
                        Rea_array_R = beta_A(m)*sqrt(Pi_R*Pi_O)*sinh_term + Rea_array_R;
        
                        for k = 1:length(Oxi_A(:,m))
                            j = Oxi_A(k,m); % Construct J and b only in a unit cell
                            if (j < eqn_ma) && (j ~= 0)
                                index1 = (ii+j)*ones(eqn_ma*N,1)-diff;
                                index2 = (1:eqn_ma*N)';
                                element = -p_mi(ii+1-diff)*s_coeff_Oxi_A(k)*Rea_array_L';
                                num = length(element);
                    
                                index1_J_Fara_A(current_index_J_Fara_A:current_index_J_Fara_A+num-1,1) = index1;
                                index2_J_Fara_A(current_index_J_Fara_A:current_index_J_Fara_A+num-1,1) = index2;
                                element_J_Fara_A(current_index_J_Fara_A:current_index_J_Fara_A+num-1,1) = element;
                                current_index_J_Fara_A = current_index_J_Fara_A+num; 
        
                                b_Fara_A(ii+j-diff) = p_mi(ii+1-diff)*s_coeff_Oxi_A(k)*Rea_array_R + b_Fara_A(ii+j-diff);
                            end
                        end

                        for k = 1:length(Red_A(:,m))
                            j = Red_A(k,m);
                            if (j < eqn_ma) && (j ~= 0)
                                index1 = (ii+j)*ones(eqn_ma*N,1)-diff;
                                index2 = (1:eqn_ma*N)';
                                element = p_mi(ii+1-diff)*s_coeff_Red_A(k)*Rea_array_L';
                                num = length(element);
                
                                index1_J_Fara_A(current_index_J_Fara_A:current_index_J_Fara_A+num-1,1) = index1;
                                index2_J_Fara_A(current_index_J_Fara_A:current_index_J_Fara_A+num-1,1) = index2;
                                element_J_Fara_A(current_index_J_Fara_A:current_index_J_Fara_A+num-1,1) = element;
                                current_index_J_Fara_A = current_index_J_Fara_A+num; 
        
                                b_Fara_A(ii+j-diff) = -p_mi(ii+1-diff)*s_coeff_Red_A(k)*Rea_array_R + b_Fara_A(ii+j-diff);
                            end
                        end
                
                    end
                end
            end
            if ~isnan(index1_J_Fara_A(1))
                J_Fara_A = sparse(index1_J_Fara_A(1:current_index_J_Fara_A-1),index2_J_Fara_A(1:current_index_J_Fara_A-1),element_J_Fara_A(1:current_index_J_Fara_A-1),N*eqn_ma,N*eqn_ma);
            else
                J_Fara_A = 0;
                b_Fara_A = 0;
            end
        end

        J_Fara_C = 0;
        b_Fara_C = 0;
        
        if Faradaic_switch ~= 0
            current_index_J_Fara_C = 1;
            b_Fara_C = zeros(N*eqn_ma,1);
                
            for n = 1:N
                if  cathode_domain_check(x_c(n))  % For cathode domain    
                    ii = (n-1)*eqn_mi; % zero index for the present unit
                    diff = (eqn_mi-eqn_ma)*(n-1);
                    for m = 1:Fara_Rea_num_C % start to loop in all kinds of surface reactions
        
                        Rea_array_L_Red = zeros(1,eqn_ma*N);
                        Rea_array_L_Oxi = zeros(1,eqn_ma*N);
                        Rea_array_L = zeros(1,eqn_ma*N);
        
                        temp = 1;
                        Rea_array_R = 0;
        
                        Pi_R = 1;
                        for j = Red_C(:,m)'
                            if j ~= 0
                                Pi_R = Pi_R*(v2(ii+j)/c_ref);
                            end
                        end
                        Pi_O = 1;
                        for j = Oxi_C(:,m)'
                            if j ~= 0
                                Pi_O = Pi_O*(v2(ii+j)/c_ref);
                            end
                        end
        
                        sinh_term = sinh(n_ele_C(m)/2* ((F/V_T/C_S)*v2(ii+eqn_ma)-OCV_rel_dimensional/V_T+1/n_ele_C(m)*log(Pi_R/Pi_O)) );
                        cosh_term = cosh(n_ele_C(m)/2* ((F/V_T/C_S)*v2(ii+eqn_ma)-OCV_rel_dimensional/V_T+1/n_ele_C(m)*log(Pi_R/Pi_O)) );
                       
                        for k = 1:length(Oxi_C(:,m))  % Construct the reaction horizontal vetcor for J & b later
                            j = Oxi_C(k,m); % Reaction species index in local J2 matrix domain
                            if j ~= 0
                                factor = s_coeff_Oxi_C(k)*(v2(ii+j)/c_ref)^(s_coeff_Oxi_C(k)-1);
                                for kk = 1:length(Oxi_C(:,m))
                                    jj = Oxi_C(kk,m);
                                    if (jj ~= j) && (jj ~= 0)
                                        factor = factor*(v2(ii+jj)/c_ref)^s_coeff_Oxi_C(kk);
                                    end
                                end
                                Rea_array_L_Oxi = factor.*Matrix_temp1_full(ii+j,:) + Rea_array_L_Oxi;
                            end
                        end
                        for k = 1:length(Red_C(:,m))  % Construct the reaction horizontal vetcor for J & b later
                            j = Red_C(k,m); % Reaction species index in local J2 matrix domain
                            if j ~= 0
                                factor = s_coeff_Red_C(k)*(v2(ii+j)/c_ref)^(s_coeff_Red_C(k)-1);
                                for kk = 1:length(Red_C(:,m))
                                    jj = Red_C(kk,m);
                                    if (jj ~= j) && (jj ~= 0)
                                        factor = factor*(v2(ii+jj)/c_ref)^s_coeff_Red_C(kk);
                                    end
                                end
                                Rea_array_L_Red = factor.*Matrix_temp1_full(ii+j,:) + Rea_array_L_Red;
                            end
                        end
                        Rea_array_L_Red = Rea_array_L_Red * beta_C(m)/2*sqrt(Pi_O/Pi_R)*(sinh_term+cosh_term);
                        Rea_array_L_Oxi = Rea_array_L_Oxi * beta_C(m)/2*sqrt(Pi_R/Pi_O)*(sinh_term-cosh_term);
                        Rea_array_L = Rea_array_L + beta_C(m)*sqrt(Pi_R*Pi_O)*cosh_term*n_ele_C(m)/2*F/V_T/C_S.*Matrix_temp1_full(ii+eqn_ma,:);   %Sigma_ele contribution
        
                        Rea_array_L = Rea_array_L+Rea_array_L_Red+Rea_array_L_Oxi;
                        Rea_array_R = beta_C(m)*sqrt(Pi_R*Pi_O)*sinh_term + Rea_array_R;
        
                        for k = 1:length(Oxi_C(:,m))
                            j = Oxi_C(k,m); % Construct J and b only in a unit cell
                            if (j < eqn_ma) && (j ~= 0)
                                index1 = (ii+j)*ones(eqn_ma*N,1)-diff;
                                index2 = (1:eqn_ma*N)';
                                element = -p_mi(ii+1-diff)*s_coeff_Oxi_C(k)*Rea_array_L';
                                num = length(element);
                    
                                index1_J_Fara_C(current_index_J_Fara_C:current_index_J_Fara_C+num-1,1) = index1;
                                index2_J_Fara_C(current_index_J_Fara_C:current_index_J_Fara_C+num-1,1) = index2;
                                element_J_Fara_C(current_index_J_Fara_C:current_index_J_Fara_C+num-1,1) = element;
                                current_index_J_Fara_C = current_index_J_Fara_C+num; 
        
                                b_Fara_C(ii+j-diff) = p_mi(ii+1-diff)*s_coeff_Oxi_C(k)*Rea_array_R + b_Fara_C(ii+j-diff);
                            end
                        end

                        for k = 1:length(Red_C(:,m))
                            j = Red_C(k,m);
                            if (j < eqn_ma) && (j ~= 0)
                                index1 = (ii+j)*ones(eqn_ma*N,1)-diff;
                                index2 = (1:eqn_ma*N)';
                                element = p_mi(ii+1-diff)*s_coeff_Red_C(k)*Rea_array_L';
                                num = length(element);
                
                                index1_J_Fara_C(current_index_J_Fara_C:current_index_J_Fara_C+num-1,1) = index1;
                                index2_J_Fara_C(current_index_J_Fara_C:current_index_J_Fara_C+num-1,1) = index2;
                                element_J_Fara_C(current_index_J_Fara_C:current_index_J_Fara_C+num-1,1) = element;
                                current_index_J_Fara_C = current_index_J_Fara_C+num; 
        
                                b_Fara_C(ii+j-diff) = -p_mi(ii+1-diff)*s_coeff_Red_C(k)*Rea_array_R + b_Fara_C(ii+j-diff);
                            end
                        end
                
                    end
                end
            end
            if ~isnan(index1_J_Fara_C(1))
                J_Fara_C = sparse(index1_J_Fara_C(1:current_index_J_Fara_C-1),index2_J_Fara_C(1:current_index_J_Fara_C-1),element_J_Fara_C(1:current_index_J_Fara_C-1),N*eqn_ma,N*eqn_ma);
            else
                J_Fara_C = 0;
                b_Fara_C = 0;
            end
        end
        %% 12. Calculate delta and update the variables

        J_total_w = J+J_Bulk_Rea+J_Bulk_Rea_mi+J_chem_Rea+p_mi./dt.*Matrix_temp1+J_Fara_A+J_Fara_C;

        b_total = b+b_Bulk_Rea+b_Bulk_Rea_mi+b_chem_Rea-p_mi./dt.*(Vector_temp1+b3)+b_Fara_A+b_Fara_C;
        
        delta = J_total_w\b_total;
 
        v = v + delta;    % v is the v_star, updated every iteration
        delta2 = J2_sum\b2_sum+J2_sum\J3*delta;

        v2 = v2+delta2;
        

%         figure(1)
%         clf
%         hold on
%         summation = 0;
%         for i = 1:s_mo
%             summation = summation + z_mo_input(i).*v(i:eqn_ma:N*eqn_ma);
%         %             plot(x_c,v(i:eqn_ma:N*eqn_ma),'-o')
%         end
%         plot(x_c,summation,'-o')
%         hold off
%         text(3e-3,10,sprintf('t=%f',t))
%         %         xlim([0.75*l_e l_e+l_sep+0.25*l_e])
%         xlabel('x')
%         ylabel('c')
%         drawnow
%             
%         figure(2)
%         clf
%         hold on
%         for i = 1:s_mo
%             plot(x_c,v2(i:eqn_mi:N*eqn_mi),'-o')
%         end    
% %         for i = eqn_ma+1:eqn_mi
% %             plot(x_c,v2(i:eqn_mi:N*eqn_mi),'-o')
% %         end    
% %         plot(x_c,v2(index_HX:eqn_mi:N*eqn_mi)+v2(index_X:eqn_mi:N*eqn_mi),'-o')
% %         plot(x_c,v2(index_HY:eqn_mi:N*eqn_mi)+v2(index_Y:eqn_mi:N*eqn_mi),'-o')
% 
%         hold off
%         xlabel('x')
%         ylabel('c_{mi}')
%         drawnow
% % % 
%         figure(3)
%         clf
%         hold on
%         plot(x_c,v(eqn_ma:eqn_ma:N*eqn_ma),'-o')
%         plot(x_c,v2(eqn_ma:eqn_mi:N*eqn_mi)*F/V_T/C_S,'-o')
%         plot(x_c,phi_el-v(eqn_ma:eqn_ma:N*eqn_ma)-v2(eqn_ma:eqn_mi:N*eqn_mi)*F/V_T/C_S,'-o')
%         plot(x_c,phi_el,'-o')
%         hold off
%         xlabel('x')
%         ylabel('\phi')
%         legend('\phi','\phi_{St}','\phi_{Donnan}','\phi_{electrode}')
%         drawnow
% % 
%     figure(4)
%     clf
%     hold on
%     plot(x_c,v2(eqn_ma:eqn_mi:N*eqn_mi),'-o')
%     hold off
%     xlabel('x')
%     ylabel('\sigma_{ele}')
%     drawnow

%     figure(7)
%     clf
%     hold on
%     for i = 3
%         plot(x_c,-log10(v2(i:eqn_mi:N*eqn_mi)/1e3),'-o')
%     end
%     plot(x_c,-log10(v2(3:eqn_mi:N*eqn_mi)/1e3) -log10(v2(4:eqn_mi:N*eqn_mi)/1e3),'-o')
%     plot(x_c,-log10(v2(3:eqn_mi:N*eqn_mi).*v2(7:eqn_mi:N*eqn_mi)./v2(6:eqn_mi:N*eqn_mi)/1e3),'-o')
%     plot(x_c,-log10(v2(3:eqn_mi:N*eqn_mi).*v2(9:eqn_mi:N*eqn_mi)./v2(8:eqn_mi:N*eqn_mi)/1e3),'-o')
% 
%     hold off
%     text(3e-3,14,sprintf('t=%.2f',t))
%     legend('pH','pH+pOH','pKa5','pKa9')
% %     ylim([12 16])
%     xlabel('x')
%     ylabel('pH & pKa')
%         drawnow


    end

    negative_value_check = 0;
    for i = 1:s_mo
        negative_value_check = sum(v(i:eqn_ma:N*eqn_ma)<0)+negative_value_check;
    end
    for i = 1:s_ch
        negative_value_check = sum(v2(eqn_ma+i:eqn_mi:N*eqn_mi)<0)+negative_value_check;
    end

    complex_value_check = ~isreal(v2) + ~isreal(v2) + ~isreal(v) + ~isreal(v);
    nan_value_check = sum(isnan(v2)) + sum(isnan(v2)) + sum(isnan(v)) + sum(isnan(v));

    disp(dt)
    if negative_value_check == 0 && complex_value_check == 0 && nan_value_check == 0
        dt_count = dt_count+1;
        t = t+dt;
        if t == t_final
            save_function(0,initial_condition_solver,v,v2,t,phi_el_value,save_name,save_all_name,save_path);
            break;
        end
        if t_plan(saved_data_index) <= t
            saved_data_index = saved_data_index+1; % for plot y vs time
            if mod(saved_data_index,1) == 0
                save_function(0,initial_condition_solver,v,v2,t,phi_el_value,save_name,save_all_name,save_path);
            end
        end
        if dt_stored_switch == 1
            dt = dt_stored;
            dt_stored_switch = 0;
        end
        if dt_count > 9 && dt < dt_max
            dt = dt*2;
            if dt > dt_max
                dt = dt_max;
            end
            dt_count = 0;
        end
        if steady_state_solver == 0
            % Reassign initial time step for unsteady state solver
            if mod(t,FCT) == 0 || mod(t,FCT) == Switch_Time
                dt = dt_min;
                dt_count = 0;
                continue
            end
        end
        % Make sure the saved time is integer
        if (t_plan(saved_data_index)-(t+dt)) < 1e-10 || t+dt > t_plan(saved_data_index)
            dt_stored = dt;
            dt_stored_switch = 1;
            dt = t_plan(saved_data_index)-t;
        end
    else
        dt_count = 0;
        dt = dt/2;
        v = v_n;
        v2 = v2_n;
    end
end
end

function save_function(save_all_index,initial_condition_solver,v,v2,t,phi_el_value,save_name,save_all_name,save_path)
    if save_all_index == 1
        save([save_path,save_all_name])
    else
        if initial_condition_solver == 1
            t = 0;
            save([save_path,save_name],'v','v2','t','phi_el_value')
        else
            save([save_path,save_name,'_t_',num2str(t)],'v','v2','t','phi_el_value')
        end
    end
end

function mesh=optimal_mesh_gen(alpha,N_eta,boundary_left,boundary_right,direction)
    %     alpha = 1;  % Decide the density, the larger the more extreme near the edge
    %     N_eta = 10; % How many elements (Mesh Center)
        if direction == 1 % Left denser
            eta = linspace(-alpha,0,N_eta+1)';
        elseif direction == 2 % Two sides denser
            eta = linspace(-alpha,alpha,N_eta+1)';
        elseif direction == 3 % Right side denser
            eta = linspace(0,alpha,N_eta+1)';
        end
        y_value = tanh(eta);
        mesh = boundary_left+(tanh(eta)-y_value(1))*(boundary_right-boundary_left)/(y_value(end)-y_value(1));
    end

function X=nonuniform_generator(length,boundary,N,x_f)
    index_X = 1;
    X(index_X) = x_f(end);
    dx = length/N;
    while 1
        index_X = index_X+1;
        X(index_X,1) = X(index_X-1,1)+dx*randi(10,1)*0.2;

        if X(index_X) >= boundary
            X(index_X) = boundary;
            break
        end
    end
    X = X(2:end);
end

function NewMatrix=FullMatrix_Extract(FullMatrix,N,eqn_ma,eqn_mi)
    for n = 1:N
        NewMatrix((n-1)*eqn_ma+1:(n-1)*eqn_ma+eqn_ma,:)=FullMatrix((n-1)*eqn_mi+1:(n-1)*eqn_mi+eqn_ma,:);
    end
end

function current_mag = current_calculator(current_collect_index,eqn_ma,s_mo,z_mo_input,D_eff_f,v,delta_x_f,F,A_c)
    flux_D_pt = 0;
    flux_E_pt = 0;
    current_index_local = current_collect_index;
    current_index_global = (current_index_local-1)*eqn_ma;% zero index
    for i = 1:s_mo
    
        flux_D_pt = flux_D_pt-z_mo_input(i)*D_eff_f(current_index_global+i).*(v(current_index_global+i+eqn_ma)-v(current_index_global+i))./delta_x_f(current_index_local);
        flux_E_pt = flux_E_pt-z_mo_input(i)*D_eff_f(current_index_global+i)*z_mo_input(i).*(v(current_index_global+i+eqn_ma)+v(current_index_global+i))/2.*(v(current_index_global+s_mo+1+eqn_ma)-v(current_index_global+s_mo+1))./delta_x_f(current_index_local);
    
    end
    %         figure(999)
    %         plot(x_f(2:end-1)',(flux_D+flux_E)*F*A_c)
    current_mag = (flux_D_pt+flux_E_pt)*F*A_c;
end
% Oxygen 8mg/L*(1/32)(mol/g)=0.00025 mol/L= 0.25mM
