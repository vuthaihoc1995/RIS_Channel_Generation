clc; clear variables; close all; 
dbclear all
L=16; % Number of RIS element
sample=1e4; % Number of channel generation
snrdB = -10:30;
%%  Pathloss
L_SR = 0.5; % From source to RIS
L_SB = 0.1; % From  source to Bob
L_RB = 0.15; % From RIS to Bob
 L_SW = 1; % From RIS to Willie 
 L_RW = 1; % From source to Bob
% Nakagami-m factor
mu_h_Ray = @(n,lambda) gamma(n/2 + 1)*(lambda)^(n/2);

% Nakagami-m factor
m_SR=2;m_RB=3;m_SB =3;m_RW=3;m_SW =1;
mu_h_Naka = @(n,m,lambda) gamma(n/2 + m)/gamma(m)*(lambda/m)^(n/2);
% Rician factor
k_SR=2;k_RB=3;k_SB =1;k_RW=3;k_SW =0;
mu_h_Rician = @(n,k,lambda) exp(-k)*gamma(n/2 +1 )*( lambda/(k+1))^(n/2)*hypergeom(n/2+1,1,k);
%% Complex Channel Generation
parfor ss=1:sample
       
        % Rayleigh enviroment
         h_Ray_SR(ss,:)  = sqrt(L_SR)*sqrt(1/2)*(randn(L,1) + 1i*randn(L,1)); 
         h_Ray_RW(ss,:) = sqrt(L_RW)*sqrt(1/2)*(randn(L,1) + 1i*randn(L,1)); 
         h_Ray_SW(ss)   = sqrt(L_SW)*sqrt(1/2)*(randn(1,1) + 1i*randn(1,1)); 
         h_Ray_RB(ss,:)  = sqrt(L_RB)*sqrt(1/2)*(randn(L,1) + 1i*randn(L,1)); 
         h_Ray_SB(ss)    = sqrt(L_SB)*sqrt(1/2)*(randn(1,1) + 1i*randn(1,1)); 

        % Nakagami-m enviroment
        h_Naka_SR(ss,:)  = sqrt(L_SR).* sqrt(gamrnd(m_SR,1/m_SR,L,1)).* exp(1i*2*pi*rand(L,1));
        h_Naka_RW(ss,:)= sqrt(L_RW).*sqrt(gamrnd(m_RW,1/m_RW,L,1)).* exp(1i*2*pi*rand(L,1));
        h_Naka_SW(ss)   = sqrt(L_SW).*sqrt(gamrnd(m_SW,1/m_SW,1,1)).* exp(1i*2*pi*rand(1,1));
        h_Naka_RB(ss,:)  = sqrt(L_RB).*sqrt(gamrnd(m_RB,1/m_RB,L,1)).* exp(1i*2*pi*rand(L,1));
        h_Naka_SB(ss)    = sqrt(L_SB).* sqrt(gamrnd(m_SB,1/m_SB,1,1)).* exp(1i*2*pi*rand(1,1));
   

        % % % Rician enviroment
         mean_SR =  sqrt(k_SR/(k_SR+1)); varSR =  sqrt(1/(k_SR+1))*sqrt(1/2);
        h_Rician_SR(ss,:)   = sqrt(L_SR).*random('Rician',mean_SR,varSR,L,1).* exp(1i*2*pi*rand(L,1));

        mean_RW =  sqrt(k_RW/(k_RW+1)); varRW =  sqrt(1/(k_RW+1))*sqrt(1/2);
        h_Rician_RW(ss,:)  = sqrt(L_RW).*random('Rician',mean_RW,varRW,L,1).* exp(1i*2*pi*rand(L,1));

        mean_SW =  sqrt(k_SW/(k_SW+1)); varSW =  sqrt(1/(k_SW+1))*sqrt(1/2);
        h_Rician_SW(ss)    = sqrt(L_SW).*random('Rician',mean_SW,varSW,1,1).* exp(1i*2*pi*rand(1,1));

        mean_RB =  sqrt(k_RB/(k_RB+1)); varRB =  sqrt(1/(k_RB+1))*sqrt(1/2);
        h_Rician_RB(ss,:)   = sqrt(L_RB).*random('Rician',mean_RB,varRB,L,1).* exp(1i*2*pi*rand(L,1));

        mean_SB =  sqrt(k_SB/(k_SB+1)); varSB =  sqrt(1/(k_SB+1))*sqrt(1/2);
        h_Rician_SB(ss)     = sqrt(L_SB).*random('Rician',mean_SB,varSB,1,1).* exp(1i*2*pi*rand(1,1));


end

%% Casecaded Complex Channel Generation at Willie (random phase-shift) + Bob (optimal phase-shift)
OnOff = 0; % OnOff = 0/1 - A factor for link from source to Willie/Bob
for ss=1:sample
    %% Rayleigh
     % Phase-shift + cascaded channel generation at Bob
        for ll=1:L
            phase_shift_config_Ray =  angle(h_Ray_SB(ss)) - angle(h_Ray_SR(ss,ll))   - angle(h_Ray_RB(ss,ll)).';
            phase_shift_complex_vector_Ray(ll) = exp(1i*phase_shift_config_Ray);
        end
        phase_shift_matrix_RIS_Ray = diag(phase_shift_complex_vector_Ray(1:L));
 
        ch_h_Ray_SB(ss) =  OnOff*h_Ray_SB(ss) + h_Ray_SR(ss,:)*phase_shift_matrix_RIS_Ray*h_Ray_RB(ss,:).';
       % cascaded channel generation at Willie 
        ch_h_Ray_SW(ss)     =  OnOff*h_Ray_SW(ss) + h_Ray_SR(ss,:)*phase_shift_matrix_RIS_Ray*h_Ray_RW(ss,:).';

        %% Nakagami-m
        % Phase-shift + cascaded channel generation at Bob
        for ll=1:L
            phase_shift_config_Naka =  angle(h_Naka_SB(ss)) - angle(h_Naka_SR(ss,ll))   -angle(h_Naka_RB(ss,ll)).';
            phase_shift_complex_vector_Naka(ll) = exp(1i*phase_shift_config_Naka);
        end
        phase_shift_matrix_RIS_Naka = diag(phase_shift_complex_vector_Naka(1:L));
 
        ch_h_Naka_SB(ss) =  OnOff*h_Naka_SB(ss) + h_Naka_SR(ss,:)*phase_shift_matrix_RIS_Naka*h_Naka_RB(ss,:).';
       % Phase-shift + cascaded channel generation at Willie
        ch_h_Naka_SW(ss)   =  OnOff*h_Naka_SW(ss) + h_Naka_SR(ss,:)*phase_shift_matrix_RIS_Naka*h_Naka_RW(ss,:).';

    %% Rician
    % Phase-shift + cascaded channel generation at Willie
    ch_h_Rician_SW(ss) =  OnOff*h_Rician_SW(ss) + h_Rician_SR(ss,:)*phase_shift_matrix_RIS_W*h_Rician_RW(:,ss);

    % Phase-shift + cascaded channel generation at Bob
        for ll=1:L
            phase_shift_config_Rician =  angle(h_Rician_SB(ss)) - angle(h_Rician_SR(ss,ll))   - angle(h_Rician_RB(ss,ll)).';
            phase_shift_complex_vector_Rician(ll) = exp(1i*phase_shift_config_Rician);
        end
        phase_shift_matrix_RIS_Rician = diag(phase_shift_complex_vector_Rician(1:L));
        ch_h_Rician_SB(ss) =  OnOff*h_Rician_SB(ss) + h_Rician_SR(ss,:)*phase_shift_matrix_RIS_Rician*h_Rician_RB(ss,:).';

        % Phase-shift + casecaded channel generation at Willie
        ch_h_Rician_SW(ss) =  OnOff*h_Rician_SW(ss) + h_Rician_SR(ss,:)*phase_shift_matrix_RIS_Rician*h_Rician_RW(ss,:).';

end


for idx=1:length(snrdB)

   snr = db2pow(snrdB(idx)); % Transmit SNR
   gamma_th =2^3-1;            % Threshold SNR
   %% At Willie
   % Rayleigh
   % Simulation OP 
   OP_sim_Ray_W(idx) = mean( abs(ch_h_Ray_SW).^2*snr < gamma_th);
   % Analysis OP
   beta_Ray_W = 1/( OnOff*mu_h_Ray(2,L_SW)   + L*mu_h_Ray(2,L_SR)*mu_h_Ray(2,L_RW)  );
   OP_ana_Ray_W(idx) = 1 - exp(- beta_Ray_W*gamma_th/snr );
  
   % Naka
   % Simulation OP   
   OP_sim_Naka_W(idx) = mean( abs(ch_h_Naka_SW).^2*snr < gamma_th);
   % Analysis OP
   beta_Naka_W = 1/(OnOff*mu_h_Naka(2,m_SW,L_SW)   + L*mu_h_Naka(2,m_SR,L_SR)*mu_h_Naka(2,m_RW,L_RW));
   OP_ana_Naka_W(idx) = 1 - exp(-beta_Naka_W*gamma_th/snr );

   % Rician
   % Simulation OP 
   OP_sim_Rician_W(idx) = mean( abs(ch_h_Rician_SW).^2*snr < gamma_th);
   % Analysis OP
   beta_Rician_W = 1/(OnOff*mu_h_Rician(2,k_SW,L_SW)   + L*mu_h_Rician(2,k_SR,L_SR)*mu_h_Rician(2,k_RW,L_RW));
   OP_ana_Rician_W(idx) = 1 - exp(-beta_Rician_W*gamma_th/snr );

%% At Bob
% Rayleigh
% Simulation OP
OP_sim_Ray_B(idx) = mean( abs(ch_h_Ray_SB).^2*snr < gamma_th); 
% Analysis OP
alpha_B = ( OnOff*mu_h_Ray(1,L_SB)   + L*mu_h_Ray(1,L_SR)*mu_h_Ray(1,L_RB) )^2/( OnOff*(mu_h_Ray(2,L_SB) -  mu_h_Ray(1,L_SB)^2)   + L*(mu_h_Ray(2,L_SR)*mu_h_Ray(2,L_RB)  - mu_h_Ray(1,L_SR)^2*mu_h_Ray(1,L_RB)^2 )  );
beta_B   = ( OnOff*mu_h_Ray(1,L_SB)   + L*mu_h_Ray(1,L_SR)*mu_h_Ray(1,L_RB) )   /(  OnOff*(mu_h_Ray(2,L_SB) -  mu_h_Ray(1,L_SB)^2)    + L*(mu_h_Ray(2,L_SR)*mu_h_Ray(2,L_RB)  - mu_h_Ray(1,L_SR)^2*mu_h_Ray(1,L_RB)^2 )  );
OP_ana_Ray_B(idx) = 1 - igamma(alpha_B,beta_B*sqrt(gamma_th/snr) )/gamma(alpha_B);

% Naka
% Simulation OP
OP_sim_Naka_B(idx) = mean( abs(ch_h_Naka_SB).^2*snr < gamma_th);
% Analysis OP
alpha_Naka_B = ( OnOff*mu_h_Naka(1,m_SB,L_SB)   + L*mu_h_Naka(1,m_SR,L_SR)*mu_h_Naka(1,m_RB,L_RB) )^2/...
    ( OnOff*(mu_h_Naka(2,m_SB,L_SB) -  mu_h_Naka(1,m_SB,L_SB)^2)   + L*(mu_h_Naka(2,m_SR,L_SR)*mu_h_Naka(2,m_RB,L_RB)  - mu_h_Naka(1,m_SR,L_SR)^2*mu_h_Naka(1,m_RB,L_RB)^2 )  );
beta_Naka_B   = ( OnOff*mu_h_Naka(1,m_SB,L_SB)   + L*mu_h_Naka(1,m_SR,L_SR)*mu_h_Naka(1,m_RB,L_RB) )/...
    ( OnOff*(mu_h_Naka(2,m_SB,L_SB) -  mu_h_Naka(1,m_SB,L_SB)^2)   + L*(mu_h_Naka(2,m_SR,L_SR)*mu_h_Naka(2,m_RB,L_RB)  - mu_h_Naka(1,m_SR,L_SR)^2*mu_h_Naka(1,m_RB,L_RB)^2 )  );
OP_ana_Naka_B(idx) = 1 - igamma(alpha_Naka_B,beta_Naka_B*sqrt(gamma_th/snr) )/gamma(alpha_Naka_B);

% Rician
% Simulation OP
OP_sim_Rician_B(idx) = mean( abs(ch_h_Rician_SB).^2*snr < gamma_th);
% Analysis OP
alpha_Rician_B = ( OnOff*mu_h_Rician(1,k_SB,L_SB)   + L*mu_h_Rician(1,k_SR,L_SR)*mu_h_Rician(1,k_RB,L_RB) )^2/...
    ( OnOff*(mu_h_Rician(2,k_SB,L_SB) -  mu_h_Rician(1,k_SB,L_SB)^2)   + L*(mu_h_Rician(2,k_SR,L_SR)*mu_h_Rician(2,k_RB,L_RB)  - mu_h_Rician(1,k_SR,L_SR)^2*mu_h_Rician(1,k_RB,L_RB)^2 )  );
beta_Rician_B   = ( OnOff*mu_h_Rician(1,k_SB,L_SB)   + L*mu_h_Rician(1,k_SR,L_SR)*mu_h_Rician(1,k_RB,L_RB) )/...
    ( OnOff*(mu_h_Rician(2,k_SB,L_SB) -  mu_h_Rician(1,k_SB,L_SB)^2)   + L*(mu_h_Rician(2,k_SR,L_SR)*mu_h_Rician(2,k_RB,L_RB)  - mu_h_Rician(1,k_SR,L_SR)^2*mu_h_Rician(1,k_RB,L_RB)^2 )  );
OP_ana_Rician_B(idx) = 1 - igamma(alpha_Rician_B,beta_Rician_B*sqrt(gamma_th/snr) )/gamma(alpha_Rician_B);

end

%% Validation Between Simulation and Analytical
% Rayleigh
figure(1)
semilogy(snrdB, OP_sim_Ray_W, 'ro'); hold on;
semilogy(snrdB, OP_ana_Ray_W, 'k-x'); hold on;

semilogy(snrdB, OP_sim_Ray_B, 'bs'); hold on;
semilogy(snrdB, OP_ana_Ray_B, 'k-+'); hold on;
lgd=legend(...
    'Willie (sim.)',...
    'Willie (ana.)',...
    'Bob (sim.)', ...
     'Bob (ana.)', ...
    'Location','se',...
    'Interpreter', 'Latex');
lgd.NumColumns = 1;
lgd.FontSize = 14;
set(gca,'fontsize',14);
title('Rayleigh Fading')
xlabel('Transmit SNR');
ylabel('Outage Probability (OP)');
axis([-Inf Inf 1e-4 10^(0)]);
% Naka
figure(2)
semilogy(snrdB, OP_sim_Naka_W, 'ro'); hold on;
semilogy(snrdB, OP_ana_Naka_W, 'k-x'); hold on;

semilogy(snrdB, OP_sim_Naka_B, 'bs'); hold on;
semilogy(snrdB, OP_ana_Naka_B, 'k-+'); hold on;
lgd=legend(...
    'Willie (sim.)',...
    'Willie (ana.)',...
    'Bob (sim.)', ...
     'Bob (ana.)', ...
    'Location','se',...
    'Interpreter', 'Latex');
lgd.NumColumns = 1;
lgd.FontSize = 14;
set(gca,'fontsize',14);
title('Nakagami-m Fading')
xlabel('Transmit SNR');
ylabel('Outage Probability (OP)');
axis([-Inf Inf 1e-4 10^(0)]);
% Rician
figure(3)
 semilogy(snrdB, OP_sim_Rician_W, 'ro'); hold on;
semilogy(snrdB, OP_ana_Rician_W, ' k-x'); hold on;

 semilogy(snrdB, OP_sim_Rician_B, 'bs'); hold on;
semilogy(snrdB, OP_ana_Rician_B, 'k-+'); hold on;
lgd=legend(...
    'Willie (sim.)',...
    'Willie (ana.)',...
    'Bob (sim.)', ...
     'Bob (ana.)', ...
    'Location','se',...
    'Interpreter', 'Latex');
lgd.NumColumns = 1;
lgd.FontSize = 14;
set(gca,'fontsize',14);
title('Rician Fading')
xlabel('Transmit SNR');
ylabel('Outage Probability (OP)');
axis([-Inf Inf 1e-4 10^(0)]);




