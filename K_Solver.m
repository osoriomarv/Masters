%Def K val and solve using previously established methodologies

% Data definitions 
T= linspace(900, 1000, 10);
T_ref=298;
P=1;
R=8.314;
P_ref=1;
%Excel sheet with shomate equation coefficients
data=readmatrix('SiO2FeOH2O.xlsx','Sheet','Sheet1','Range','B2:H10');
A=data(1,:);
B=data(2,:);
C=data(3,:);
D=data(4,:);
E=data(5,:);
F=data(6,:);
G=data(7,:);
H_Kj=data(8,:); %Kj/mol
H_jmol=data(8,:).*1000; %J/mol
S_298=data(9,:); %J/mol/k

%Pressure
        %FeO    Fe     SiO2
Molar_v=[1.2040 0.7092 2.2608]; %J/bar

%Number of Species (MUST BE THE SAME VALUE)
num_p = numel(P);
num_T = numel(T);
num_a = numel(A);

%Prellocation
VdP_FeO = zeros(1,num_T);
VdP_SiO2 = zeros(1,num_T);
VdP_H2O = zeros(1,num_T);
t_inv=zeros(1, num_T);
Cp_T=zeros(num_T,num_a);
dH_species=zeros(num_T,num_a);
S_Species=zeros(num_T,num_a);
delHrxnlogFeO=zeros(1,num_T);
delSrxnlogFeO=zeros(1,num_T);
Grxn_FeO_Jmol=zeros(1,num_T);
log_K_FeO=zeros(1,num_T);
delHrxnlogSiO2=zeros(1,num_T);
delSrxnlogSiO2=zeros(1,num_T);
Grxn_SiO2_Jmol=zeros(1,num_T);
log_K_SiO2=zeros(1,num_T);
delHrxnlogH2O=zeros(1,num_T);
delSrxnlogH2O=zeros(1,num_T);
Grxn_H2O_Jmol=zeros(1,num_T);
log_K_H2O=zeros(1,num_T);
Grxn_FeO_KJmol=zeros(1,num_T);
Grxn_SiO2_KJmol=zeros(1,num_T);
Grxn_H2O_KJmol=zeros(1,num_T);
%K_val=zeros(4,10)

for i=1:num_T
for p=1:num_p   
    %Fe + .5O2 = FeO
    VdP_FeO(i) = Molar_v(1).*P(p) - (Molar_v(2).*P(p));
    %SiO2 = SiO + .5O2
    VdP_SiO2(i) = (.5*(R.*T(i).*log(P(p)./P_ref)) + (R.*T(i).*log(P(p)./P_ref)))  - (Molar_v(3).*P(p)); %J/mol
    %H2 + .5O2 = H2O
    VdP_H2O(i) = (R.*T(i).*log(P(p)./P_ref)) - (.5.*(R.*T(i).*log(P(p)./P_ref))  + R.*T(i).*log(P(p)./P_ref)); %J/mol
    %SiO2 + H2 = H2O + SiO
    VdP_H2O_SiO2(i) = R.*T(i).*log(P(p)./P_ref) + (Molar_v(3).*P(p)); %J/mol

    %Inverse Temperature
t_inv(i)=T(i)/1000;
    %Heat Capacity
Cp_T(i,:)= A + B*t_inv(i) + C*t_inv(i).^2 + D.*t_inv(i).^3 + E/t_inv(i).^2;

%Enthalpy and Entropy
dH_species(i,:) = (A*t_inv(i)+ (B*t_inv(i).^2)/2 + (C*t_inv(i).^3)/3 + (D*t_inv(i).^4)/4 - E./t_inv(i)+ F - H_Kj).*1000 %J/mol
S_Species(i,:) = (A*log(t_inv(i)) + B*t_inv(i)+ C*t_inv(i).^2/2 + D*t_inv(i).^3/3 - E./(2*t_inv(i).^2) + G) %J/mol/k

%Fe + .5O2 = FeO
Hrxn_29815_FeO = (H_jmol(7) - H_jmol(6) - .5.*H_jmol(2)); %J/mol
Srxn_29815_FeO = (S_298(7) - S_298(6) - .5.*S_298(2));%J/mol/K
Grxn_29815_FeO = Hrxn_29815_FeO - 298.15.*(Srxn_29815_FeO);%J/mol
K_29815_feO = (-Grxn_29815_FeO./(2.2303.*R.*T_ref));
%Temperature and Pressure corrections
Hrxn_FeO = (Hrxn_29815_FeO + (dH_species(i,7) -dH_species(i,6) - .5*dH_species(i,2)));%j/mol
delHrxnlogFeO(i)=Hrxn_FeO; %Logs the Hrxn FeO calculations
Srxn_FeO = S_Species(i,7) - S_Species(i,6) - .5*S_Species(i,2); %J/mol
delSrxnlogFeO(i)=Srxn_FeO; %Logs the entropy 
Grxn_FeO_Jmol(i)= (Hrxn_FeO - T(i).*(Srxn_FeO) + VdP_FeO(i));% J/mol
log_K_FeO(i) = (-Grxn_FeO_Jmol(i)./(2.303.*R.*T(i)));

% Iron_Wustite Buffer Calculation
Grxn_FeO_Jmol_IW(i)= (Hrxn_FeO - T(i).*(Srxn_FeO));

IW_FO2(i)=(exp((Grxn_FeO_Jmol_IW(i))./(-R.*T(i)))).^-2;

IW_FO2_Picked = IW_FO2(i);

%SiO2 = SiO + .5O2
Hrxn_29815_SiO2 = (.5.*H_jmol(2) + H_jmol(1)) - H_jmol(3); %J/mol
Srxn_29815_SiO2 = (.5.*S_298(2)+ S_298(1)) - S_298(3); %J/Mol
Grxn_29815_SiO2 = Hrxn_29815_SiO2 - 298.15.*(Srxn_29815_SiO2);%J/mol
K_29815_SiO2 = (-Grxn_29815_SiO2./(2.2303.*R.*T_ref));

%Correction for Temperature and Pressure
Hrxn_SiO2 = (Hrxn_29815_SiO2 + ((.5.*dH_species(i,2) + dH_species(i,1)) - dH_species(i,3)));%j/mol
delHrxnlogSiO2(i)=Hrxn_SiO2; %logs the enthalpy
Srxn_SiO2 = (.5.*S_Species(i,2) + S_Species(i,1)) - S_Species(i,3); %J/mol
delSrxnlogSiO2(i)=Srxn_SiO2; %Logs the entropy
Grxn_SiO2_Jmol(i)= (Hrxn_SiO2 - T(i).*(Srxn_SiO2) + VdP_SiO2(i));%J/mol
%Grxn_SiO2_Jmol(i)= (Hrxn_SiO2 - T(i).*(Srxn_SiO2));%J/mol
log_K_SiO2(i) = (-Grxn_SiO2_Jmol(i)./(2.303.*R.*T(i)));

%H2 + .5O2 = H2O
Hrxn_29815_H2O =  (H_jmol(4) - (.5.*H_jmol(2) + H_jmol(5))); %J/mol
Srxn_29815_H2O =  S_298(4) - (.5.*S_298(2)+ S_298(5)); %J/Mol
Grxn_29815_H2O = Hrxn_29815_H2O - 298.15.*(Srxn_29815_H2O);%J/mol
K_29815_H2O = (-Grxn_29815_H2O./(2.2303.*R.*T_ref));

%Correction for Temperature and Pressure
Hrxn_H2O = Hrxn_29815_H2O + (dH_species(i,4) - dH_species(i,5) - 0.5.*dH_species(i,2));%J/mol
delHrxnlogH2O(i)=Hrxn_H2O;
Srxn_H2O =  S_Species(i,4) - S_Species(i,5) - 0.5.*S_Species(i,2);
delSrxnlogH2O(i)=Srxn_H2O;
Grxn_H2O_Jmol(i)= (Hrxn_H2O - T(i).*(Srxn_H2O) + VdP_H2O(i)); %kj/mol
log_K_H2O(i) = (-Grxn_H2O_Jmol(i)./(2.303.*R.*T(i)));

% SiO2 + H2 = H2O + SiO
Hrxn_29815_H2O_SiO2 =  H_jmol(4) + H_jmol(1) - H_jmol(5) - H_jmol(3); %J/mol
Srxn_29815_H2O_SiO2 =  S_298(4) + S_298(1) - S_298(5) - S_298(3); %J/Mol
Grxn_29815_H2O_SiO2 = Hrxn_29815_H2O_SiO2 - 298.15.*(Srxn_29815_H2O_SiO2);%J/mol
K_29815_H2O = (-Grxn_29815_H2O_SiO2./(2.2303.*R.*T_ref));

%Correction for Temperature and Pressure
Hrxn_H2O_SiO2 = Hrxn_29815_H2O_SiO2 + (dH_species(i,4) + dH_species(i,1) - dH_species(i,5) - dH_species(i,3));%J/mol
delHrxnlogH2O_SiO2(i)=Hrxn_H2O_SiO2;
Srxn_H2O_SiO2 =  S_Species(i,4) + S_Species(i,1) - S_Species(i,5) - S_Species(i,3);
delSrxnlogH2O_SiO2(i)=Srxn_H2O_SiO2;
Grxn_H2O_SiO2_Jmol(i)= (Hrxn_H2O_SiO2 - T(i).*(Srxn_H2O_SiO2) + VdP_H2O_SiO2(i)); %kj/mol
log_K_H2O_SiO2(i) = (-Grxn_H2O_SiO2_Jmol(i)./(2.303.*R.*T(i)));

%Gibbs Free energies in 
Grxn_FeO_KJmol(i)=Grxn_FeO_Jmol(i)./1000; %J/mol
Grxn_SiO2_KJmol(i)=Grxn_SiO2_Jmol(i)./1000; %J/mol
Grxn_H2O_KJmol(i)=Grxn_H2O_Jmol(i)./1000; %J/mol

%Def K Val and passing to main function
log_K_SiO2_picked = log_K_SiO2(i);
log_K_FeO_picked = log_K_FeO(i);
log_K_H2O_picked = log_K_H2O(i);
log_K_H2O_SiO2_picked = log_K_H2O_SiO2(i);


%K_val = [log_K_SiO2(:) log_K_FeO(:) log_K_H2O(:) IW_FO2(:)];
%K_val = [log_K_SiO2; log_K_FeO; log_K_H2O; IW_FO2];

end
end