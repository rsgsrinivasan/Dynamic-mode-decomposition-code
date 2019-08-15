% SREE ANJINEYA
% OM GAM GANAPATIYE NAMAHA
%% Dynamic mode Decomposition and DMD reconstruction code
% Developed by G Srinivasan, Priyam Chakraborty and Prof K.P.Sinhamapatra
% Indian Institute of Technology Kharagpur India

%% Flushing memory, deleting old output files and creating necessary folders 
clc;
clear all;
system('rm *.plt');
system('rm Modes/*.plt');
system('rm Recon_field/*.plt');
%% INPUT FILES AND SNAPSHOTS DETAILS
input_files_path = '../M219_DMD_';
coordinates_file = '../COORDINATES_FILE.plt';
connectivity = '../connectivity.txt';
NUMNP = 3378000;
NELEM = 3280256;
file_start = 200000;
Last_file = 1100000;
delta = 1000;
Num_files = (Last_file - file_start)/delta;       %dim_time is total time in seconds
dim_time = Last_file*0.005;                       %delt_given is delta_t of files in seconds
delt_given = 0.005;
dt = 1000*0.005;


U_Process = 1;
V_Process = 0;
W_Process = 0;
UVW_Process = 0;
Tecplot_mode_out = 11;                             % Tecplot output files to write
NUM_MODES = (Last_file - file_start)/delta;        % Takes all modes
Phase_shift_Recon = 0;								% Enter 1 for phase angle reconstruction
Recon_Modes = 12;
mode_user_def = 0;
% 
tic
tsart = tic;
tn =1;
for temp = 0:Num_files
i = file_start + delta*temp;
filename = sprintf('%s%d.plt',input_files_path,i);
A = importdata(filename); 
  if UVW_Process == 1
      C1 = A(:,2:4);
      C1 = C1';
      V(:,tn) = C1(:); 
  end
  if U_Process == 1
      Pressure = A(:,1); 
	  S_Entropy = A(:,4); 
	  
	  V(:,tn) = A(:,1); 
  end
  if V_Process == 1
      V(:,tn) = A(:,1); 
  end    
if W_Process == 1
      V(:,tn) = A(:,1); 
  end 
tn =tn+1;
fprintf('file %d, %s\n',tn-1,filename);
end
clear A C1;
time_taken = toc(tsart);
fprintf('file read time_taken = g mins\n\n',time_taken/60.0);

Vmean = sum(V,2)/size(V,2);

%% SNAPSHOTS PREPARATION
N = size(V,2);    
V1 = V(:,1:N-1);
V2 = V(:,2:N);
tn = Num_files+1;
r = NUM_MODES ;        % COMPUTED_MODES 

%% Compute full matrix S and pullout dynamic modes and spectrum
clear U E W S X D dm Phi;
[U,E,W] = svd(V1,'econ');
U_tr = U(:,1:NUM_MODES);
E_tr = E(1:NUM_MODES, 1:NUM_MODES);
W_tr = W(:,1:NUM_MODES);
clear U E W;
S = conj(U_tr)' * V2 * W_tr * pinv(E_tr);
[X,D] = eig(S);
clear S;
lambda = diag(D);
lam = log(diag(D))/dt;
for j = 1:1:size(V1,2)
	for i = 1:1:size(V1,2)
		Vand(i,j) = (lambda(i))^(j-1);
	end
end
dm = V2 * W_tr * pinv(E_tr) * X;                  % Exact DMD Tu et al 
%dm = U_tr*X;                                       % Projected DMD Schmid 2010
Mode_Ampl = Vand*(W_tr*(inv(E_tr))*X);                                      
Mode_Amp = (diag(inv(Mode_Ampl)));

%dm = dm*inv(Mode_Ampl);
%% Compute DMD mode amplitudes for modes identification and reconstruction
energy = (abs((Mode_Amp)')).^2; 
AMP = abs((inv((conj(X)'*X).*(conj(Vand*conj(Vand)')))*(conj(diag(Vand*W_tr*conj(E_tr)'*X))))');
clear X U_tr E_tr W_tr Mode_Ampl Vand

%% OUTPUT FILES
% REAL_vs_IMAG_SPECTRUM.plt    ==> MODE SPECTRUM PLOT (LAMBDA IMAG vs REAL)
% REAL_vs_IMAG_UNIT_CIRCLE.plt ==> LAMBDA REAL vs IMAG ON UNIT CIRCLE PLOT
% FREQUENCY_vs_AMPLITUDE.plt   ==> FREQUENCY vs L1 NORM PLOT
% FREQUENCY_vs_ENERGY.plt      ==> FREQUENCY vs GLOBAL ENERGY NORM PLOT (ENERGY --> INNER PRODUCT OF MODES)

fid = fopen('REAL_vs_IMAG_SPECTRUM.plt', 'w');
%fid=fopen('.plt','w');
for i=1:r
fprintf(fid,'%g %g\n',imag(lam(i))/(2.0*pi),real(lam(i))/(2.0*pi));
end
fclose(fid);
fid = fopen('REAL_vs_IMAG_UNIT_CIRCLE.plt', 'w');
hh = diag(D);
%fid=fopen('.plt','w');
for i=1:r
fprintf(fid,'%g %g\n',real(hh(i,1)),imag(hh(i,1)));
end
fclose(fid);
fid = fopen('FREQUENCY_vs_AMPLITUDE.plt', 'w');
%fid=fopen('.plt','w');
for i=1:r
%fprintf(fid,'%g %g\n',imag(lam(i)/(2.0*pi)),abs(real(energy(i))) );
fprintf(fid,'%g %g %g\n',imag(lam(i)/(2.0*pi)),AMP(i), abs(lambda(i)) );
end
fclose(fid);
fid = fopen('FREQUENCY_vs_ENERGY.plt', 'w');
for i=1:r
fprintf(fid,'%g\t%g\t%g\n',imag(lam(i)/(2.0*pi)),abs(Mode_Amp(i)), abs(Mode_Amp(i))^2 );
end
fclose(fid);

%% MATLAB FIGURES OUTPUT
figure (1)
plot(imag(lam(1:r)),real(lam(1:r)),'b.','markersize',20)
%set(gca,'DataAspectRatio',[1 25 1])
xlabel(['\lambda','_i'])
ylabel(['\lambda','_r'])
grid on

figure(2)
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(D)),imag(diag(D)),'ok')
axis([-1.1 1.1 -1.1 1.1]);

%% MODE NUMBERS ARRANGEMENT IN ASCENDING ORDER (RUN THIS IF NOT ABLE TO FIND MODES VISUALLY)

[TEMP_VAL, IND] = sort(AMP,'descend');
IND(1,size(IND,2)+1) = IND(1,size(IND,2)-1);
clear mode Phi omega DMDb time_dynamics Xdmd;
k =1;
mm = 1;
for i = 1:1:size(IND,2)-1
	if abs(lambda(IND(1,i))) > 0.995 && imag(lam(IND(1,i))) ~= 0 && (lam(IND(1,i)) == conj(lam(IND(1,i+1))) || lam(IND(1,i)) == conj(lam(IND(1,i-1))))  
		mode(k,1) =  IND(1,i);
		k = k+1;
	end
	if imag(lam(IND(1,i))) == 0
		Mean_mode(mm,1) =  IND(1,i);
		mm = mm+1;
    end
end
Fluc_modes = k-1;
fid = fopen('MODES_EIGENVALUES_LIST.plt', 'w');
for i=1:size(mode,1)
fprintf(fid,'%d\t%d\t%g\t%g\n', i, mode(i,1), imag(lam(mode(i,1))/(2.0*pi)),imag(lam(mode(i,1))/(2.0*pi)) );
end
fclose(fid);

TEMP_MODE = mode;
clear mode;
for i =1:1:Recon_Modes;
	mode(i,1) = TEMP_MODE(i,1);
end
i = i+1;
for j =1:1:1;
	mode(i,1) = Mean_mode(j,1);
	i = i+1;
end
fid = fopen('MODES_EIGENVALUES_LIST2.plt', 'w');
for i=1:size(mode,1)
fprintf(fid,'%d\t%d\t%g\t%g\n', i, mode(i,1), imag(lam(mode(i,1))/(2.0*pi)),imag(lam(mode(i,1))/(2.0*pi)) );
end
fclose(fid);

if mode_user_def == 1
    clear mode;
     for i = 1:1:Recon_Modes
         mode(i,1) = input('mode numbers ');
     end
     if Phase_shift_Recon == 1 && mode_user_def == 1
         Freq = input('mode frequency ');
         Phase_angles = input('Enter number of  equally spaced snapshots to reconstruct ');
         mode_dt = 1.0/Freq;
     end
end
Recon_Modes= size(mode,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time varying amplitude
Admd = zeros(size(V1,1),size(V1,2));
for i = 1:1:Recon_Modes
    Phi(:,i) = dm(:,mode(i,1));
    omega(i,1) = lam(mode(i,1));
	rec_lambda(i,1) = lambda(mode(i,1));
end
D_ji = pinv(Phi)*V1;
for i = 1:1:Recon_Modes
	Admd = Admd + Phi(:,i)*D_ji(i,:);
end
mm1 = size(V1, 2);
time_dynamics = zeros(Recon_Modes, mm1);
t = (0:mm1-1)*(dt); % time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DMD reconstruction for single mode
% MODES LIST SEGREGATION AND MODE AMPLITUDE CALC FOR RECONSTRUCTION
if Phase_shift_Recon == 1 && mode_user_def == 1
	for i = 1:1:Recon_Modes
		Phi(:,i) = dm(:,mode(i,1));
		omega(i,1) = lam(mode(i,1));
		rec_lambda(i,1) = lambda(mode(i,1));
		%DMDb(i,1) = AMP(mode(i,1));
	end
	x1 = V1(:, 1);
	DMDb = pinv(Phi)*x1;
	mm1 = size(V1, 2); % mm1 = m - 1
	if Phase_shift_Recon == 1 && mode_user_def == 1
	   mm1 =  Phase_angles+1;
	end
	time_dynamics = zeros(Recon_Modes, mm1);
	t = (0:mm1-1)*(dt); % time vector
	%t = t+file_start*delt_given;
	if Phase_shift_Recon == 1 && mode_user_def == 1
	   for i = 0:1:Phase_angles-1
		   Phase(i+1,1) = i*(2*pi/Phase_angles);
	   end
	   Phase(i+2,1) = 2*pi;
	   t(1) = 0.0; 
	   for i = 2:1:size(Phase,1)+1
		   t(i) = (t(i-1)+1.0);
	   end
	   t = t*abs(atan(imag(rec_lambda(1,1))/real(rec_lambda(1,1))));
	end
	%t = t.*dt;
	for iter = 1:1:mm1
	  time_dynamics(:,iter) = (DMDb.*(exp(omega*t(iter))));
	  if Phase_shift_Recon == 1 && mode_user_def == 1
		  time_dynamics(:,iter) = (DMDb.*(exp(sqrt(-1)*t(iter))));
	  end
    end	 
	Admd = Phi * time_dynamics;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
if UVW_Process == 1
    for i = 1:1:size(Admd,2)
        TEMP_VAR = reshape(Admd(:,i),[3,(size(Admd,1)/3)])';
        Xdmd(:,i) = TEMP_VAR(:,1);
        Ydmd(:,i) = TEMP_VAR(:,2);
		Ydmd(:,i) = TEMP_VAR(:,3);

        TEMP_VAR = reshape(V(:,i),[3,(size(V,1)/3)])';
        U1(:,i) = TEMP_VAR(:,1);
        U2(:,i) = TEMP_VAR(:,2);
		U3(:,i) = TEMP_VAR(:,3);
    end
    for i = 1:1:Recon_Modes
        TEMP_VAR = reshape(Phi(:,i),[3,(size(Phi,1)/3)])';
        UPhi(:,i) = TEMP_VAR(:,1);
        VPhi(:,i) = TEMP_VAR(:,2);
		WPhi(:,i) = TEMP_VAR(:,3);
    end
end
if U_Process == 1 || V_Process == 1 || W_Process == 1
   Xdmd = Admd; 
end
clear Admd;

%% Tecplot output files for modes visualisation
A = importdata(coordinates_file);
XYZ(:,1) = A(:,1);
XYZ(:,2) = A(:,2);
XYZ(:,3) = A(:,3);
CONNECT = importdata(connectivity);
for MODES = 1:1:Tecplot_mode_out+2
	filename =  sprintf('Modes/MODE%d_%d.plt',MODES,mode(MODES));
	fp= fopen(filename,'w');
	fprintf(fp,'TITLE = \"Node file\"\n');
	if UVW_Process == 1
        fprintf(fp, 'VARIABLES = \"X\", \"Y\", \"Z\", \"U_REAL_LAMDA\", \"U_IMAG_LAMDA\",\"V_REAL_LAMDA\", \"V_IMAG_LAMDA\", \"W_REAL_LAMDA\", \"W_IMAG_LAMDA\",\n');
    end
    if U_Process == 1
        fprintf(fp, 'VARIABLES = \"X\", \"Y\", \"Z\", \"U_REAL_LAMDA\", \"U_IMAG_LAMDA\",\n');
    end
    if V_Process == 1
        fprintf(fp, 'VARIABLES = \"X\", \"Y\", \"Z\", \"V_REAL_LAMDA\", \"V_IMAG_LAMDA\",\n');
    end
	if W_Process == 1
        fprintf(fp, 'VARIABLES = \"X\", \"Y\", \"Z\", \"W_REAL_LAMDA\", \"W_IMAG_LAMDA\",\n');
    end
	fprintf(fp, 'ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n',NUMNP,NELEM);
	for i=1:NUMNP
		if UVW_Process == 1
			fprintf(fp,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n', XYZ(i,1), XYZ(i,2), XYZ(i,3), real(UPhi(i,MODES)), imag(UPhi(i,MODES)), real(VPhi(i,MODES)), imag(VPhi(i,MODES)), real(WPhi(i,MODES)), imag(WPhi(i,MODES)) );
		end
		if U_Process == 1 || V_Process == 1 || W_Process
			fprintf(fp,'%g\t%g\t%g\t%g\t%g\n', XYZ(i,1), XYZ(i,2), XYZ(i,3), real(Phi(i,MODES)), imag(Phi(i,MODES)));
		end
	end
	fprintf(fp,'\n\n\n');
	for i=1:NELEM
		fprintf(fp,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', CONNECT(i,1), CONNECT(i,2), CONNECT(i,3), CONNECT(i,4), CONNECT(i,5), CONNECT(i,6), CONNECT(i,7), CONNECT(i,8) );
	end		
	fclose(fp);
end

%% Tecplot output files of reconstructed field
if mode_user_def == 1
	directory = sprintf('mkdir Recon_field/Mode_%d',mode(1,1));
	system(directory);
end
if Phase_shift_Recon ~= 1
	directory = sprintf('mkdir Recon_field/Recon_Multiple_Modes');
	system(directory);
	system('rm Recon_field/Recon_Multiple_Modes/*');
end
A = importdata(coordinates_file);
XYZ(:,1) = A(:,1);
XYZ(:,2) = A(:,2);
XYZ(:,3) = A(:,3);
CONNECT = importdata(connectivity);
if UVW_Process == 1
     TEMP_VAR = reshape(Vmean(:,1),[3,(size(Admd,1)/3)])';
     c(:,4) = TEMP_VAR(:,1);
     c(:,5) = TEMP_VAR(:,2);
     c(:,6) = TEMP_VAR(:,3);
end
if U_Process == 1 || V_Process == 1 || W_Process == 1
    c(:,4) = Vmean; 
end

if Phase_shift_Recon ~= 1
    hop =1;
end
if Phase_shift_Recon == 1
    hop = (pi/4)*(1/abs(atan(imag(rec_lambda(1,1))/real(rec_lambda(1,1)))));
    hop = floor(hop);
end

for TIMESTEP = 1:hop:size(time_dynamics,2)
	if Phase_shift_Recon ~= 1
		filename =  sprintf('Recon_field/Recon_Multiple_Modes/Recon_%d.plt',TIMESTEP);
	end
	if mode_user_def == 1 && Phase_shift_Recon == 1
		filename =  sprintf('Recon_field/Mode_%d/Recon_%d.plt',mode(1,1),TIMESTEP);
	end
	fp= fopen(filename,'w');
	fprintf(fp,'TITLE = \"Node file\"\n');
    if UVW_Process == 1
        fprintf(fp, 'VARIABLES = \"X\", \"Y\", \"Z\",\"Actual_U\", \"Actual_V\", \"Actual_W\", \"Reconstructed_U\", \"Reconstructed_V\", \"Reconstructed_W\", \"MEAN_U\", \"MEAN_V\", \"MEAN_W\",\n');
    end
    if U_Process == 1 
        fprintf(fp, 'VARIABLES = \"X\", \"Y\", \"Z\",\"Instantaneous_U\", \"Reconstructed_U\", \"MEAN_U\",\n');
    end
    if V_Process == 1
        fprintf(fp, 'VARIABLES = \"X\", \"Y\", \"Z\",\"Instantaneous_V\", \"Reconstructed_V\",  \"MEAN_V\",\n');
    end
	if W_Process == 1
        fprintf(fp, 'VARIABLES = \"X\", \"Y\", \"Z\",\"Instantaneous_W\", \"Reconstructed_W\", \"MEAN_W\", \n');
    end
	if Phase_shift_Recon ~= 1
		fprintf(fp, 'ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK, SOLUTIONTIME=%g\n',NUMNP,NELEM,TIMESTEP*dt);
	end
	if Phase_shift_Recon == 1 && mode_user_def == 1
		fprintf(fp, 'ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK, SOLUTIONTIME=%g\n',NUMNP,NELEM,TIMESTEP);
	end
	for i=1:NUMNP
		if UVW_Process == 1
			fprintf(fp,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n', XYZ(i,1), XYZ(i,2), XYZ(i,3), U1(i,TIMESTEP), U2(i,TIMESTEP), U3(i,TIMESTEP), real(Xdmd(i,TIMESTEP)), real(Ydmd(i,TIMESTEP)), real(Zdmd(i,TIMESTEP)), UMean(i), VMean(i), WMean(i));
		end
		if U_Process == 1 || V_Process == 1 || W_Process
			fprintf(fp,'%g\t%g\t%g\t%g\t%g\t%g\n', XYZ(i,1), XYZ(i,2), XYZ(i,3), V(i,TIMESTEP), real(Xdmd(i,TIMESTEP)), Vmean(i));
		end	
	end
	fprintf(fp,'\n\n\n');
	for i=1:NELEM
		fprintf(fp,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', CONNECT(i,1), CONNECT(i,2), CONNECT(i,3), CONNECT(i,4), CONNECT(i,5), CONNECT(i,6), CONNECT(i,7), CONNECT(i,8) );
	end		
	fclose(fp);
end


