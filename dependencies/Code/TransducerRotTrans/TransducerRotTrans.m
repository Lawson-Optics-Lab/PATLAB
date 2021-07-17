function [transducer_pos_new] = TransducerRotTrans( eff_x,eff_y,eff_z,eff_phi,eff_theta )
%%version 1
%Takes the log file coordinates (one at a time, so loop if you want this
%for multiple scan points), calculates the centre of rotation, and applies
%the appropriate rotations to match the transducer positions happening in
%the scan.

R = 254.25; %Radius of the robot rotation, this could be made into a variable in the future
%%this section calculates the centre of rotation
cor_x = eff_x - R.*sind(eff_theta).*cosd(eff_phi);
cor_y = eff_y - R.*sind(eff_theta).*sind(eff_phi);
cor_z = eff_z - R.*cosd(eff_theta);
centre_rot = [cor_x cor_y cor_z]; % Centre of rotation

%%This applies the transforms
load('20200311_Transducer_Position_Home.mat','transducer_pos_home','centre_rot_home'); % Load the mat file containing the calibrated transducer locations, but rotated back to phi and theta = 0
delta_cor = centre_rot - centre_rot_home; %calculate diff between centre of rotation of current scan data and home position
transducer_pos = transducer_pos_home; 
transducer_pos = transducer_pos + delta_cor; %add difference to each scan point
transducer_pos = transducer_pos';% transpose the array to match that needed by AxelRot
u = [0 0 1]; %set phi rotation axis
[transducer_desired_afterphi, ~, ~] = AxelRot(transducer_pos, eff_phi, u, centre_rot); %Rotate points for phi, going through centre of rotation
[u2x,u2y,u2z] = sph2cart(deg2rad(eff_phi-90),0,1); % Create vector that is 90 deg from the current position of the array to set as rotation axis
[transducer_desired_aftertheta, ~, ~] = AxelRot(transducer_desired_afterphi, -eff_theta, [u2x u2y u2z], centre_rot); %Rotate points for theta, going through centre of rotation
transducer_pos_new = transducer_desired_aftertheta';%transpose back to original orientation