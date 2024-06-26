% exMD is compatible with MATLAB and GNU Octave (www.octave.org).     
% The script plots the step response of a 2nd-order mass-damper-spring
% system.
%
% Author:    Thor I. Fossen
% Date:      16th June 2001
% Revisions: 

wn = 1;

figure(gcf)
subplot(211)
t = 0:0.01:20;
z = 0.5; sys = tf(wn*wn,[1 2*z*wn wn*wn]); step(sys,t)
hold on
z = 1.0; sys = tf(wn*wn,[1 2*z*wn wn*wn]); step(sys,t)
z = 2.0; sys = tf(wn*wn,[1 2*z*wn wn*wn]); step(sys,t)
hold off
legend('\zeta = 0.5','\zeta = 1.0','\zeta = 2.0')

subplot(212)
t = 0:0.01:50;
z = 0.1; sys = tf(wn*wn,[1 2*z*wn wn*wn]); step(sys,t)
hold on
sys = tf(wn*wn,[1 0 wn*wn]); step(sys,t)
hold off
legend('\zeta = 0.1','\zeta = 0.0')

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
