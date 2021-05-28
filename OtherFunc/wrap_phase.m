function phase_out = wrap_phase(phase_in);
% wrap the phase in [-pi pi]
% Usage: phase_out=wrap_phase(phase_in);

phase_out=mod(phase_in+pi,2*pi)-pi;
