function force = periodic_force_at_time(time_s, Omega, force_samples)

n = size(force_samples,2);
phase_index = mod(Omega*time_s,2*pi)/(2*pi)*n;
left = floor(phase_index)+1;
weight = phase_index-floor(phase_index);
right = left+1;
if right > n, right = 1; end
force = (1-weight)*force_samples(:,left)+weight*force_samples(:,right);
end
