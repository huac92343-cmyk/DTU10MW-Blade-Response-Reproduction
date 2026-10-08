function history = houbolt_integrate(M, C, K, force_function, time, q0, v0)

arguments
    M double
    C double
    K double
    force_function (1,1) function_handle
    time double
    q0 double
    v0 double
end

time = time(:).';
if numel(time) < 4 || any(diff(time) <= 0)
    error('DTU10MW:TimeGrid', 'Houbolt integration needs at least four increasing times.');
end
dt_all = diff(time);
dt = dt_all(1);
if max(abs(dt_all-dt)) > 100*eps(max(time))
    error('DTU10MW:TimeGrid', 'Houbolt integration requires a uniform time step.');
end
n = numel(time);
ndof = numel(q0);
q = zeros(ndof,n);
v = zeros(ndof,n);
a = zeros(ndof,n);
q(:,1) = q0(:);
v(:,1) = v0(:);
a(:,1) = M\(force_function(time(1))-C*v(:,1)-K*q(:,1));

beta = 1/4;
gamma = 1/2;
K_newmark = K+gamma/(beta*dt)*C+1/(beta*dt^2)*M;
solver_newmark = decomposition(K_newmark,'lu');
for j = 1:2
    rhs = force_function(time(j+1))+ ...
        M*(q(:,j)/(beta*dt^2)+v(:,j)/(beta*dt)+ ...
        (1/(2*beta)-1)*a(:,j))+ ...
        C*(gamma*q(:,j)/(beta*dt)+(gamma/beta-1)*v(:,j)+ ...
        dt*(gamma/(2*beta)-1)*a(:,j));
    q(:,j+1) = solver_newmark\rhs;
    a(:,j+1) = (q(:,j+1)-q(:,j))/(beta*dt^2)- ...
        v(:,j)/(beta*dt)-(1/(2*beta)-1)*a(:,j);
    v(:,j+1) = v(:,j)+dt*((1-gamma)*a(:,j)+gamma*a(:,j+1));
end

K_houbolt = K+11/(6*dt)*C+2/dt^2*M;
solver_houbolt = decomposition(K_houbolt,'lu');
for j = 3:n-1
    rhs = force_function(time(j+1))+ ...
        C*(3*q(:,j)/dt-3*q(:,j-1)/(2*dt)+q(:,j-2)/(3*dt))+ ...
        M*(5*q(:,j)/dt^2-4*q(:,j-1)/dt^2+q(:,j-2)/dt^2);
    q(:,j+1) = solver_houbolt\rhs;
    v(:,j+1) = (11*q(:,j+1)-18*q(:,j)+9*q(:,j-1)-2*q(:,j-2))/(6*dt);
    a(:,j+1) = (2*q(:,j+1)-5*q(:,j)+4*q(:,j-1)-q(:,j-2))/dt^2;
end
history.time_s = time;
history.displacement = q;
history.velocity = v;
history.acceleration = a;
history.dt_s = dt;
history.startup_method = 'Newmark average acceleration for first two steps';
history.main_method = 'third-order Houbolt backward difference';
end
