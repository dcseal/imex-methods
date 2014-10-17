% Split Integrator.  This routine takes a classical RK4 step, following by a
% step of trapezoid (Crank-Nicholson) followed by a half step of RK4.
function ynp1 = split_integrator( tn, dt, yn, de )

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initial step
    dt = 0.5*dt;
    k1 = fE( tn, yn );
    k2 = fE( tn + 0.5*dt, yn + 0.5*dt*k1 );
    k3 = fE( tn + 0.5*dt, yn + 0.5*dt*k2 );
    k4 = fE( tn + 1.0*dt, yn + 1.0*dt*k3 );
    yn = yn + dt/6*(k1 + 2*(k2+k3) + k4 );

    % trapezoidal step
    dt = 2*dt;
    ys = yn + 0.5 * dt * fI( tn, yn );
    yn = implicit_solve( ys, tn, 0.5*dt );

    % final step
    dt = 0.5*dt;
    tn = tn + dt;
    k1 = fE( tn, yn );
    k2 = fE( tn + 0.5*dt, yn + 0.5*dt*k1 );
    k3 = fE( tn + 0.5*dt, yn + 0.5*dt*k2 );
    k4 = fE( tn + 1.0*dt, yn + 1.0*dt*k3 );
    yn = yn + dt/6*(k1 + 2*(k2+k3) + k4 );


    ynp1 = yn;


end
