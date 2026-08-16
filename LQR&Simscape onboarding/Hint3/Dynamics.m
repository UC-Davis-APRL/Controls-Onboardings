function x_dot = Dynamics(x, F)

    % Extract states
    pos   = x(1);
    vel   = x(2);
    theta = x(3);
    omega = x(4);

    % Parameters
    M = 1;
    m = 1;
    l = 1;
    r = 0.1;
    lc = l/2;

    % Moment of inertia of the cylindrical rod about its center of mass
    I_com = (1/12) * m * (3*r^2 + l^2);

    % Moment of inertia about the pivot
    I = I_com + m*lc^2;

    g = 9.81;

    A = [
        M + m,             m*lc*cos(theta)
        m*lc*cos(theta),   I
    ];

    b = [
        F + m*lc*sin(theta)*omega^2
        m*g*lc*sin(theta)
    ];

    % Solve for accelerations
    accelerations = A\b;

    x_ddot = accelerations(1);
    theta_ddot = accelerations(2);

    % Return state derivative
    x_dot = [
        vel
        x_ddot
        omega
        theta_ddot
    ];

end