clear 
clc
syms pos vel theta omega real

syms F_in real

syms M m l r lc g real

    % Moment of inertia of the cylindrical rod about its center of mass
    I_com = (1/12) * m * (3*r^2 + l^2);

    % Moment of inertia about the pivot
    I = I_com + m*lc^2;

    I_mat = [
        M + m,             -m*lc*cos(theta)
        -m*lc*cos(theta),   I
    ];

    F_vec = [
        F_in + m*lc*sin(theta)*omega^2
        m*g*lc*sin(theta)
    ];

    % Solve for accelerations
    A_vec = I_mat\F_vec;

    x_ddot = A_vec(1);
    theta_ddot = A_vec(2);

    x = [pos; vel; theta; omega]
   
    x_dot = [
        vel
        x_ddot
        omega
        theta_ddot
    ];

    A_lin = jacobian(x_dot, x);
    B_lin = jacobian(x_dot, F_in);

    A_lin = simplify(A_lin);
    B_lin = simplify(B_lin);

    A_num = subs(A_lin, ...
        [M m l r lc g], ...
        [1 1 1 0.1 0.5 9.81]);

    B_num = subs(B_lin, ...
        [M m l r lc g], ...
        [1 1 1 0.1 0.5 9.81]);

    A_num = simplify(A_num);
    B_num = simplify(B_num);

    A0 = double(subs(A_num, ...
        [pos vel theta omega F_in], ...
        [0 0 0 0 0]))
    B0 = double(subs(B_num, ...
        [pos vel theta omega F_in], ...
        [0 0 0 0 0]))