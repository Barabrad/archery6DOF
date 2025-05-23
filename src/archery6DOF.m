function diff_coords = archery6DOF(L_shaft, OR_shaft, IR_shaft, ...
    GPI_shaft, W_nock, W_tip, ang_tip, L_tip_in, L_fin, h_fin, W_fin, ...
    N_fins, ang_fin, v_arrow, h_arrow, ang_aim_lr, ang_aim_ud, ...
    d_target, h_target, v_wind, dir_wind, rho_air)
% This function will run a 6DOF simulation of the firing of an arrow.
% Inputs:
% - L_shaft: length of arrow shaft [in]
%   - Must be positive
% - OR_shaft: outer radius of arrow shaft [in]
%   - Must be positive
% - IR_shaft: inner radius of arrow shaft [in]
%   - Must be positive (or 0 for a solid arrow)
% - GPI_shaft: grains per inch of arrow shaft [GPI]
%   - 1 grain = 1/7000 lbf. Grains per inch is the industry standard
%   - Must be positive
% - W_nock: weight of arrow nock [grain]
%   - Must be nonnegative
% - W_tip: weight of arrow tip [grain]
%   - Must be positive
% - ang_tip: angle made by a cross-section of the tip [deg]
%   - Range: 15º to 45º
% - L_tip_in: length of the arrow tip's insert [in]
%   - Range: 0 in to L_shaft
% - L_fin: length of a single tail fin along the shaft [in]
%   - Must be nonnegative
% - h_fin: maximum height of a single tail fin above the shaft [in]
%   - Must be nonnegative
% - ang_fin: off-shaft angle of a single tail fin [deg]
%   - 0º: fins are perfectly aligned with shaft -> no spin during flight
%   - Range: -10º to 10º (positive angle will cause the arrow's spin vector
%     to be in the same direction as its pointing vector)
% - W_fin: weight of a single fin [grain]
%   - Must be nonnegative
% - N_fins: number of equally-spaced fins [-]
%   - Range: 3 to 4
% - v_arrow: speed of the arrow upon leaving the bow [ft/s]
%   - There are a lot of parameters needed to calculate this from bow
%     specifications, so I might add them in later
%   - Range: 100 ft/s to 500 ft/s
% - h_arrow: height of the arrow's tail upon leaving the bow [ft]
%   - Must be high enough such that the tip of the arrow isn't touching
%     the ground (h = 0)
% - ang_aim_lr: left-right aiming angle [deg]
%   - Angle between the aimed arrow's pointing vector and the plane
%     perpendicular to the ground that contains the archer and target
%   - Positive angles go right
% - ang_aim_ud: up-down aiming angle [deg]
%   - Angle between the aimed arrow's pointing vector and the ground plane
%   - Positive angles go up
% - d_target: distance of the target [ft]
%   - Should be positive, but I'll allow negative for trick shots that take
%     advantage of strong winds
% - h_target: height of the target's bulls-eye [ft]
%   - Must be positive
% - v_wind: wind speed [ft/s]
%   - Must be nonnegative
%   - Default value is 0 ft/s
% - dir_wind: a 2D column vector that gives the direction of wind [-]
%   - Let y-hat be the vector pointing from the archer to the target
%   - Let x-hat be perpendicular to y-hat and pointing to the right of the
%     archer
%     - This makes z-hat = cross(x-hat, y-hat), which isn't included in
%       this 2D vector
%   - Doesn't have to be a unit vector
%   - Default vector is [0;0]
% - rho_air: air density [lbm/ft^3] (optional)
%   - I've seen bow specs measured with Imperial Units, so as much as it
%     annoys me, I won't use SI for consistency's sake
%   - Must be positive
%   - Default value is 0.0752 lbm/ft^3 (1.204 kg/m^3)
% Outputs:
% - diff_coords: a 3D column vector giving the differences in coordinates
%   of the arrow relative to the target when either 1) the arrow crosses
%   the plane of the target (containing x-hat and z-hat) OR 2) the arrow
%   hits the ground [ft]
%   - The coordinate frame of the vector will be based on the three hat
%     vectors defined earlier, so getting a negative y value means that you
%     undershot the target
% Assumptions:
% - Constant air density, wind speed, and gravity over the arrow's flight
% - The archer and target are in an inertial frame
% - Arrow is a streamlined body for the purposes of drag
%   (https://www.engineeringtoolbox.com/drag-coefficient-d_627.html)
% - Moments from drag on a component can be modeled by having the drag
%   force act on the component's center of mass
% - The archer won't impart any rotation on the arrow when firing, such
%   that the initial velocity vector matches the arrow's pointing vector
% - Fins are right triangles with negligible thickness with their bases at
%   at the end of the arrow shaft
% - The nock's center of mass is located at the end of the shaft (since
%   there is an inserted part)
% - The nock and tip are close enough to axisymmetric
% - The nock's length is negligible for purposes of calculating drag

    %% Account for missing arguments and validate inputs
    if (nargin < 22)
        rho_air = 0.0752; % lbm/ft^3
    else
        checkIfPositive(rho_air);
    end
    if (nargin < 21)
        dir_wind = [0;0];
    else
        % Make sure it's a normalized 2D column vector
        dir_wind = makeUnitVec([dir_wind(1); dir_wind(2)]);
    end
    if (nargin < 20)
        v_wind = 0; % ft/s
    else
        checkIfPositive(v_wind, true);
    end
    checkIfPositive(h_target);
    checkIfInRange(ang_fin, -10, 10);
    checkIfInRange(N_fins, 3, 4);
    checkIfInRange(v_arrow, 100, 500);
    checkIfPositive(v_wind, true);
    checkIfPositive(L_fin, true);
    checkIfPositive(h_fin, true);
    checkIfPositive(L_shaft);
    checkIfInRange(L_tip_in, 0, L_shaft);
    checkIfInRange(ang_tip, 15, 45);
    checkIfPositive(W_tip);
    checkIfPositive(W_nock, true);
    checkIfPositive(IR_shaft, true);
    checkIfPositive(OR_shaft);
    checkIfPositive(GPI_shaft);

    %% Globals
    ft2in = 12;
    lbf2grain = 7000;
    g_imp = 32.174; % ft/s^2

    %% Find moments of inertia
    [pI_vec, z_arrow_COM, z_tip_COM, m_shaft, m_nock, m_tip, m_fin, ...
        L_arrow] = getArrowMOI(L_shaft, OR_shaft, IR_shaft, GPI_shaft, ...
        W_nock, W_tip, ang_tip, L_tip_in, L_fin, h_fin, W_fin, N_fins, ...
        g_imp, lbf2grain);
    m_arrow = m_shaft + m_nock + m_tip + (m_fin*N_fins);
    [I1, I2, I3] = struct('x', num2cell(pI_vec)).x;
    fprintf("[I1, I2, I3] = [%f, %f, %f] lb•in^2\n", I1, I2, I3);
    fprintf("COM is %f inches from the back of the arrow\n", z_arrow_COM);
    fprintf("The arrow's mass is %f lbm\n", m_arrow);
    

    %% Set up ODE
    % Convert lengths to ft
    [I1, I2, I3] = struct('x', num2cell(pI_vec/(ft2in^2))).x;
    z_arrow_COM = z_arrow_COM/ft2in;
    z_tip_COM = z_tip_COM/ft2in;
    L_arrow = L_arrow/ft2in;
    L_shaft = L_shaft/ft2in;
    L_tip = L_arrow - L_shaft;
    L_fin = L_fin/ft2in;
    h_fin = h_fin/ft2in;
    OR_shaft = OR_shaft/ft2in;
    L_com2tip = (L_arrow - z_arrow_COM);
    % Find the inital state (x, v, and rot)
    aim_dir = [sind(ang_aim_lr)*cosd(ang_aim_ud); ...
               cosd(ang_aim_lr)*cosd(ang_aim_ud); ...
               sind(ang_aim_ud)];
    com2tip0_F = L_com2tip*aim_dir;
    com2tip0_C = [0; 0; L_com2tip];
    tip0 = [0; 0; h_arrow] + L_arrow*aim_dir;
    % Check if tip is above the ground
    if (tip0(3) <= 0)
        error("Arrow tip must be above the ground at time of firing");
    end
    x0 = tip0 - com2tip0_F;
    v0 = v_arrow*aim_dir;
    w0 = [0;0;0]; % No spin initially
    % We want the initial quaternion to be the rotation from the arrow's
    % shaft axis to the inertial frame orientation
    % (These quaternions can also convert body-fixed angular momentum to
    % inertial angular momentum)
    % We can use principal-axis rotation matrices to describe the rotation
    % from the inertial frame to the body-fixed frame: z(-90º), y(-90º),
    % x(-aim_lr), then y(aim_ud)
    F2C = Rz(-90) * Ry(-90) * Rx(-ang_aim_lr) * Ry(ang_aim_ud);
    e0 = rotMat2Quat(F2C);
    y0 = [x0; v0; w0; e0];
    % Set ODE tolerances
    [rel_tol, abs_tol] = deal(1e-8, 1e-10);
    fprintf("Relative Tolerance: %e\n", rel_tol);
    fprintf("Absolute Tolerance: %e\n", abs_tol);
    options = odeset('RelTol',rel_tol, 'AbsTol',abs_tol, 'Events',...
        @(t,y) eventsFun(t, y, L_com2tip, d_target));
    % Set time bounds
    % Estimate the time it will take for the arrow to land using kinematics
    % 0 = h_arrow + v0_z*t - 0.5*g*t^2
    % t = -v0_z +/- sqrt(v0_z^2 + 2*h_arrow*g), but t >= 0
    v0_z = v0(3);
    t_max = (sqrt(v0_z^2 + 2*h_arrow*g_imp) - v0_z)*2; % Margin for drag
    t_bnds = [0; t_max]; % t_max can overshoot since events will stop ode45
    fprintf("Time bound: %f s\n", t_max);
    % Solve ODEs
    % Make anonymous function that deals with other parameters
    odefun_p = @(t,y) odefun(t, y, I1, I2, I3, dir_wind, v_wind, ...
        ang_fin, L_fin, h_fin, N_fins, z_arrow_COM, z_tip_COM, ...
        OR_shaft, L_shaft, L_tip, rho_air, g_imp, m_tip, m_shaft, ...
        m_fin, m_nock);
    [t_vec, y_vecs, t_events] = ode45(odefun_p, t_bnds, y0, options);

    if (isempty(t_events))
        warning("Arrow was still in flight when alloted time was " + ...
            "reached. Check your function inputs.")
    else
        fprintf("Simulation ended at t = %f s\n", t_vec(end));
    end

    % Return coordinate difference
    target_coords = [0; d_target; h_target];
    state_final = y_vecs(end,:)';
    com_coords_final = state_final(1:3);
    e_final = state_final(10:13);
    com2tip1 = rotVecWithQuat(com2tip0_C, e_final);
    tip_coords_final = com_coords_final + com2tip1;
    diff_coords = tip_coords_final - target_coords;

    % Find tip location
    [x_vecs, v_vecs, w_vecs, e_vecs] = deal(y_vecs(:,1:3), ...
        y_vecs(:,4:6), y_vecs(:,7:9), y_vecs(:,10:13));
    tip_vecs = zeros(size(x_vecs));
    for i = 1:length(t_vec)
        x_i = x_vecs(i,:)';
        e_i = e_vecs(i,:)';
        tip_vecs(i,:) = (x_i + rotVecWithQuat(com2tip0_C, e_i))';
    end
    % Plot
    leg_text = ["Arrow COM"; "Arrow Tip"];
    fig1 = graph3D_wComps(x_vecs, t_vec, "x", "ft", "s", true);
    graph3D_wComps(tip_vecs, t_vec, "x", "ft", "s", false, fig1, leg_text);
    graph3D_wComps(v_vecs, t_vec, "v", "ft/s");
    graph3D_wComps(w_vecs, t_vec, "\omega", "rad/s");
    graphFinalState(diff_coords, e_final, L_arrow);
    
end


%% ODE functions

function [value, isterminal, direction] = eventsFun(~, y, L_com2tip, ...
    y_target)
    [value, isterminal, direction] = deal([0;0]); % Preallocate
    % y = [x; v; w; e]
    % Find current location of tip
    e1 = y(10:13);
    com2tip1 = rotVecWithQuat([0;0;L_com2tip], e1);
    y_tip = y(2) + com2tip1(2);
    z_tip = y(3) + com2tip1(3);
    % Configure event for crossing target plane
    value(1) = y_tip - y_target; % The value that we want to be zero
    isterminal(1) = 1;  % Halt integration
    direction(1) = 0;   % The zero can be approached from either direction
    % Configure event for hitting ground
    value(2) = z_tip; % The value that we want to be zero
    isterminal(2) = 1;  % Halt integration
    direction(2) = 0;   % The zero can be approached from either direction
end


function dydt = odefun(~, y, I1, I2, I3, dir_wind, v_wind, ang_fin, ...
    L_fin, h_fin, N_fins, z_arrow_COM, z_tip_COM, OR_shaft, L_shaft, ...
    L_tip, rho_air, g, m_tip, m_shaft, m_fin, m_nock)
    % y has 13 elements: 3 for x, 3 for v, 3 for wdot, and 4 for edot
    dydt = zeros(size(y));
    [v, w, e] = deal(y(4:6), y(7:9), y(10:13));
    % xdot = v
    dydt(1) = v(1);
    dydt(2) = v(2);
    dydt(3) = v(3);
    % vdot = a, which includes drag and gravity
    [arrow_DF, arrow_MC] = findDragFM(e, w, v, dir_wind, v_wind, ...
        ang_fin, L_fin, h_fin, N_fins, z_arrow_COM, OR_shaft, L_shaft, ...
        L_tip, rho_air);
    dydt(4) = arrow_DF(1);
    dydt(5) = arrow_DF(2);
    dydt(6) = arrow_DF(3) - g;
    % wdot requires us to solve for torques
    g_MC = findGravityM(e, g, z_arrow_COM, z_tip_COM, L_shaft, L_fin, ...
        N_fins, m_tip, m_shaft, m_fin, m_nock);
    M = arrow_MC + g_MC;
    % Break up the Euler equations:
    % M1 = I1*wdot1 + w2*w3*(I3-I2)
    % M2 = I2*wdot2 + w1*w3*(I1-I3)
    % M3 = I3*wdot3 + w1*w2*(I2-I1)
    dydt(7) = (M(1) + w(2)*w(3)*(I2-I3))/I1;
    dydt(8) = (M(2) + w(1)*w(3)*(I3-I1))/I2;
    dydt(9) = (M(3) + w(1)*w(2)*(I1-I2))/I3;
    % Break up the ERP matrix equation:
    % edot = 0.5*M*omega
    % e1dot = (e4*w1 - e3*w2 + e2*w3)/2
    % e2dot = (e3*w1 + e4*w2 - e1*w3)/2
    % e3dot = (-e2*w1 + e1*w2 + e4*w3)/2
    % e4dot = (-e1*w1 - e2*w2 - e3*w3)/2
    dydt(10) = (e(4)*w(1) - e(3)*w(2) + e(2)*w(3))/2;
    dydt(11) = (e(3)*w(1) + e(4)*w(2) - e(1)*w(3))/2;
    dydt(12) = (-e(2)*w(1) + e(1)*w(2) + e(4)*w(3))/2;
    dydt(13) = (-e(1)*w(1) - e(2)*w(2) - e(3)*w(3))/2;
end


%% Input validation functions

function checkIfFinite(x)
% This function will be used for input validation for finite numbers
    if (isinf(x))
        error(inputname(1) + " must be a finite number");
    end
end


function checkIfPositive(x, allowZero, allowInf)
% This function will be used for input validation for positive numbers
    % Account for missing arguments
    if (nargin < 3)
        allowInf = false;
    end
    if (nargin < 2)
        allowZero = false;
    end
    % Check if nonnegative or positive
    if (allowZero)
        criteria = "nonnegative";
    else
        criteria = "positive";
    end
    valid = or((x > 0), and(allowZero, (x == 0)));
    if (~valid)
        error(inputname(1) + " must be " + criteria + "!");
    end
    % Check if infinite (based on inputs)
    if (~allowInf)
        checkIfFinite(x);
    end
end


function checkIfInRange(x, l, u)
% This function will be used for input validation for bounded numbers
    if or((x < l), (x > u))
        error(inputname(1) + " must be between " + l + " and " + u);
    end
end


%% Helper functions

function graphFinalState(rel_tip_final, e_final, L_arrow)
% This function will graph a 3D plot of the arrow's final state relative to
% the target bulls-eye
    % Find arrow vector, remembering that the arrow's z is along the shaft
    arrow_vec = rotVecWithQuat([0;0;L_arrow], e_final);
    rel_tail_final = rel_tip_final - arrow_vec;
    figure();
    scatter3(0,0,0, 36,"green","filled") % Target
    hold on;
    quiver3(rel_tail_final(1),rel_tail_final(2),rel_tail_final(3), ...
        arrow_vec(1),arrow_vec(2),arrow_vec(3), "off", "Color","red");
    hold off;
    legend(["Bulls-eye"; "Arrow"], 'Location','northeast');
    % Adjust xlim so arrowhead doesn't look weird
    c_xlim = xlim; % Current
    n_xlim = rel_tip_final(1) + [-1, 1]; % Narrow
    xlim([min(n_xlim(1),c_xlim(1)), max(n_xlim(2),c_xlim(2))]);
    % Adjust ylim so plot ends at target
    c_ylim = ylim;
    ylim([c_ylim(1), min(0,c_ylim(2))]);
    % Adjust zlim so plot ends at ground if it falls short of the target
    if (rel_tip_final(2) < 0)
        c_zlim = zlim;
        zlim([rel_tip_final(3), c_zlim(2)])
    end
    xlabel("x (Inertial - Target) [ft]");
    ylabel("y (Inertial - Target) [ft]");
    zlabel("z (Inertial - Target) [ft]");
end


function fig = graph3D_wComps(x_vecs, t_vec, x_name, x_units, t_units, ...
    doHold, fig, leg_text)
% This function will graph the components of angular velocity over time as
% well as the path in 3D space made by the vector over time
    % Account for missing arguments
    use_leg = (nargin >= 8);
    if (nargin < 7)
        fig = NaN;
    end
    if (nargin < 6)
        doHold = false;
    end
    if (nargin < 5)
        t_units = "s";
    end
    if (nargin < 4)
        x_units = "ft";
    end
    if (nargin < 3)
        x_name = "x";
    end
    do_time_plot = (nargin >= 2);
    % Split x_vecs into its 3 components
    [x1, x2, x3] = deal(x_vecs(:,1), x_vecs(:,2), x_vecs(:,3));
    % Plot 3D path
    % Check if fig is actually a figure
    % (https://stackoverflow.com/questions/14812999/check-if-matlab-handle-
    % is-a-figure-handle)
    if (ishandle(fig) && strcmp(get(fig,'type'),'figure'))
        figure(fig);
    else
        fig = figure();
    end
    if (do_time_plot)
        subplot(3,2,[1,3,5])
    end
    plot3(x1, x2, x3, "-");
    title("3D Path")
    xlabel(x_name + "_1 [" + x_units + "]");
    ylabel(x_name + "_2 [" + x_units + "]");
    zlabel(x_name + "_3 [" + x_units + "]");
    if (doHold)
        hold on;
    else
        hold off;
    end
    if (use_leg)
        legend(leg_text, 'Location','northeast');
    end
    if (do_time_plot)
        for i = 1:3
            subplot(3,2,i*2)
            plot(t_vec, x_vecs(:,i), "-");
            title(x_name + "_" + i + " Over Time")
            xlabel("t [" + t_units + "]");
            ylabel(x_name + "_" + i + " [" + x_units + "]");
            xlim([t_vec(1), t_vec(end)]);
            if (doHold)
                hold on;
            else
                hold off;
            end
        end
    end
end


function R_mat = Rx(theta, inDeg)
% This function performs a rotation about the x axis, given the angle (in
% degrees by default)
    if (nargin < 2)
        inDeg = true;
    end
    if (~inDeg)
        theta = rad2deg(theta);
    end
    [c, s] = deal(cosd(theta), sind(theta));
    R_mat = [1, 0, 0;...
             0, c, -s;...
             0, s, c];
end


function R_mat = Ry(theta, inDeg)
% This function performs a rotation about the y axis, given the angle (in
% degrees by default)
    if (nargin < 2)
        inDeg = true;
    end
    if (~inDeg)
        theta = rad2deg(theta);
    end
    [c, s] = deal(cosd(theta), sind(theta));
    R_mat = [c, 0, s;...
             0, 1, 0;...
             -s, 0, c];
end


function R_mat = Rz(theta, inDeg)
% This function performs a rotation about the z axis, given the angle (in
% degrees by default)
    if (nargin < 2)
        inDeg = true;
    end
    if (~inDeg)
        theta = rad2deg(theta);
    end
    [c, s] = deal(cosd(theta), sind(theta));
    R_mat = [c, -s, 0;...
             s, c, 0;...
             0, 0, 1];
end


function u = makeUnitVec(v)
% This function makes a unit vector from the given vector
% NOTE: If v is a vector of zeros, the returned vector will be v
    v_mag = norm(v);
    if (v_mag == 0)
        u = v;
    else
        u = v/norm(v);
    end
end


function q = vecAng2Quat(vec, ang, inDeg)
% This function converts a vector-angle to quaternion
% Note that the vector must be a unit vector!
    % Account for missing arguments
    if (nargin < 3)
        inDeg = true;
    end
    if (~inDeg)
        ang = rad2deg(ang);
    end
    % Initialize quat
    q = zeros(4,1);
    % Populate quat
    q(1:3) = vec*sind(ang/2);
    q(4) = cosd(ang/2);
end


function q = rotMat2Quat(M)
% This function finds the quaternion that represents the rotation matrix M
    % Now convert to vector-angle pair
    ang = acosd((trace(M) - 1)/2);
    s_ang = sind(ang);
    if (s_ang == 0)
        V_cross = zeros(3,3);
    else
        V_cross = 0.5*(M - M')/sind(ang);
    end
    vec = [V_cross(3,2); V_cross(1,3); V_cross(2,1)];
    % Now convert to a quaternion
    q = vecAng2Quat(vec, ang);
end


function res = multTwoQuats(q1, q2)
% This function multiplies two quaternions
    [q1_vec, q2_vec, q1_4, q2_4] = deal(q1(1:3), q2(1:3), q1(4), q2(4));
    res = zeros(size(q1));
    res(1:3) = q1_4*q2_vec + q2_4*q1_vec + cross(q1_vec, q2_vec);
    res(4) = q1_4*q2_4 - dot(q1_vec, q2_vec);
end


function q_T = invQuat(q)
% This function inverts a quaternion
    q_T = q;
    q_T(1:3) = -q_T(1:3);
end


function res = rotVecWithQuat(vec, quat)
% This function rotates a vector using a quaternion
    % Make the vector able to be multiplied by a quaternion
    vec(4) = 0;
    % Make the inverse of quat
    q_T = invQuat(quat);
    % Perform rotation: q*v*qT
    res = multTwoQuats(multTwoQuats(quat, vec), q_T);
    res = res(1:3);
end


function g_MC = findGravityM(quat, g, z_COM, z_tip_COM, L_shaft, L_fin, ...
    N_fins, m_tip, m_shaft, m_fin, m_nock)
% This function finds the moment on the arrow from gravity, in the
% body-fixed frame (C)
    % Find gravity in C
    gC = rotVecWithQuat([0;0;-g], invQuat(quat));
    % Prepare for r x m*g
    z_axis_body = [0;0;1];
    r_tip_COM = (z_tip_COM - z_COM)*z_axis_body;
    r_shaft_COM = (L_shaft/2 - z_COM)*z_axis_body;
    r_tail_COM = (-z_COM)*z_axis_body;
    % Since the fins have the same mass (m), we know that (f_1 x m*g) +
    % (f_2 x m*g) + ... = (f_1 + f_2 + ...) x m*g
    % We also know that the fins are equally spaced, so the off-z
    % components of that sum will cancel out
    r_fin_COM_sum = N_fins*(L_fin/3 - z_COM)*z_axis_body;
    % (r1 x m1*g) + (r2 x m2*g) + ... = (m1*r1 x g) + (m2*r2 x g) + ... = 
    % (m1*r1 + m2*r2 + ...) x g, which minimizes number of cross products
    mr_sum = m_tip*r_tip_COM + m_shaft*r_shaft_COM + ...
        m_nock*r_tail_COM + m_fin*r_fin_COM_sum;
    g_MC = cross(mr_sum, gC);
end


function D = calcDrag(Cd, A, rho, v)
% This function calculates the drag force acting on an area, accounting for
% the sign from the velocity (of the object):
% - v(i) > 0 -> D(i) < 0
% - v(i) < 0 -> D(i) > 0
% This way, D is easier to use in vectors
    v_mag = norm(v);
    if (v_mag == 0)
        D = zeros(size(v));
    else
        % Drag acts in the opposite direction of velocity
        D = -0.5*Cd*A*rho*dot(v,v)*(v/v_mag);
    end
end


function [fin_DF, fin_MC] = findDragFM(quat, omega, v_arrow, ...
    dir_wind, v_wind, ang_fin, L_fin, h_fin, N_fins, z_COM, OR_shaft, ...
    L_shaft, L_tip, rho_air)
% Finds the drag forces and moments (about the arrow's COM) in the
% inertial (F) and body-fixed (C) frames, respectively.
% In F, y is along the archer-target vector, and z is the altitude vector
% In C, z is along the arrow shaft, and y is to the right of the arrow when
% looking at it from the tail
    % Prepare to find the total drag and moment vectors in body frame
    [fin_DC, fin_MC] = deal(zeros(3,1));
    % From the front and back, the arrow's frontal area is pi*OR_shaft^2
    % From https://www.engineeringtoolbox.com/drag-coefficient-d_627.html,
    % we get the following drag coefficients:
    % - Streamlined body (arrow front): 0.04
    % - Wires and cables (arrow shaft side): 1.0-1.3 -> 1.15
    % - Solid hemisphere flow normal to curved side (side of tip): 0.47
    % - Solid hemisphere flow normal to flat side (arrow back): 1.17
    % - Squared flat plate at 90 deg (fins normal to wind): 1.17
    % We need to figure out the relative velocity of the arrow relative to
    % the wind, then figure out the areas affected by drag. For instance,
    % if the arrow is traveling slower than the wind, drag will push the
    % arrow forward.
    v_rel = v_arrow - v_wind*[dir_wind; 0]; % Wind is 2D
    % Find the arrow's orientation
    z_axis_body = [0;0;1];
    % Find the drag acting on each arrow component
    % What we could do is transform the relative velocity into the
    % body-fixed frame, so that we could not only calculate the torques in
    % the body frame, but also simplify the code that figures out the
    % directions of the applied forces:
    % 1. Use the inverse of e to get v_rel to body-fixed
    % 2. Figure out the drag vector on each arrow component
    % 3. Transform those vectors back to inertial using e
    % 4. Use the x, y, and z components for the acceleration calculations
    v_rel_bf = rotVecWithQuat(v_rel, invQuat(quat));
    % Now, v_rel_bf is along the axes of [nsym1, nsym2, sym]
    cs_area_arrow = pi*OR_shaft^2;
    % However, we need to account for the velocity of the tip due to the
    % rotation of the arrow relative to inertial space
    r_tip_COM = (L_shaft - z_COM + L_tip/3)*z_axis_body;
    r_tail_COM = (-z_COM)*z_axis_body;
    vrbf_tip = v_rel_bf + cross(omega, r_tip_COM);
    vrbft_f = vrbf_tip(3); % Tip's forwards relative velocity in body frame
    if (vrbft_f >= 0)
        % Arrow is traveling faster than wind in its forward direction
        Cd_arrow = 0.04;
        v_Cd_arrow = vrbft_f;
    else
        Cd_arrow = 1.17;
        % We need the velocity of the arrow's tail
        v_Cd_arrow = dot(v_rel_bf + cross(omega, r_tail_COM), z_axis_body);
    end
    D_arrow_3 = calcDrag(Cd_arrow, cs_area_arrow, rho_air, v_Cd_arrow);
    D_arrow = [0; 0; D_arrow_3]; % Produces no moment
    fin_DC = fin_DC + D_arrow;
    % For the shaft, both the nsym axes are important
    cs_area_shaft = 2*L_shaft*OR_shaft;
    Cd_shaft = 1.15;
    r_shaft_COM = (L_shaft/2 - z_COM)*z_axis_body;
    vrbfs = (v_rel_bf.*[1;1;0]) + cross(omega, r_shaft_COM);
    D_shaft = calcDrag(Cd_shaft, cs_area_shaft, rho_air, vrbfs);
    fin_DC = fin_DC + D_shaft;
    fin_MC = fin_MC + cross(r_shaft_COM, D_shaft);
    % Also account for the tip
    cs_area_tip = 0.5*(2*OR_shaft)*L_tip;
    Cd_tip_side = 0.47;
    D_tip = calcDrag(Cd_tip_side, cs_area_tip, rho_air, vrbf_tip(1:2));
    D_tip(3) = 0; % Already accounted for in D_arrow
    fin_DC = fin_DC + D_tip;
    fin_MC = fin_MC + cross(r_tip_COM, D_tip);
    % Now it's time for the fins
    z_fin_COM = L_fin/3 - z_COM;
    r_fin_COM_sp = h_fin/3 + OR_shaft; % Scalar and in plane of fin COMs
    % First, account for drag caused by forwards motion (if angled fins)
    Cd_face = 1.17; % A flat plane normal to flow
    s_af = sind(ang_fin);
    c_af = cosd(ang_fin);
    % The fins exert a normal force perpendicular to the surface, so
    % that means there's a component pushing back. The normal force (n)
    % comes from vrbf_for*sind(ang_fin), and we can further decompose
    % that normal force into a sideways force (s = -n*cosd(ang_fin)) and
    % backwards force (b = -n*sind(ang_fin))
    % We can assume drag D = -n
    A_fin = L_fin*h_fin/2; % Right triangle
    % Now we need the fin COM relative to arrow COM
    fin_COM_1 = [r_fin_COM_sp; 0; z_fin_COM]; % True for 3 and 4 fins
    ang_rot = 360/N_fins;
    % First find COMs and Ns
    % Remember that positive fin angle causes spin || pointing
    fin_N_1 = [0;-c_af;s_af];
    % Now find the drag vectors
    for i = 1:N_fins
        q_i = vecAng2Quat(z_axis_body, (i-1)*ang_rot);
        r_COM_i = rotVecWithQuat(fin_COM_1, q_i);
        N_i = rotVecWithQuat(fin_N_1, q_i);
        % Accounting for relative velocity from rotation at fin's COM
        fin_V_i = v_rel_bf + cross(omega, r_COM_i);
        % We need only to examine the velocity component along the fin's
        % normal vector (which is a unit vector)
        V_n_i = dot(fin_V_i, N_i)*N_i;
        D_i = calcDrag(Cd_face, A_fin, rho_air, V_n_i);
        M_i = cross(r_COM_i, D_i);
        fin_DC = fin_DC + D_i;
        fin_MC = fin_MC + M_i;
    end
    % Although the moments are in the body-fixed frame, we need drag in the
    % inertial frame
    fin_DF = rotVecWithQuat(fin_DC, quat);
end


function [pI_vec, z_arrow_COM, z_tip_COM, m_shaft, m_nock, m_tip, ...
    m_fin, L_arrow] = getArrowMOI(L_shaft, OR_shaft, IR_shaft, ...
    GPI_shaft, W_nock, W_tip, ang_tip, L_tip_in, L_fin, h_fin, W_fin, ...
    N_fins, g_imp, lbf2grain)
% This function calculates the arrow's moment of inertia in lbm•in^2
% Inputs are mostly covered in the main function's description
% g_imp and lbf2grain are constants defined in the main function

    % Shaft (hollow cylinder)
    W_shaft = (GPI_shaft * L_shaft)/lbf2grain; % lbf
    m_shaft = W_shaft/g_imp; % lbm
    sumRSq = OR_shaft^2 + IR_shaft^2; % in^2
    I3_shaft = 0.5*m_shaft*sumRSq; % lbm•in^2
    I1_shaft = (m_shaft/12)*(3*sumRSq + L_shaft^2); % lbm•in^2

    % Nock (assumed point mass along shaft axis of symmetry)
    m_nock = W_nock/(lbf2grain*g_imp); % lbm

    % Tip (solid cone with a solid cylinder insert)
    m_tip = W_tip/(lbf2grain*g_imp); % lbm
    L_tip_cone = OR_shaft/tand(ang_tip/2); % in
    % Find density to get m_tip_cone and m_tip_in
    V_tip_cone = pi*(OR_shaft^2)*L_tip_cone/3; % in^3
    V_tip_in = pi*(IR_shaft^2)*L_tip_in;
    rho_tip = m_tip/(V_tip_cone + V_tip_in); % lbm/in^3
    m_tip_cone = rho_tip * V_tip_cone;
    m_tip_in = rho_tip * V_tip_in;
    % Now calculate MOIs
    I3_tip_cone = (3/10)*m_tip_cone*(OR_shaft^2); % lbm•in^2
    I1_tip_cone = (3*m_tip_cone/80)*(4*(OR_shaft^2) + L_tip_cone^2);
    I3_tip_in = 0.5*m_tip_in*(IR_shaft^2);
    I1_tip_in = (m_tip_in/12)*(3*(IR_shaft^2) + L_tip_in^2);
    
    % Fins (point masses at triangle COMs)
    % Since the fins act as point masses at h_fin/3 above the shaft, we can
    % find the MOI of the whole fin system using summations:
    % I3 = Izz = ∑(m*(x^2 + y^2))
    % I1 = Ixx = ∑(m*(z^2 + y^2))
    % I2 = Iyy = ∑(m*(z^2 + x^2))
    % Ixy = ∑(m*(x*y))
    % z is the axis along the arrow's shaft
    % Since z = 0, no need to check for Ixz and Iyz
    % Using axes where one of the point-masses is on an axis guarantees
    % principal axes (Ixy = 0) for N_fin in [3, 4]
    % Also, I1 = I2
    m_fin = W_fin/(lbf2grain*g_imp); % lbm
    R_fin = OR_shaft + (h_fin/3);
    I3_fins = N_fins*m_fin*(R_fin^2); % lbm•in^2
    I1_fins = I3_fins/2; % This is how the math worked out

    % Assemble the arrow
    % Find COM first (we know it's on the I3 axis, so focus only on that)
    m_fins = N_fins*m_fin;
    m_arrow = m_nock + m_fins + m_shaft + m_tip;
    L_arrow = L_shaft + L_tip_cone;
    % Let the nock (end of shaft) be the reference 0
    z_tip_in_COM = L_shaft - L_tip_in/2;
    z_tip_cone_COM = L_shaft + L_tip_cone/3;
    z_arrow_COM = (m_fins*(L_fin/3) + m_shaft*(L_shaft/2) + ...
        m_tip_in*(z_tip_in_COM) + m_tip_cone*(z_tip_cone_COM))/m_arrow;
    % Find MOIs wrt total COM
    % Since the fins are the only components with a restriction on the
    % principal axes' orientation, we can just add the moments
    I3_arrow = I3_fins + I3_shaft + I3_tip_in + I3_tip_cone; % lbm•in^2
    I1_arrow = m_nock*(z_arrow_COM)^2 + ...
        (I1_fins + m_fins*(z_arrow_COM - L_fin/3)^2) + ...
        (I1_shaft + m_shaft*(z_arrow_COM - L_shaft/2)^2) + ...
        (I1_tip_in + m_tip_in*(z_arrow_COM - (z_tip_in_COM))^2) + ...
        (I1_tip_cone + m_tip_cone*(z_arrow_COM - (z_tip_cone_COM))^2);
    % Arrow is close enough to axisymmetric
    pI_vec = [I1_arrow; I1_arrow; I3_arrow];

    % Don't forget to find the COM of the tip (with insert)
    z_tip_COM = (m_tip_in*z_tip_in_COM + m_tip_cone*z_tip_cone_COM)/m_tip;
end