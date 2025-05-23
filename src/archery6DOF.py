import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def archery6DOF(L_shaft, OR_shaft, IR_shaft, GPI_shaft, W_nock, W_tip, ang_tip, L_tip_in, L_fin, \
        h_fin, W_fin, N_fins, ang_fin, v_arrow, h_arrow, ang_aim_lr, ang_aim_ud, d_target, \
        h_target, v_wind=0.0, dir_wind=np.zeros(2), rho_air=0.0752):
    # This function will run a 6DOF simulation of the firing of an arrow.
    # Inputs:
    # - L_shaft: length of arrow shaft [in]
    #   - Must be positive
    # - OR_shaft: outer radius of arrow shaft [in]
    #   - Must be positive
    # - IR_shaft: inner radius of arrow shaft [in]
    #   - Must be positive (or 0 for a solid arrow)
    # - GPI_shaft: grains per inch of arrow shaft [GPI]
    #   - 1 grain = 1/7000 lbf. Grains per inch is the industry standard
    #   - Must be positive
    # - W_nock: weight of arrow nock [grain]
    #   - Must be nonnegative
    # - W_tip: weight of arrow tip [grain]
    #   - Must be positive
    # - ang_tip: angle made by a cross-section of the tip [deg]
    #   - Range: 15º to 45º
    # - L_tip_in: length of the arrow tip's insert [in]
    #   - Range: 0 in to L_shaft
    # - L_fin: length of a single tail fin along the shaft [in]
    #   - Must be nonnegative
    # - h_fin: maximum height of a single tail fin above the shaft [in]
    #   - Must be nonnegative
    # - ang_fin: off-shaft angle of a single tail fin [deg]
    #   - 0º: fins are perfectly aligned with shaft -> no spin during flight
    #   - Range: -10º to 10º (positive angle will cause the arrow's spin vector
    #     to be in the same direction as its pointing vector)
    # - W_fin: weight of a single fin [grain]
    #   - Must be nonnegative
    # - N_fins: number of equally-spaced fins [-]
    #   - Range: 3 to 4
    # - v_arrow: speed of the arrow upon leaving the bow [ft/s]
    #   - There are a lot of parameters needed to calculate this from bow
    #     specifications, so I might add them in later
    #   - Range: 100 ft/s to 500 ft/s
    # - h_arrow: height of the arrow's tail upon leaving the bow [ft]
    #   - Must be high enough such that the tip of the arrow isn't touching
    #     the ground (h = 0)
    # - ang_aim_lr: left-right aiming angle [deg]
    #   - Angle between the aimed arrow's pointing vector and the plane
    #     perpendicular to the ground that contains the archer and target
    #   - Positive angles go right
    # - ang_aim_ud: up-down aiming angle [deg]
    #   - Angle between the aimed arrow's pointing vector and the ground plane
    #   - Positive angles go up
    # - d_target: distance of the target [ft]
    #   - Should be positive, but I'll allow negative for trick shots that take
    #     advantage of strong winds
    # - h_target: height of the target's bulls-eye [ft]
    #   - Must be positive
    # - v_wind: wind speed [ft/s]
    #   - Must be nonnegative
    #   - Default value is 0 ft/s
    # - dir_wind: a 2D column vector that gives the direction of wind [-]
    #   - Let y-hat be the vector pointing from the archer to the target
    #   - Let x-hat be perpendicular to y-hat and pointing to the right of the
    #     archer
    #     - This makes z-hat = cross(x-hat, y-hat), which isn't included in
    #       this 2D vector
    #   - Doesn't have to be a unit vector
    #   - Default vector is [0;0]
    # - rho_air: air density [lbm/ft^3] (optional)
    #   - I've seen bow specs measured with Imperial Units, so as much as it
    #     annoys me, I won't use SI for consistency's sake
    #   - Must be positive
    #   - Default value is 0.0752 lbm/ft^3 (1.204 kg/m^3)
    # Outputs:
    # - diff_coords: a 3D column vector giving the differences in coordinates
    #   of the arrow relative to the target when either 1) the arrow crosses
    #   the plane of the target (containing x-hat and z-hat) OR 2) the arrow
    #   hits the ground [ft]
    #   - The coordinate frame of the vector will be based on the three hat
    #     vectors defined earlier, so getting a negative y value means that you
    #     undershot the target
    # Assumptions:
    # - Constant air density, wind speed, and gravity over the arrow's flight
    # - The archer and target are in an inertial frame
    # - Arrow is a streamlined body for the purposes of drag
    #   (https://www.engineeringtoolbox.com/drag-coefficient-d_627.html)
    # - Moments from drag on a component can be modeled by having the drag
    #   force act on the component's center of mass
    # - The archer won't impart any rotation on the arrow when firing, such
    #   that the initial velocity vector matches the arrow's pointing vector
    # - Fins are right triangles with negligible thickness with their bases at
    #   at the end of the arrow shaft
    # - The nock's center of mass is located at the end of the shaft (since
    #   there is an inserted part)
    # - The nock and tip are close enough to axisymmetric
    # - The nock's length is negligible for purposes of calculating drag

    # Unit conversions and constants
    ft2in = 12.0
    lbf2grain = 7000.0
    g_imp = 32.174  # ft/s^2

    # Normalize wind direction
    dir_wind = dir_wind / np.linalg.norm(dir_wind) if np.linalg.norm(dir_wind) > 0 else np.zeros(2)

    # Validate inputs
    check_positive(h_target)
    check_range(ang_fin, -10, 10)
    check_range(N_fins, 3, 4)
    check_range(v_arrow, 100, 500)
    check_positive(v_wind, allow_zero=True)
    check_positive(L_fin, allow_zero=True)
    check_positive(h_fin, allow_zero=True)
    check_positive(L_shaft)
    check_range(L_tip_in, 0, L_shaft)
    check_range(ang_tip, 15, 45)
    check_positive(W_tip)
    check_positive(W_nock, allow_zero=True)
    check_positive(IR_shaft, allow_zero=True)
    check_positive(OR_shaft)
    check_positive(GPI_shaft)
    check_positive(rho_air)

    # Get moments, masses, length, COM
    (pI_vec, z_arrow_COM, z_tip_COM, m_shaft, m_nock, m_tip, m_fin, L_arrow) = get_arrow_moi( \
        L_shaft, OR_shaft, IR_shaft, GPI_shaft, W_nock, W_tip, ang_tip, L_tip_in, L_fin, h_fin, \
        W_fin, N_fins, g_imp, lbf2grain)
    m_arrow = m_shaft + m_nock + m_tip + m_fin * N_fins
    I1, I2, I3 = pI_vec
    print(f"[I1, I2, I3] = [{I1:.6f}, {I2:.6f}, {I3:.6f}] lb•in^2")
    print(f"COM is {z_arrow_COM:.6f} inches from the back of the arrow")
    print(f"The arrow's mass is {m_arrow:.6f} lbm")
    # Convert to ft
    I1, I2, I3 = pI_vec / (ft2in ** 2)
    z_arrow_COM /= ft2in
    z_tip_COM /= ft2in
    L_arrow /= ft2in
    L_shaft /= ft2in
    L_tip = L_arrow - L_shaft
    L_fin /= ft2in
    h_fin /= ft2in
    OR_shaft /= ft2in
    L_com2tip = L_arrow - z_arrow_COM

    # Initial orientation and position
    aim_dir = np.array([
        np.sin(np.radians(ang_aim_lr)) * np.cos(np.radians(ang_aim_ud)),
        np.cos(np.radians(ang_aim_lr)) * np.cos(np.radians(ang_aim_ud)),
        np.sin(np.radians(ang_aim_ud))
    ])
    com2tip0_F = L_com2tip * aim_dir
    com2tip0_C = np.array([0.0, 0.0, L_com2tip])
    tip0 = np.array([0.0, 0.0, h_arrow]) + L_arrow * aim_dir
    if tip0[2] <= 0:
        raise ValueError("Arrow tip must be above ground at launch.")
    # Summarize in state vector
    x0 = tip0 - com2tip0_F
    v0 = v_arrow * aim_dir
    w0 = np.zeros(3) # No spin at first
    e0 = rot_mat_to_quat(Rz(-90) @ Ry(-90) @ Rx(-ang_aim_lr) @ Ry(ang_aim_ud))
    y0 = np.concatenate([x0, v0, w0, e0])
    
    # Time setup
    v0_z = v0[2]
    t_max = 2 * (np.sqrt(v0_z**2 + 2*h_arrow*g_imp) - v0_z)
    t_span = (0, t_max)
    
    # Solve ODE
    odelam = lambda t, y: odefun(t, y, I1, I2, I3, dir_wind, v_wind, ang_fin, L_fin, h_fin, \
                                N_fins, z_arrow_COM, z_tip_COM, OR_shaft, L_shaft, L_tip, \
                                rho_air, g_imp, m_tip, m_shaft, m_fin, m_nock)
    events_lam = lambda t, y: events_fun(t, y, L_com2tip, d_target)
    events_lam.terminal = True
    sol = solve_ivp(odelam, t_span, y0, method="RK45", rtol=1e-8, atol=1e-10, events=events_lam)

    if (sol.status == 0):
        print("Warning: Arrow still in flight and short of target plane when time limit reached.")
    
    # Return coordiante difference
    final_state = sol.y[:, -1]
    com_coords_final = final_state[0:3]
    e_final = final_state[9:13]
    com2tip1 = rot_vec_with_quat(com2tip0_C, e_final)
    tip_coords_final = com_coords_final + com2tip1
    target_coords = np.array([0.0, d_target, h_target])
    diff_coords = tip_coords_final - target_coords

    # Assume y_vecs is a NumPy array of shape (N, 13)
    x_vecs = sol.y[0:3, :]
    v_vecs = sol.y[3:6, :]
    w_vecs = sol.y[6:9, :]
    e_vecs = sol.y[9:13, :]

    # Compute tip positions
    tip_vecs = np.zeros_like(x_vecs)
    for i in range(len(sol.t)):
        x_i = x_vecs[:, i]
        e_i = e_vecs[:, i]
        tip_vecs[:, i] = x_i + rot_vec_with_quat(com2tip0_C, e_i)

    # Plot
    plt.ion()
    leg_text = ["Arrow COM", "Arrow Tip"]
    graph3D_w_comps([x_vecs, tip_vecs], sol.t, "x", "ft", "s", leg_text)
    graph3D_w_comps([v_vecs], sol.t, "v", "ft/s")
    graph3D_w_comps([w_vecs], sol.t, r"$\omega$", "rad/s")

    # Final visualization
    graph_final_state(diff_coords, e_final, L_arrow)

    # Show all plots
    plt.ioff()
    plt.show()

    return diff_coords


# ODE functions

def events_fun(t, y, L_com2tip, y_target):
    # Extract quaternion from y ([x; v; w; e]), then find current location of tip
    e1 = y[9:13]
    com2tip1 = rot_vec_with_quat(np.array([0, 0, L_com2tip]), e1)
    y_tip = y[1] + com2tip1[1]
    z_tip = y[2] + com2tip1[2]
    # Event 1: crossing the target plane
    value1 = y_tip - y_target
    # Event 2: hitting the ground
    value2 = z_tip
    # Since both events are terminal, and can be approached from either direction, combine them
    # to minimize number of quaternion calculations
    return value1*value2


def odefun(t, y, I1, I2, I3, dir_wind, v_wind, ang_fin, L_fin, h_fin, N_fins, z_arrow_COM, \
           z_tip_COM, OR_shaft, L_shaft, L_tip, rho_air, g, m_tip, m_shaft, m_fin, m_nock):
    # y has 13 elements: 3 for x, 3 for v, 3 for wdot, and 4 for edot
    dydt = np.zeros_like(y)
    v = y[3:6]
    w = y[6:9]
    e = y[9:13]
    # xdot = v
    dydt[0:3] = v
    # vdot = a from drag + gravity
    arrow_DF, arrow_MC = find_drag_FM(e, w, v, dir_wind, v_wind, ang_fin, L_fin, h_fin, N_fins, \
                                  z_arrow_COM, OR_shaft, L_shaft, L_tip, rho_air)
    dydt[3:6] = arrow_DF
    dydt[5] -= g # Subtract gravity from vertical acceleration (remember 0-indexing)
    # Angular acceleration, which requires us to solve for torques
    g_MC = find_gravity_M(e, g, z_arrow_COM, z_tip_COM, L_shaft, L_fin, N_fins, m_tip, m_shaft, \
                          m_fin, m_nock)
    sum_MC = arrow_MC + g_MC
    # Break up Euler equations
    dydt[6] = (sum_MC[0] + w[1]*w[2]*(I2 - I3)) / I1
    dydt[7] = (sum_MC[1] + w[0]*w[2]*(I3 - I1)) / I2
    dydt[8] = (sum_MC[2] + w[0]*w[1]*(I1 - I2)) / I3
    # Quaternion derivative
    # Break up the ERP matrix equation: edot = 0.5*M*omega
    dydt[9]  = 0.5 * ( e[3]*w[0] - e[2]*w[1] + e[1]*w[2] )
    dydt[10] = 0.5 * ( e[2]*w[0] + e[3]*w[1] - e[0]*w[2] )
    dydt[11] = 0.5 * (-e[1]*w[0] + e[0]*w[1] + e[3]*w[2] )
    dydt[12] = 0.5 * (-e[0]*w[0] - e[1]*w[1] - e[2]*w[2] )
    # Return state derivative
    return dydt


# Input verification functions

def check_positive(x, allow_zero=False):
    if not np.isfinite(x):
        raise ValueError("Input must be finite.")
    if (x < 0) or (not allow_zero and x == 0):
        raise ValueError("Input must be positive" + (" or zero." if allow_zero else "."))


def check_range(x, lo, hi):
    if x < lo or x > hi:
        raise ValueError(f"Input must be in range [{lo}, {hi}].")


# Helper functions

def graph_final_state(rel_tip_final, e_final, L_arrow, show_plot=False):
    # This function will graph a 3D plot of the arrow's final state relative to
    # the target bulls-eye

    # Find arrow vector, remembering that the arrow's z is along the shaft
    arrow_vec = rot_vec_with_quat([0, 0, L_arrow], e_final)
    rel_tail_final = rel_tip_final - arrow_vec
    # Initialize axes
    fig = plt.figure(layout="constrained")
    ax = fig.add_subplot(111, projection='3d')
    # Plot the target
    ax.scatter(0, 0, 0, color='green', s=100, label="Bulls-eye")
    # Plot the arrow as a vector from tail to tip
    ax.quiver(
        rel_tail_final[0], rel_tail_final[1], rel_tail_final[2],
        arrow_vec[0], arrow_vec[1], arrow_vec[2],
        color='red', label="Arrow", arrow_length_ratio=0.1
    )
    # Adjust x-limits
    current_xlim = ax.get_xlim()
    new_xlim = [rel_tip_final[0] - 1, rel_tip_final[0] + 1]
    ax.set_xlim([min(new_xlim[0], current_xlim[0]), max(new_xlim[1], current_xlim[1])])
    # Adjust y-limits
    current_ylim = ax.get_ylim()
    ax.set_ylim([current_ylim[0], min(0, current_ylim[1])])
    # Adjust z-limits if arrow hit below target height
    if (rel_tip_final[1] < 0):
        current_zlim = ax.get_zlim()
        ax.set_zlim([rel_tip_final[2], current_zlim[1]])
    # Labels and legend
    ax.zaxis.set_rotate_label(False)
    ax.set_xlabel("x (Inertial - Target) [ft]")
    ax.set_ylabel("y (Inertial - Target) [ft]")
    ax.set_zlabel("z (Inertial - Target) [ft]")
    ax.legend(loc='upper right')
    if not show_plot:
        plt.show()


def graph3D_w_comps(x_vecs_list, t_vec=None, x_name="x", x_units="ft", t_units="s", \
                    leg_text=None, show_plot=False):
    # This function will graph the components of angular velocity over time as
    # well as the path in 3D space made by the vector over time

    do_time_plot = t_vec is not None
    use_leg = leg_text is not None
    num_x_vecs = len(x_vecs_list)

    # Initialize figure and axes depending on inputs
    fig = plt.figure(layout="constrained")
    if do_time_plot:
        gs = GridSpec(3, 2, figure=fig)
        ax_3d = fig.add_subplot(gs[:, 0], projection='3d')
    else:
        ax_3d = fig.add_subplot(111, projection='3d')

    # 3D path plot
    for i in range(num_x_vecs):
        x_vecs = x_vecs_list[i]
        ax_3d.plot(x_vecs[0,:], x_vecs[1,:], x_vecs[2,:], label=leg_text[i] if use_leg else None)
    # Format plot
    ax_3d.set_title("3D Path")
    ax_3d.zaxis.set_rotate_label(False)
    ax_3d.set_xlabel(f"{x_name}_1 [{x_units}]")
    ax_3d.set_ylabel(f"{x_name}_2 [{x_units}]")
    ax_3d.set_zlabel(f"{x_name}_3 [{x_units}]")
    if use_leg:
        ax_3d.legend(loc='upper right')

    # Component plots over time
    if do_time_plot:
        for i in range(3):
            ax = fig.add_subplot(gs[i, 1])
            for ii in range(num_x_vecs):
                x_vecs = x_vecs_list[ii]
                ax.plot(t_vec, x_vecs[i, :])
            # Format plot
            ax.set_title(f"{x_name}_{i+1} Over Time")
            ax.set_xlabel(f"t [{t_units}]")
            ax.set_ylabel(f"{x_name}_{i+1} [{x_units}]")
            ax.set_xlim(t_vec[0], t_vec[-1])
    # Show plot if requested
    if not show_plot:
        plt.show()


def Rx(theta_deg):
    theta = np.radians(theta_deg)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def Ry(theta_deg):
    theta = np.radians(theta_deg)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def Rz(theta_deg):
    theta = np.radians(theta_deg)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def vec_ang_to_quat(vec, ang, in_deg=True):
    # This function converts a vector-angle to quaternion
    # Note that the vector must be a unit vector!
    if in_deg:
        ang = np.radians(ang)
    q = np.zeros(4)
    q[:3] = vec * np.sin(ang/2)
    q[3] = np.cos(ang/2)
    return q


def rot_mat_to_quat(M):
    # This function finds the quaternion that represents the rotation matrix M
    # Convert to vector-angle pair since there's a function that converts that to a quaternion
    ang = np.arccos((np.trace(M) - 1) / 2)
    sin_ang = np.sin(ang)
    if np.isclose(sin_ang, 0):
        V_cross = np.zeros((3, 3)) # Avoid division by zero or instability
    else:
        V_cross = 0.5 * (M - M.T) / sin_ang
    vec = np.array([V_cross[2,1], V_cross[0,2], V_cross[1,0]])
    return vec_ang_to_quat(vec, ang, False)


def mult_two_quats(q1, q2):
    # This function multiplies two quaternions
    q1_vec = np.array(q1[:3])
    q2_vec = np.array(q2[:3])
    q1_4 = q1[3]
    q2_4 = q2[3]
    # Initialize then populate result
    res = np.zeros_like(q1)
    res[:3] = q1_4*q2_vec + q2_4*q1_vec + np.cross(q1_vec, q2_vec)
    res[3] = q1_4*q2_4 - np.dot(q1_vec, q2_vec)
    return res


def inv_quat(q):
    # This function inverts a quaternion
    # q_T = [-q1, -q2, -q3, q4]
    return np.concatenate((-q[:3], [q[3]]))


def rot_vec_with_quat(vec, quat):
    # This function rotates a vector using a quaternion
    vec_q = np.append(vec, 0) # Convert vector to quaternion format (x, y, z, 0)
    q_inv = inv_quat(quat) # Make the inverse of quat
    # Perform rotation: q*v*qT
    res_q = mult_two_quats(mult_two_quats(quat, vec_q), q_inv)
    return res_q[:3]


def find_gravity_M(quat, g, z_COM, z_tip_COM, L_shaft, L_fin, N_fins, m_tip, m_shaft, m_fin, \
                   m_nock):
    # This function finds the moment on the arrow from gravity, in the
    # body-fixed frame (C)

    # Gravity vector in C
    gC = rot_vec_with_quat(np.array([0, 0, -g]), inv_quat(quat))
    # Prepare for r x m*g
    z_axis_body = np.array([0, 0, 1])
    r_tip_COM = (z_tip_COM - z_COM) * z_axis_body
    r_shaft_COM = (L_shaft / 2 - z_COM) * z_axis_body
    r_tail_COM = (-z_COM) * z_axis_body
    # Since the fins have the same mass (m), we know that (f_1 x m*g) + (f_2 x m*g) + ... =
    # (f_1 + f_2 + ...) x m*g
    # We also know that the fins are equally spaced, so the off-z
    # components of that sum will cancel out
    r_fin_COM_sum = N_fins * (L_fin / 3 - z_COM) * z_axis_body
    # (r1 x m1*g) + (r2 x m2*g) + ... = (m1*r1 x g) + (m2*r2 x g) + ... =
    # (m1*r1 + m2*r2 + ...) x g, which minimizes number of cross products
    mr_sum = m_tip*r_tip_COM + m_shaft*r_shaft_COM + m_nock*r_tail_COM + m_fin*r_fin_COM_sum
    return np.cross(mr_sum, gC)



def calc_drag(Cd, A, rho, v):
    # This function calculates the drag force acting on an area, accounting for
    # the sign from the velocity (of the object):
    # - v(i) > 0 -> D(i) < 0
    # - v(i) < 0 -> D(i) > 0
    # This way, D is easier to use in vectors
    if (not isinstance(v, np.ndarray)):
        v = np.array(v) # In case a scalar or list is passed in
    v_mag = np.linalg.norm(v)
    if v_mag == 0:
        return np.zeros_like(v)
    else:
        # Drag acts in the opposite direction of velocity
        return -0.5 * Cd * A * rho * np.dot(v, v) * (v / v_mag)


def find_drag_FM(quat, omega, v_arrow, dir_wind, v_wind, ang_fin, L_fin, h_fin, N_fins, \
        z_COM, OR_shaft, L_shaft, L_tip, rho_air):
    # Finds the drag forces and moments (about the arrow's COM) in the
    # inertial (F) and body-fixed (C) frames, respectively.
    # In F, y is along the archer-target vector, and z is the altitude vector
    # In C, z is along the arrow shaft, and y is to the right of the arrow when
    # looking at it from the tail

    # From the front and back, the arrow's frontal area is pi*OR_shaft^2
    # From https://www.engineeringtoolbox.com/drag-coefficient-d_627.html,
    # we get the following drag coefficients:
    # - Streamlined body (arrow front): 0.04
    # - Wires and cables (arrow shaft side): 1.0-1.3 -> 1.15
    # - Solid hemisphere flow normal to curved side (side of tip): 0.47
    # - Solid hemisphere flow normal to flat side (arrow back): 1.17
    # - Squared flat plate at 90 deg (fins normal to wind): 1.17

    # Prepare to find the total drag and moment vectors in body frame
    fin_DC = np.zeros(3)
    fin_MC = np.zeros(3)

    # We need to figure out the relative velocity of the arrow relative to
    # the wind, then figure out the areas affected by drag. For instance,
    # if the arrow is traveling slower than the wind, drag will push the
    # arrow forward.
    v_rel = v_arrow - v_wind * np.append(dir_wind, 0) # Wind is 2D

    # Body z-axis
    z_axis_body = np.array([0, 0, 1])

    # Find the drag acting on each arrow component
    # What we could do is transform the relative velocity into the
    # body-fixed frame, so that we could not only calculate the torques in
    # the body frame, but also simplify the code that figures out the
    # directions of the applied forces:
    # 1. Use the inverse of e to get v_rel to body-fixed
    # 2. Figure out the drag vector on each arrow component
    # 3. Transform those vectors back to inertial using e
    # 4. Use the x, y, and z components for the acceleration calculations
    v_rel_bf = rot_vec_with_quat(v_rel, inv_quat(quat))
    # Now, v_rel_bf is along the axes of [nsym1, nsym2, sym]
    cs_area_arrow = np.pi * OR_shaft**2

    # However, we need to account for the velocity of the tip due to the
    # rotation of the arrow relative to inertial space
    r_tip_COM = (L_shaft - z_COM + L_tip/3) * z_axis_body
    r_tail_COM = -z_COM * z_axis_body
    vrbf_tip = v_rel_bf + np.cross(omega, r_tip_COM)
    vrbft_f = vrbf_tip[2] # Tip's forwards relative velocity in body frame
    if (vrbft_f >= 0):
        # Arrow is traveling faster than wind in its forward direction
        Cd_arrow = 0.04
        v_Cd_arrow = vrbft_f
    else:
        Cd_arrow = 1.17
        # We need the velocity of the arrow's tail
        v_Cd_arrow = np.dot(v_rel_bf + np.cross(omega, r_tail_COM), z_axis_body)
    D_arrow_3 = calc_drag(Cd_arrow, cs_area_arrow, rho_air, v_Cd_arrow)
    D_arrow = np.array([0, 0, D_arrow_3]) # Produces no moment
    fin_DC += D_arrow

    # For the shaft, both the nsym axes are important
    cs_area_shaft = 2 * L_shaft * OR_shaft
    Cd_shaft = 1.15
    r_shaft_COM = (L_shaft/2 - z_COM) * z_axis_body
    vrbfs = v_rel_bf * np.array([1, 1, 0]) + np.cross(omega, r_shaft_COM)
    D_shaft = calc_drag(Cd_shaft, cs_area_shaft, rho_air, vrbfs)
    fin_DC += D_shaft
    fin_MC += np.cross(r_shaft_COM, D_shaft)

    # Also account for the tip
    cs_area_tip = 0.5*(2*OR_shaft)*L_tip
    Cd_tip_side = 0.47
    D_tip = calc_drag(Cd_tip_side, cs_area_tip, rho_air, vrbf_tip[:2])
    D_tip = np.append(D_tip, 0) # Dimension 3 is already accounted for in D_arrow
    fin_DC += D_tip
    fin_MC += np.cross(r_tip_COM, D_tip)

    # Fins
    z_fin_COM = L_fin/3 - z_COM
    r_fin_COM_sp = h_fin/3 + OR_shaft # Scalar and in plane of fin COMs
    # First, account for drag caused by forwards motion (if angled fins)
    Cd_face = 1.17 # A flat plane normal to flow
    s_af = np.sin(np.radians(ang_fin))
    c_af = np.cos(np.radians(ang_fin))
    # The fins exert a normal force perpendicular to the surface, so
    # that means there's a component pushing back. The normal force (n)
    # comes from vrbf_for*sind(ang_fin), and we can further decompose
    # that normal force into a sideways force (s = -n*cosd(ang_fin)) and
    # backwards force (b = -n*sind(ang_fin))
    # We can assume drag D = -n
    A_fin = L_fin * h_fin / 2 # Right triangle

    # Now we need the fin COM relative to arrow COM
    fin_COM_1 = np.array([r_fin_COM_sp, 0, z_fin_COM]) # True for 3 and 4 fins
    # First find COMs and Ns
    # Remember that positive fin angle causes spin || pointing
    fin_N_1 = np.array([0, -c_af, s_af])
    ang_rot = 360 / N_fins
    # Now find the drag vectors
    for i in range(N_fins):
        q_i = vec_ang_to_quat(z_axis_body, (i) * ang_rot)
        r_COM_i = rot_vec_with_quat(fin_COM_1, q_i)
        N_i = rot_vec_with_quat(fin_N_1, q_i)
        # Accounting for relative velocity from rotation at fin's COM
        fin_V_i = v_rel_bf + np.cross(omega, r_COM_i)
        # We need only to examine the velocity component along the fin's
        # normal vector (which is a unit vector)
        V_n_i = np.dot(fin_V_i, N_i) * N_i
        D_i = calc_drag(Cd_face, A_fin, rho_air, V_n_i)
        M_i = np.cross(r_COM_i, D_i)
        fin_DC += D_i
        fin_MC += M_i
    # Although the moments are in the body-fixed frame, we need drag in the
    # inertial frame
    fin_DF = rot_vec_with_quat(fin_DC, quat)

    return fin_DF, fin_MC


def get_arrow_moi(L_shaft, OR_shaft, IR_shaft, GPI_shaft, W_nock, W_tip, ang_tip, L_tip_in, \
        L_fin, h_fin, W_fin, N_fins, g_imp, lbf2grain):
    # This function calculates the arrow's moment of inertia in lbm•in^2
    # Inputs are mostly covered in the main function's description
    # g_imp and lbf2grain are constants defined in the main function

    # Shaft (hollow cylinder)
    W_shaft = (GPI_shaft * L_shaft) / lbf2grain # lbf
    m_shaft = W_shaft/g_imp # lbm
    sum_r_sq = OR_shaft**2 + IR_shaft**2 # in^2
    I3_shaft = 0.5 * m_shaft * sum_r_sq # lbm•in^2
    I1_shaft = (m_shaft/12) * (3*sum_r_sq + L_shaft**2) # lbm•in^2

    # Nock (assumed point mass along shaft axis of symmetry)
    m_nock = W_nock / (lbf2grain * g_imp) # lbm

    # Tip (solid cone with a solid cylinder insert)
    m_tip = W_tip / (lbf2grain * g_imp) # lbm
    L_tip_cone = OR_shaft / np.tan(np.radians(ang_tip/2)) # in
    # Find density to get m_tip_cone and m_tip_in
    V_tip_cone = np.pi * (OR_shaft**2) * L_tip_cone / 3 # in^3
    V_tip_in = np.pi * (IR_shaft**2) * L_tip_in # in^3
    rho_tip = m_tip / (V_tip_cone + V_tip_in) # lbm/in^3
    m_tip_cone = rho_tip * V_tip_cone
    m_tip_in = rho_tip * V_tip_in
    # Now calculate MOIs
    I3_tip_cone = (3/10) * m_tip_cone * (OR_shaft**2)
    I1_tip_cone = (3 * m_tip_cone / 80) * (4 * (OR_shaft**2) + L_tip_cone**2)
    I3_tip_in = 0.5 * m_tip_in * (IR_shaft**2)
    I1_tip_in = (m_tip_in/12) * (3 * (IR_shaft**2) + L_tip_in**2)

    # Fins (point masses at triangle COMs)
    # Since the fins act as point masses at h_fin/3 above the shaft, we can find the MOI of the
    # whole fin system using summations:
    # I3 = Izz = ∑(m*(x^2 + y^2))
    # I1 = Ixx = ∑(m*(z^2 + y^2))
    # I2 = Iyy = ∑(m*(z^2 + x^2))
    # Ixy = ∑(m*(x*y))
    # z is the axis along the arrow's shaft
    # Since z = 0, no need to check for Ixz and Iyz
    # Using axes where one of the point-masses is on an axis guarantees principal axes (Ixy = 0)
    # for N_fin in [3, 4]
    # Also, I1 = I2
    m_fin = W_fin / (lbf2grain * g_imp) # lbm
    R_fin = OR_shaft + (h_fin/3)
    I3_fins = N_fins * m_fin * (R_fin**2)
    I1_fins = I3_fins / 2 # This is how the math worked out

    # Assemble the arrow
    # Find COM first (we know it's on the I3 axis, so focus only on that)
    m_fins = N_fins * m_fin
    m_arrow = m_nock + m_fins + m_shaft + m_tip
    L_arrow = L_shaft + L_tip_cone
    # Let the nock (end of shaft) be the reference 0
    z_tip_in_COM = L_shaft - (L_tip_in/2)
    z_tip_cone_COM = L_shaft + (L_tip_cone/3)
    z_arrow_COM = (m_fins * (L_fin/3) + \
                   m_shaft * (L_shaft/2) + \
                   m_tip_in * z_tip_in_COM + \
                   m_tip_cone * z_tip_cone_COM) / m_arrow
    # Find MOIs wrt total COM
    # Since the fins are the only components with a restriction on the
    # principal axes' orientation, we can just add the moments
    I3_arrow = I3_fins + I3_shaft + I3_tip_in + I3_tip_cone
    I1_arrow = (m_nock * (z_arrow_COM**2) + \
                (I1_fins + m_fins * (z_arrow_COM - (L_fin/3))**2) + \
                (I1_shaft + m_shaft * (z_arrow_COM - (L_shaft/2))**2) + \
                (I1_tip_in + m_tip_in * (z_arrow_COM - z_tip_in_COM)**2) + \
                (I1_tip_cone + m_tip_cone * (z_arrow_COM - z_tip_cone_COM)**2))
    # Arrow is close enough to axisymmetric
    pI_vec = np.array([I1_arrow, I1_arrow, I3_arrow])
    # Don't forget to find the COM of the tip (with insert)
    z_tip_COM = ((m_tip_in * z_tip_in_COM) + (m_tip_cone * z_tip_cone_COM)) / m_tip

    return pI_vec, z_arrow_COM, z_tip_COM, m_shaft, m_nock, m_tip, m_fin, L_arrow


# Main

def main():
    # Run test case to check with MATLAB version
    diff_coords = archery6DOF(29, 0.295, 0.245, 8.2, 10, 100, 30, 0.5, 1, 0.5, 6, 3, 5, 125, 4, \
                              0, 3.4, 60, 4, 2, np.array([1,2]))
    print("Python:", diff_coords, sep="\n")
    # Print expected results from MATLAB
    diff_coords_M = np.array([-0.006523698913424, 0, 0.033347616796036])
    print("MATLAB:", diff_coords_M, sep="\n")


if __name__ == "__main__":
    main()