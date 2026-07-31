import math
import numpy as np
from numpy.typing import NDArray
from warnings import warn
from typing import Callable, Optional

from ..utils import common as cmn

"""
This Abkowitz-type dynamics module follows the seminar paper
by Youjun Yang and Ould el Moctar (2024)
https://doi.org/10.1016/j.oceaneng.2024.116927
"""

# Type aliases to make function signatures clearer
h_over_T = float
unit_type = float

def _shallow_water_correction_factory(
    alpha: float, 
    a: float, 
    b: float, 
    c: float) -> Callable[[h_over_T, unit_type], float]:

    """
    Factory from Eq (8) to create the shallow water correction functions
    of Table 3 from Youjun Yang and Ould el Moctar (2024).
    """
    def correction(h_over_T: float, x: float) -> float:
        """Shallow water correction polynomial

        Args:
            h_over_T (float): Water depth to vessel draft ratio
            x (float): Generic variable (e.g., speed)

        Returns:
            float: Correction factor
        """
        return (a*x**2 + b*x + c) * h_over_T**(-alpha) + 1.0

    return correction

# Create shallow water correction functions
f_Xu_dot = _shallow_water_correction_factory(alpha=1.717, a=0.0,      b=0.0,     c=3.552)
f_Xu     = _shallow_water_correction_factory(alpha=3.760, a=0.186,    b=1.013,   c=1.423)
f_Xv     = _shallow_water_correction_factory(alpha=3.661, a=-664.667, b=0.0,     c=21.096)
f_Xvr    = _shallow_water_correction_factory(alpha=2.615, a=-52.017,  b=-18.278, c=10.088)
f_Xd     = _shallow_water_correction_factory(alpha=1.760, a=-0.369,   b=-0.865,  c=-0.264)

f_Yv_dot = _shallow_water_correction_factory(alpha=2.809, a=0.0,      b=0.0,     c=4.933)
f_Yr_dot = _shallow_water_correction_factory(alpha=2.847, a=0.0,      b=0.0,     c=4.297)
f_Yv     = _shallow_water_correction_factory(alpha=4.315, a=96.490,   b=24.697,  c=8.843)
f_Yr     = _shallow_water_correction_factory(alpha=6.785, a=17.440,   b=30.944,  c=7.680)
f_Yvr    = _shallow_water_correction_factory(alpha=6.003, a=-239.093, b=218.306, c=148.738)
f_Yd     = _shallow_water_correction_factory(alpha=4.288, a=-1.107,   b=-2.781,  c=1.178)

f_Nv_dot = _shallow_water_correction_factory(alpha=2.846, a=0.0,      b=0.0,     c=4.292)
f_Nr_dot = _shallow_water_correction_factory(alpha=2.724, a=0.0,      b=0.0,     c=2.936)
f_Nv     = _shallow_water_correction_factory(alpha=2.966, a=47.870,   b=0.0,     c=6.268)
f_Nr     = _shallow_water_correction_factory(alpha=4.353, a=1.267,    b=2.962,   c=3.615)
f_Nvr    = _shallow_water_correction_factory(alpha=4.721, a=0.0,      b=0.0,     c=68.509)
f_Nd     = _shallow_water_correction_factory(alpha=4.720, a=-1.218,   b=-3.050,  c=0.167)

class AbkowitzModel(cmn.Maneuvervable):
    """
    Wrapper for an Abkowitz-type dynamics model 
    following Yang and el Moctar (2024) for shallow water 
    conditions for inland vessels.
    
    :var vessel: AbkowitzVessel
        Vessel parameters according to the Abkowitz model.
    """

    def __init__(self, vessel: cmn.AbkowitzVessel):
        self.vessel = vessel

    def dynamics(
        self,
        X: np.ndarray, 
        psi:float,
        delta: float, 
        h: float,
        fl_psi: float | None, 
        fl_vel: float | None) -> np.ndarray:
        """System of ODEs after Yang and el Moctar (2024)

        Args:
            X (np.ndarray): 
                Initial values (u, v, r) at current timestep
            psi (float): 
                vessel heading in the global frame
            delta (float): 
                Rudder angle [rad]
            h (float): 
                Water depth [m]
            fl_psi (float): 
                Attack angle of current relative
                longitudinal axis of motion [rad]
            fl_vel (float): 
                Velocity of current [m/s]
            w_vel (float): 
                Wind velocity [m/s] (inactive, only for signature matching)
            beta_w (float): 
                Wind attack angle relative to heading [rad] (inactive, only for signature matching)

        Returns:
            np.ndarray: Motion derivatives (u_dot, v_dot, r_dot)
        """
        

        # Shorten the parameter dict to avoid clutter
        p = self.vessel


        u, v, r  = X[0], X[1], X[2]

        # Correct for current if given
        if fl_vel is not None and fl_psi is not None:
            u_c = fl_vel * math.cos(fl_psi - psi - math.pi)
            v_c = fl_vel * math.sin(fl_psi - psi - math.pi)
            
            # Relative velocities
            u_r = u - u_c
            v_r = v - v_c

            u = u_r
            v = v_r
        
        # Difference from refernence speed
        du = u - p.U_0

        # Water depth to draft ratio
        h_prime = min(h / p.T, 8.00)  # Cap at max value from Yang and el Moctar (2024)
        
        # Instantaneous speed
        U = math.sqrt(u**2 + v**2) + 1e-6 # Add small number to avoid division by zero in corrections at zero speed
        
        # Non-dimensionalized quantities
        v_prime = v / U
        r_prime = r * p.L / U
        du_prime = du / U
        
        # Model forces from Eq (7)
        # Surge force without added mass couplings
        F_X = (
            (p.X_u * du_prime + p.X_uu * du_prime**2) * f_Xu(h_prime, du_prime) +
            (p.X_vv * v_prime**2 + p.X_vvvv * v_prime**4) * f_Xv(h_prime, abs(v_prime)) +
            (p.X_rr * r_prime**2 + p.X_vr * v_prime * r_prime + p.X_vvrr * v_prime**2 * r_prime**2) * f_Xvr(h_prime, abs(v_prime)) +
            (p.X_dd * delta**2 + p.X_udd * du_prime * delta**2) * f_Xd(h_prime, du_prime)
        )
        
        # Sway force without added mass couplings
        F_Y = (
            (p.Y_v * v_prime + p.Y_vvv * v_prime**3) * f_Yv(h_prime, abs(v_prime)) +
            (p.Y_r * r_prime + p.Y_rrr * r_prime**3) * f_Yr(h_prime, abs(r_prime)) +
            (p.Y_vrr * v_prime * r_prime**2 + p.Y_vvr * v_prime**2 * r_prime) * f_Yvr(h_prime, abs(v_prime)) +
            (p.Y_d * delta + p.Y_ddd * delta**3 + p.Y_ud * du_prime * delta + p.Y_uddd * du_prime * delta**3) * f_Yd(h_prime, du_prime)
        )
        
        # Yaw moment without added mass couplings
        F_N = (
            (p.N_v * v_prime + p.N_vvv * v_prime**3 + p.N_vrr * v_prime * r_prime**2) * f_Nv(h_prime, abs(v_prime)) +
            (p.N_r * r_prime + p.N_rrr * r_prime**3) * f_Nr(h_prime, abs(r_prime)) +
            (p.N_vvr * v_prime**2 * r_prime) * f_Nvr(h_prime, abs(v_prime)) +
            (p.N_d * delta + p.N_ddd * delta**3 + p.N_ud * du_prime * delta + p.N_uddd * du_prime * delta**3) * f_Nd(h_prime, du_prime)
        )

        # Dimensionalizers
        D_FORCE  = 0.5 * p.rho * U**2 * p.L**2
        D_MOMENT = 0.5 * p.rho * U**2 * p.L**3

        
        # Re-dimensionalize forces and moments
        F_X = F_X * D_FORCE
        F_Y = F_Y * D_FORCE
        F_N = F_N * D_MOMENT
        
        # Moment of inertia in yaw from radius of gyration
        I_zz = p.m * p.r_zz**2
        
        # Re-dimensionalize added mass couplings
        am_surge  = p.X_udot * 0.5 * p.rho * (p.L**3) * f_Xu_dot(h_prime, 0)
        am_sway_v = p.Y_vdot * 0.5 * p.rho * (p.L**3) * f_Yv_dot(h_prime, 0)
        am_sway_r = p.Y_rdot * 0.5 * p.rho * (p.L**4) * f_Yr_dot(h_prime, 0)
        am_yaw_v  = p.N_vdot * 0.5 * p.rho * (p.L**4) * f_Nv_dot(h_prime, 0)  
        am_yaw_r  = p.N_rdot * 0.5 * p.rho * (p.L**5) * f_Nr_dot(h_prime, 0)
        
        # Left-hand side matrix from Eq (6) plus added mass terms
        # minus all non-derivative forces and correction for shallow water
        M = np.array(
            [
                [p.m - am_surge, 0.0,                    0.0                    ],
                [0.0,            p.m - am_sway_v,        p.m * p.x_G - am_sway_r],
                [0.0,            p.m * p.x_G - am_yaw_v, I_zz - am_yaw_r        ]
            ]
        )
        
        # Right-hand side from Eq (6) plus all non-derivative forces
        F = np.array(
            [
                F_X + p.m * v * r + p.m * p.x_G * (r**2),
                F_Y - p.m * u * r,
                F_N - p.m * p.x_G * u * r
            ]
        )
        
        return np.linalg.solve(M, F)

    def step(
        self,
        *,
        X: np.ndarray,
        delta: float,
        psi: float,
        water_depth: float,
        fl_psi: Optional[float] = None,
        fl_vel: Optional[float] = None,
    ) -> np.ndarray:
        """Compute body-fixed accelerations for one instantaneous vessel state.

        This method evaluates the ODE right-hand side only (no time integration).

        Args:
            X (np.ndarray): Body-fixed state vector ``[u, v, r]``.
            delta (float): Rudder angle in radians.
            psi (float): Vessel heading in the earth-fixed frame, radians.
            water_depth (float): Water depth in meters.
            fl_psi (Optional[float]): Current direction in earth-fixed frame, radians.
            fl_vel (Optional[float]): Current speed in m/s.

        Raises:
            cmn.LogicError:
                If ``fl_vel`` is set but ``fl_psi`` is missing, or if
                ``water_depth < vessel.T``.

        Returns:
            np.ndarray: ``[u_dot, v_dot, r_dot]`` in the body-fixed frame.
        """

        if fl_vel is not None and fl_psi is None:
            raise cmn.LogicError(
                "No current direction specified. Use the `fl_psi` keyword to set it."
            )
        if fl_vel is None and fl_psi is not None:
            warn(
                "Current attack angle is given but current velocity is "
                "turned off. Attack angle is ignored."
            )

        # Correct for shallow water if a water depth is given.
        # If none is given, open water with infinite depth is assumed
        if water_depth is not None:
            if water_depth < self.vessel.T:
                raise cmn.LogicError(
                    "Water depth cannot be less than ship draft.\n"
                    f"Water depth: {np.round(water_depth, 3)} | Ship draft: {self.vessel.T}"
                )

        # Develop equation of motion
        uvr_dot = self.dynamics(
            X=X,
            delta=delta,
            h=water_depth,
            psi=psi,
            fl_psi=fl_psi,
            fl_vel=fl_vel,
        )

        # Return derivatives
        return uvr_dot

    def pstep(
        self,
        *,
        X: np.ndarray,
        pos: cmn.EarthFixedPos,
        dT: float,
        delta: float,
        psi: float,
        water_depth: float,
        fl_psi: Optional[float] = None,
        fl_vel: Optional[float] = None,
        mode = cmn.IntegrationMode.TRAPEZOIDAL
    ) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
        """Advance velocity and position by one time step.

        Integrates body-fixed dynamics and updates earth-fixed pose.

        Args:
            X (np.ndarray): Body-fixed velocity state ``[u, v, r]``.
            pos (cmn.EarthFixedPos): Current earth-fixed position ``[N, E]``.
            dT (float): Time step in seconds.
            delta (float): Rudder angle in radians.
            psi (float): Heading angle in radians.
            water_depth (float): Water depth in meters.
            fl_psi (Optional[float]): Current direction in radians.
            fl_vel (Optional[float]): Current speed in m/s.
            mode (cmn.IntegrationMode): Integration method (TRAPEZOIDAL or RK4).

        Returns:
            tuple[NDArray[np.float32], NDArray[np.float32]]:
                ``(new_uvr, new_eta)`` where ``new_uvr = [u, v, r]`` and
                ``new_eta = [N, E, psi]``.
        """
        if mode is cmn.IntegrationMode.TRAPEZOIDAL:
            return self._impl_pstep_trapezoidal(
                X=X,
                pos=pos,
                dT=dT,
                delta=delta,
                psi=psi,
                water_depth=water_depth,
                fl_psi=fl_psi,
                fl_vel=fl_vel
            )
        elif mode is cmn.IntegrationMode.RK4:
            return self._impl_pstep_rk4(
                X=X,
                pos=pos,
                dT=dT,
                delta=delta,
                psi=psi,
                water_depth=water_depth,
                fl_psi=fl_psi,
                fl_vel=fl_vel
            )
        else:
            raise ValueError(f"Unknown integration mode: {mode}")

    def _impl_pstep_trapezoidal(
        self,
        *,
        X: np.ndarray,
        pos: cmn.EarthFixedPos,
        dT: float,
        delta: float,
        psi: float,
        water_depth: float,
        fl_psi: Optional[float] = None,
        fl_vel: Optional[float] = None,
    ) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
        """Trapezoidal position update with Euler velocity update.

        Scheme:
            1) Compute ``uvr_dot`` from the current state.
            2) Update ``[u, v, r]`` with forward Euler.
            3) Update heading and earth-fixed position with trapezoidal averaging.

        Returns:
            tuple[NDArray[np.float32], NDArray[np.float32]]:
                Updated ``(new_uvr, new_eta)``, with ``new_eta = [N, E, psi]``.
        """
        old_uvr = X.copy()
        
        # 1. Get body-fixed accelerations
        # (Assuming self.step now correctly returns strictly uvr_dot)
        uvr_dot = self.step(
            X=X,
            delta=delta,
            psi=psi,
            fl_psi=fl_psi,
            fl_vel=fl_vel,
            water_depth=water_depth,
        )

        # 2. Integrate body velocities (Forward Euler)
        new_uvr = old_uvr + uvr_dot * dT
        
        # 3. Integrate heading (Trapezoidal)
        r_old = old_uvr[2]
        r_new = new_uvr[2]
        new_psi = psi + 0.5 * (r_old + r_new) * dT

        # 4. Calculate Earth-fixed velocities
        # Abkowitz rotpsi returns [d_xi (East), d_eta (North), d_psi]
        old_earth_dot = np.dot(cmn.rotpsi(psi,np.pi/2), old_uvr)
        new_earth_dot = np.dot(cmn.rotpsi(new_psi,np.pi/2), new_uvr)

        # Map Abkowitz outputs to North and East
        # xi is horizontal (East), eta is vertical (North)
        dot_E_old, dot_N_old = old_earth_dot[0], old_earth_dot[1]
        dot_E_new, dot_N_new = new_earth_dot[0], new_earth_dot[1]

        # 5. Integrate Earth-fixed position (Trapezoidal)
        new_N = pos[0] + 0.5 * (dot_N_old + dot_N_new) * dT
        new_E = pos[1] + 0.5 * (dot_E_old + dot_E_new) * dT
        
        # 6. Assemble state and normalize angle
        new_eta = np.array([new_N, new_E, new_psi], dtype=np.float32)
        new_eta[2] = cmn.angle_to_two_pi(new_eta[2])

        return new_uvr.astype(np.float32), new_eta
    
    def create_6dof_dynamic(
        self, 
        delta: float, 
        fl_psi: float, 
        fl_vel: float, 
        water_depth: float
    ) -> Callable[[NDArray[np.float32]], NDArray[np.float32]]:
        """Create a closure for 6-state RK4 integration.

        The generated function maps
        ``state = [N, E, psi, u, v, r]`` to
        ``state_dot = [N_dot, E_dot, psi_dot, u_dot, v_dot, r_dot]``
        under fixed environmental/control inputs.
        """
        def frozen_dynamic(state: NDArray[np.float32]) -> NDArray[np.float32]:
            # 1. Unpack the state
            N, E, psi, u, v, r = state
            uvr = np.array([u, v, r], dtype=np.float32)

            # 2. Get Kinetics (Body-fixed accelerations: u_dot, v_dot, r_dot)
            uvr_dot = self.dynamics(
                X=uvr, 
                delta=delta, 
                h=water_depth, 
                psi=psi, 
                fl_psi=fl_psi, 
                fl_vel=fl_vel
            )

            # 3. Get Kinematics (Earth-fixed velocities)
            # Remember: Abkowitz rotpsi returns [d_East, d_North, d_psi]
            earth_dot = np.dot(cmn.rotpsi(psi, np.pi/2), uvr)
            
            dot_E = earth_dot[0]
            dot_N = earth_dot[1]
            dot_psi = earth_dot[2]

            # 4. Return the combined state derivative
            # Must match the order of the input state: [N, E, psi, u, v, r]
            return np.array([dot_N, dot_E, dot_psi, uvr_dot[0], uvr_dot[1], uvr_dot[2]], dtype=np.float32)

        return frozen_dynamic

    def _impl_pstep_rk4(
        self,
        *,
        X: np.ndarray,
        pos: cmn.EarthFixedPos,
        dT: float,
        delta: float,
        psi: float,
        water_depth: float,
        fl_psi: Optional[float] = None,
        fl_vel: Optional[float] = None,
    ) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
        """RK4 integration of combined pose and velocity state.

        Builds a 6-state vector ``[N, E, psi, u, v, r]``, integrates one RK4
        step, normalizes heading to ``[0, 2π)``, and splits back into
        ``(new_uvr, new_eta)``.
        """

        # 0. Concatenate eta and nu
        full_state = np.hstack([pos,psi,X])
    
        # 1. Generate the frozen function with current environmental params
        ext_dyn = self.create_6dof_dynamic(delta, fl_psi, fl_vel, water_depth)
        
        # 2. Perform the RK4 integration on the entire 6-DOF state
        new_state = self.rk4_step(ext_dynamic=ext_dyn, full_state=full_state, dT=dT)
        
        # 3. Normalize the heading angle (index 2 is psi)
        new_state[2] = cmn.angle_to_two_pi(new_state[2])
        
        return new_state[3:], new_state[:3] # Returns the updated [N, E, psi, u, v, r]

    def rk4_step(
        self,*, 
        ext_dynamic: Callable, 
        full_state: NDArray[np.float32],
        dT: float) -> tuple[NDArray[np.float32]]:
        """Run one classical RK4 step for an arbitrary state vector.

        Args:
            ext_dynamic (Callable): Function ``f(x)`` returning ``dx/dt``.
            full_state (NDArray[np.float32]): Current state vector.
            dT (float): Time step in seconds.

        Returns:
            NDArray[np.float32]: State after one RK4 increment.
        """
        k1 = ext_dynamic(full_state)
        k2 = ext_dynamic(full_state + 0.5 * dT * k1)
        k3 = ext_dynamic(full_state + 0.5 * dT * k2)
        k4 = ext_dynamic(full_state + dT * k3)
        return full_state + (dT / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4) 
