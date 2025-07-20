import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

# Constants
G = 6.67430e-11
c = 299792458
m1 = 30 * 1.9885e30
m2 = 30 * 1.9885e30
m = m1 + m2
nu = (m1 * m2) / m**2
R = 500e6 * 3.086e16
p0 = 20 * G * m / c**2
p_isco_wave = 3 * G * m / c**2
p_isco_orbit = 0 * p_isco_wave

color1 = 'tab:orange'
color2 = 'tab:blue'

def dpdt(t, p, e):
    return -(64/5) * nu * c * (G * m / (c**2 * p))**3 * (1 - e**2)**1.5 * (1 + (7/8) * e**2)

def dedt(t, p, e):
    return -(304/15) * nu * c * e / p * (G * m / (c**2 * p))**3 * (1 - e**2)**1.5 * (1 + (121/304) * e**2)

def dpsidt(t, psi, p, e):
    return np.sqrt(G * m / p**3) * (1 + e * np.cos(psi))**2

def h_plus(t, psi, e, h0):
    return -h0 * (2 * np.cos(2 * psi) + e * np.cos(psi) + 2 * e * np.cos(psi)**3 + e**2)

def solve_with_isco(e0, t_max, dt, p_isco_stop):
    def odes(t, y):
        p, e, psi = y
        return [dpdt(t, p, e), dedt(t, p, e), dpsidt(t, psi, p, e)]

    def isco_event(t, y):
        return y[0] - p_isco_stop
    isco_event.terminal = True
    isco_event.direction = -1

    y0 = [p0, e0, 0]
    t_eval = np.arange(0, t_max, dt)
    sol = solve_ivp(odes, (0, t_max), y0, t_eval=t_eval, events=isco_event, rtol=1e-9, atol=1e-9)
    return sol

# Parameters
e0 = 0
t_max = 400
dt = 0.0005

# Solve orbits and waveforms
sol_orbit = solve_with_isco(e0, t_max, dt, p_isco_orbit)
p_orb, e_orb, psi_orb = sol_orbit.y
t_orb = sol_orbit.t
r_orb = p_orb / (1 + e_orb * np.cos(psi_orb))
x_orb = r_orb * np.cos(psi_orb)
y_orb = r_orb * np.sin(psi_orb)

sol_wave = solve_with_isco(e0, t_max, dt, p_isco_wave)
p_wave, e_wave, psi_wave = sol_wave.y
t_wave = sol_wave.t
h0_wave = 2 * G**2 * m1 * m2 / (c**4 * R * p_wave)
hplus = h_plus(t_wave, psi_wave, e_wave, h0_wave)

# Prepare figure
fig, axs = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [2, 1]})

# Orbit axes
axs[0].set_title(f'Orbit with decaying trail for $e_0 = {e0}$')
axs[0].set_xlabel(r'$x~[Gm/c^2]$')
axs[0].set_ylabel(r'$y~[Gm/c^2]$')
axs[0].set_aspect('equal')
axs[0].grid(True)

scale = G * m / c**2
axs[0].set_xlim(-1.2 * np.max(x_orb) / scale, 1.2 * np.max(x_orb) / scale)
axs[0].set_ylim(-1.2 * np.max(x_orb) / scale, 1.2 * np.max(x_orb) / scale)

# Waveform axes
axs[1].set_title(r'Gravitational Waveform $h_+(t)$')
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel(r'$h_+(t)$')
axs[1].grid(True)
axs[1].set_xlim(0, 1.05 * t_wave[-1])
axs[1].set_ylim(1.1 * np.min(hplus), 1.1 * np.max(hplus))

# Initialize orbit trails and points
trail_length = 500  # Number of points in the trail

trail1_line, = axs[0].plot([], [], color=color1, alpha=0.5)
trail2_line, = axs[0].plot([], [], color=color2, alpha=0.5)
star1_point, = axs[0].plot([], [], 'o', color=color1, markersize=8)
star2_point, = axs[0].plot([], [], 'o', color=color2, markersize=8)

# Initialize waveform line
wave_line, = axs[1].plot([], [], color='tab:orange')

def init():
    trail1_line.set_data([], [])
    trail2_line.set_data([], [])
    star1_point.set_data([], [])
    star2_point.set_data([], [])
    wave_line.set_data([], [])
    return [trail1_line, trail2_line, star1_point, star2_point, wave_line]

def update(frame):
    idx_start = max(0, frame - trail_length)
    idx_end = min(frame, len(x_orb))
    idxs = np.arange(idx_start, idx_end)

    x_scaled = x_orb[idxs] / scale
    y_scaled = y_orb[idxs] / scale

    trail1_line.set_data(0.5 * x_scaled, 0.5 * y_scaled)
    trail2_line.set_data(-0.5 * x_scaled, -0.5 * y_scaled)

    if frame < len(x_orb):
        star1_point.set_data([+0.5 * x_orb[frame] / scale], [+0.5 * y_orb[frame] / scale])
        star2_point.set_data([-0.5 * x_orb[frame] / scale], [-0.5 * y_orb[frame] / scale])
        star1_point.set_color(color1)
        star2_point.set_color(color2)
        star1_point.set_markersize(8)
        star2_point.set_markersize(8)
    elif frame == len(x_orb):
        # Collision moment: show black hole at center
        star1_point.set_data([0], [0])
        star2_point.set_data([], [])
        star1_point.set_color('black')
        star1_point.set_markersize(25)
    else:
        # After merger: keep only the black hole at center
        star1_point.set_data([0], [0])
        star2_point.set_data([], [])
        star1_point.set_color('black')
        star1_point.set_markersize(25)

    if frame < len(t_wave):
        wave_line.set_data(t_wave[:frame], hplus[:frame])

    return [trail1_line, trail2_line, star1_point, star2_point, wave_line]



frames = max(len(x_orb), len(hplus)) + trail_length
ani = FuncAnimation(fig, update, frames=frames, init_func=init, interval=1, blit=True, repeat = False)

plt.tight_layout()
plt.show()
