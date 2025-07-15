import datascout as ds
import numpy as np
import matplotlib.pyplot as plt
import time

# --- Load and reshape data ---
bpm_id = 0  # choose BPM source
d = ds.parquet_to_dict('2025.04.03.11.01.21.904741.parquet')
N_bunches = d['BPLOFSBA5/GetCapData']['value']['nbOfCapBunches']
N_turns = d['BPLOFSBA5/GetCapData']['value']['nbOfCapTurns']
v_ids = d['BPLOFSBA5/GetCapData']['value']['verBunchId'][bpm_id]
v_pos = d['BPLOFSBA5/GetCapData']['value']['verPosition'][bpm_id]
h_ids = d['BPLOFSBA5/GetCapData']['value']['horBunchId'][bpm_id]
h_pos = d['BPLOFSBA5/GetCapData']['value']['horPosition'][bpm_id]

v_ids = np.reshape(v_ids, (N_turns, N_bunches))
h_ids = np.reshape(h_ids, (N_turns, N_bunches))

v_pos = np.reshape(v_pos, (N_turns, N_bunches))
h_pos = np.reshape(h_pos, (N_turns, N_bunches))

# --- Animation parameters ---
n_turns_per_frame = 1        # how many turns to show per frame
frame_rate = N_turns/20              # frames per second (sleep time = 1/frame_rate)
pause_time = 1.0 / frame_rate
trace_every = 10  # trace every n turns

# --- Set up the plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))
line1, = ax1.plot([], [], 'b')
line2, = ax2.plot([], [], 'r')

ax1.set_xlim(np.min(v_ids), np.max(v_ids))
ax1.set_ylim(np.min(v_pos), np.max(v_pos))
ax1.set_ylabel("Vertical Position [a.u.]")
ax1.set_title("Vertical")

ax2.set_xlim(np.min(h_ids), np.max(h_ids))
ax2.set_ylim(np.min(h_pos), np.max(h_pos))
ax2.set_ylabel("Horizontal Position [a.u.]")
ax2.set_xlabel("Bunch ID")
ax2.set_title("Horizontal")

plt.ion()  # Turn on interactive mode
plt.show()

# --- Main animation loop ---
for turn in range(0, N_turns, n_turns_per_frame):
    line1.set_data(v_ids[turn], v_pos[turn])
    line2.set_data(h_ids[turn], h_pos[turn])

    if turn % trace_every == 0:
        ax1.plot(v_ids[turn], v_pos[turn], 'b', lw=0.7, alpha=0.2,)
        ax2.plot(h_ids[turn], h_pos[turn], 'r', lw=0.7, alpha=0.2,)

    fig.suptitle(f"Bunch Position Evolution, turn = {turn}")
    fig.canvas.draw()
    fig.canvas.flush_events()
    time.sleep(pause_time)
