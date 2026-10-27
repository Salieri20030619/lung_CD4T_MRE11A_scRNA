from pathlib import Path
import os,signal
p=Path('/proc/2899665/cmdline')
if p.exists() and b'/three_cohort_four_axis_20261017_v16/code/intersect.py' in p.read_bytes():os.kill(2899665,signal.SIGTERM);print('Stopped own integration process for I/O cache optimization')