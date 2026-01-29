import os
import csv
import matplotlib.pyplot as plt
import glob


def _resolve_body_csv(output_dir, ts_file):
    """
    Resolve the path to temperature_body_*.csv based on the provided ts_file.
    - If ts_file is an existing file path, use it directly.
    - Otherwise, treat ts_file as the timestep/tag and compose the path.
    - If ts_file looks like a filename (e.g., temperature_body_100.0.csv),
      extract the tag (100.0) and compose the path under output_dir/per_body/.
    """
    # If given a concrete, existing file path, return it as-is
    if ts_file and os.path.isfile(ts_file):
        return ts_file

    tag = ts_file
    if ts_file:
        base = os.path.basename(ts_file)
        if base.endswith('.csv'):
            base_no_ext = base[:-4]
            if base_no_ext.startswith('temperature_body_'):
                tag = base_no_ext[len('temperature_body_'):]

    # New layout (breaking change): body CSVs live under output/RC/per_body.
    # Here, output_dir is expected to be the RC directory.
    cand_new = os.path.join(output_dir, 'per_body', f'temperature_body_{tag}.csv')
    if os.path.exists(cand_new):
        return cand_new

    # Legacy fallbacks (best-effort)
    cand_legacy_per_body = os.path.join(output_dir, 'output', 'per_body', f'temperature_body_{tag}.csv')
    if os.path.exists(cand_legacy_per_body):
        return cand_legacy_per_body
    cand_legacy_flat = os.path.join(output_dir, 'output', f'temperature_body_{tag}.csv')
    return cand_legacy_flat


def plot_final_body_temperatures(output_dir, ts_file, time_s, dt):
    body_csv = _resolve_body_csv(output_dir, ts_file)
    if not os.path.exists(body_csv):
        print(f"[plot] body CSV not found at expected path: {body_csv}")
        # fallback: pick the latest temperature_body_*.csv (new layout first)
        cand = sorted(glob.glob(os.path.join(output_dir, 'per_body', 'temperature_body_*.csv')))
        if not cand:
            cand = sorted(glob.glob(os.path.join(output_dir, 'output', 'per_body', 'temperature_body_*.csv')))
        if not cand:
            cand = sorted(glob.glob(os.path.join(output_dir, 'output', 'temperature_body_*.csv')))
        if not cand:
            print(f"[plot] no temperature_body_*.csv files found under {output_dir}")
            return
        body_csv = cand[-1]
        print(f"[plot] using fallback body CSV: {body_csv}")

    # determine column index for requested time
    if dt <= 0:
        col_idx = -1
    else:
        col_idx = int(time_s/float(dt))

    names = []
    values = []
    with open(body_csv, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            name = row[0]
            temps = [t for t in row[1:] if t != '']
            if not temps:
                continue
            # clamp index
            idx = col_idx if 0 <= col_idx < len(temps) else len(temps)-1
            try:
                val_kelvin = float(temps[idx])
            except ValueError:
                continue
            names.append(name)
            values.append(val_kelvin - 273.15)

    if not names:
        print('[plot] parsed body CSV but found no valid rows to plot')
        return

    plots_dir = os.path.join(output_dir, 'plots')
    if not os.path.exists(plots_dir):
        try:
            os.makedirs(plots_dir)
        except Exception as e:
            print(f"[plot] failed to create plots directory {plots_dir}: {e}")
            return

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(range(len(names)), values)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha='right')
    ax.set_ylabel('Temperature (°C)')
    ax.set_title(f'Final Body Temperatures at t={time_s}s')
    fig.tight_layout()
    out_path = os.path.join(plots_dir, f'final_body_temperatures_t{time_s}.png')
    fig.savefig(out_path, dpi=300)
    print(f"[plot] saved final body temperature bar chart to: {out_path}")
    plt.close(fig)


def plot_body_temperatures_over_time(output_dir, ts_file, total_duration, dt):
    body_csv = _resolve_body_csv(output_dir, ts_file)
    if not os.path.exists(body_csv):
        print(f"[plot] body CSV not found at expected path: {body_csv}")
        # fallback to latest
        cand = sorted(glob.glob(os.path.join(output_dir, 'per_body', 'temperature_body_*.csv')))
        if not cand:
            cand = sorted(glob.glob(os.path.join(output_dir, 'output', 'per_body', 'temperature_body_*.csv')))
        if not cand:
            cand = sorted(glob.glob(os.path.join(output_dir, 'output', 'temperature_body_*.csv')))
        if not cand:
            print(f"[plot] no temperature_body_*.csv files found under {output_dir}")
            return
        body_csv = cand[-1]
        print(f"[plot] using fallback body CSV: {body_csv}")

    plots_dir = os.path.join(output_dir, 'plots')
    if not os.path.exists(plots_dir):
        try:
            os.makedirs(plots_dir)
        except Exception as e:
            print(f"[plot] failed to create plots directory {plots_dir}: {e}")
            return

    fig, ax = plt.subplots(figsize=(10, 6))
    
    with open(body_csv, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            name = row[0]
            try:
                temps_k = [float(t) for t in row[1:] if t]
                temps_c = [tk - 273.15 for tk in temps_k]
            except ValueError:
                continue
            
            time_points = [i * dt for i in range(len(temps_c))]
            ax.plot(time_points, temps_c, label=name)

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Temperature (°C)')
    ax.set_title('Body Temperatures Over Time')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    out_path = os.path.join(plots_dir, 'body_temperatures_over_time.png')
    fig.savefig(out_path, dpi=300)
    print(f"[plot] saved body temperature line plot to: {out_path}")
    plt.close(fig)


