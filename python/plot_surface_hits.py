import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "surface_hits.csv"
    out_png = sys.argv[2] if len(sys.argv) > 2 else "surface.png"
    df = pd.read_csv(path)

    if len(df) == 0:
        print("No rows in file.")
        return

    # 只保留真正打到输出平面的光子
    df = df[df["hit_surface"] == 1].copy()
    if len(df) == 0:
        print("No surface hits.")
        return

    u = df["u_m"].to_numpy(float)
    v = df["v_m"].to_numpy(float)
    w = df["weight"].to_numpy(float) * df["relative_efficiency"].to_numpy(float)

    r2 = u * u + v * v
    rms = np.sqrt(np.sum(w * r2) / np.sum(w))

    print(f"N hits = {len(df)}")
    print(f"Weighted RMS = {rms:.6e} m")

    # 图 1：spot 图
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(u * 100.0, v * 100.0, s=4, alpha=0.6)
    ax.set_xlabel("u [cm]")
    ax.set_ylabel("v [cm]")
    ax.set_title("Optical surface spot")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)

    # 图 2：按镜片编号着色的 spot
    if "mirror_id" in df.columns:
        fig2, ax2 = plt.subplots(figsize=(6, 6))
        sc = ax2.scatter(u * 100.0, v * 100.0, c=df["mirror_id"], s=4, alpha=0.7)
        ax2.set_xlabel("u [cm]")
        ax2.set_ylabel("v [cm]")
        ax2.set_title("Surface hits colored by mirror_id")
        ax2.set_aspect("equal", adjustable="box")
        ax2.grid(True, alpha=0.3)
        plt.colorbar(sc, ax=ax2, label="mirror_id")
        
    plt.savefig(out_png)
    print(f"Saved plot = {out_png}")
    if "agg" not in plt.get_backend().lower():
        plt.show()

if __name__ == "__main__":
    main()
