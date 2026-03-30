import numpy as np
import matplotlib.pyplot as plt


def plot_ellipsoids(centers, directions, a, c, rve_side_len=1.0):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Create base sphere coordinates
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x = a * np.outer(np.cos(u), np.sin(v))
    y = a * np.outer(np.sin(u), np.sin(v))
    z = c * np.outer(np.ones_like(u), np.cos(v))

    for center, target_dir in zip(centers, directions):
        target_dir = target_dir / np.linalg.norm(target_dir)

        z_axis = np.array([0, 0, 1])
        v_cross = np.cross(z_axis, target_dir)
        s = np.linalg.norm(v_cross)
        cos_c = np.dot(z_axis, target_dir)

        if s == 0:
            R = np.eye(3) if cos_c > 0 else -np.eye(3)
        else:
            vx = np.array(
                [
                    [0, -v_cross[2], v_cross[1]],
                    [v_cross[2], 0, -v_cross[0]],
                    [-v_cross[1], v_cross[0], 0],
                ]
            )
            R = np.eye(3) + vx + (vx @ vx) * ((1 - cos_c) / (s**2))
        points = np.stack([x.flatten(), y.flatten(), z.flatten()])
        rotated_points = R @ points

        x_p = rotated_points[0, :].reshape(x.shape) + center[0]
        y_p = rotated_points[1, :].reshape(y.shape) + center[1]
        z_p = rotated_points[2, :].reshape(z.shape) + center[2]

        ax.plot_surface(x_p, y_p, z_p, alpha=1.0, color="blue", edgecolors=None, linewidth=None)
        start = center + c*target_dir
        end = center + a*target_dir
        ax.plot((start[0], end[0]), (start[1], end[1]),(start[2], end[2]), color="red")
    ax.set_box_aspect((1,1,1))
    ax.set_xlim(0.0, rve_side_len)
    ax.set_ylim(0.0, rve_side_len)
    ax.set_zlim(0.0, rve_side_len)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()


def main():
    n = 10
    centers = np.random.rand(n, 3)
    directions = np.random.standard_normal((n, 3))  # Random orientations for c-axis
    plot_ellipsoids(centers, directions, a=0.05, c=0.0125)


if __name__ == "__main__":
    main()
