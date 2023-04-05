from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from pyvrp import ProblemData


def plot_demands(
    data: ProblemData,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
):
    """
    Plots demands for clients, as vertical bars sorted by demand.

    Parameters
    ----------
    data
        Data instance.
    title, optional
        Title to add to the plot.
    ax, optional
        Axes object to draw the plot on. One will be created if not provided.
    """
    if not ax:
        _, ax = plt.subplots()

    dim = data.num_clients + 1
    # Exclude depot
    demand = np.array([data.client(client).demand for client in range(1, dim)])
    demand = np.sort(demand)

    ax.bar(np.arange(1, dim), demand)

    if title is None:
        vehicle_caps = np.array(
            [data.route(i).vehicle_capacity for i in range(data.num_vehicles)]
        )
        title = (
            f"Demands (avg. cap = {vehicle_caps.mean()}, "
            + f"{vehicle_caps.mean() / demand.mean():.2f} stops/route)"
        )

    ax.set_title(title)
    ax.set_xlabel("Client (sorted by demand)")
    ax.set_ylabel("Demand")
