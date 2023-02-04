from typing import Optional, Union

from pyvrp.OptimisationTarget import OptimisationTarget


class NoImprovement:
    """
    Criterion that stops if the best solution has not been improved for a fixed
    number of iterations.

    Parameters
    ----------
    max_iterations
        The maximum number of non-improving iterations.
    """

    def __init__(self, max_iterations: int):
        if max_iterations < 0:
            raise ValueError("max_iterations < 0 not understood.")

        self._max_iterations = max_iterations
        self._target: Optional[Union[int, float]] = None
        self._counter = 0

    def __call__(self, best: OptimisationTarget) -> bool:
        if self._target is None or best.cost() < self._target:
            self._target = best.cost()
            self._counter = 0
        else:
            self._counter += 1

        return self._counter >= self._max_iterations
