from typing import overload

class PenaltyBooster:
    def __init__(self, *args, **kwargs) -> None: ...
    def __enter__(self) -> PenaltyBooster: ...
    def __exit__(
        self, type: object, value: object, traceback: object
    ) -> None: ...

class PenaltyManager:
    @overload
    def __init__(
        self, vehicle_capacity: int, params: PenaltyParams
    ) -> None: ...
    @overload
    def __init__(self, vehicle_capacity: int) -> None: ...
    def get_penalty_booster(self) -> PenaltyBooster: ...
    def load_penalty(self, load: int) -> int: ...
    def tw_penalty(self, time_warp: int) -> int: ...
    def update_capacity_penalty(self, curr_feas_pct: float) -> None: ...
    def update_time_warp_penalty(self, curr_feas_pct: float) -> None: ...

class PenaltyParams:
    def __init__(
        self,
        init_capacity_penalty: int = ...,
        init_time_warp_penalty: int = ...,
        repair_booster: int = ...,
        penalty_increase: float = ...,
        penalty_decrease: float = ...,
        target_feasible: float = ...,
    ) -> None: ...
    @property
    def init_capacity_penalty(self) -> int: ...
    @property
    def init_time_warp_penalty(self) -> int: ...
    @property
    def penalty_decrease(self) -> float: ...
    @property
    def penalty_increase(self) -> float: ...
    @property
    def repair_booster(self) -> int: ...
    @property
    def target_feasible(self) -> float: ...
