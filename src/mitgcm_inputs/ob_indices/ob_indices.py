import logging
from collections.abc import Mapping
from collections.abc import Sequence
from enum import Enum
from itertools import product as cart_prod
from typing import Any

import numpy as np
from bitsea.commons.mask import MaskWithRivers
from bitsea.components.component_mask_2d import ComponentMask2D

LOGGER = logging.getLogger(__name__)


class Side(Enum):
    NORTH = "N"
    SOUTH = "S"
    WEST = "W"
    EAST = "E"

    def __str__(self):
        if self == Side.NORTH:
            return "north"
        elif self == Side.SOUTH:
            return "south"
        elif self == Side.WEST:
            return "west"
        elif self == Side.EAST:
            return "east"
        else:
            raise ValueError(f"Invalid side: {self}")


def cut_at_side(mask: np.ndarray, side: Side, amount: int):
    if side == Side.NORTH:
        return mask[:, -amount:, :]
    elif side == Side.SOUTH:
        return mask[:, :amount, :]
    elif side == Side.WEST:
        return mask[:, :, :amount]
    elif side == Side.EAST:
        return mask[:, :, -amount:]
    else:
        raise ValueError(f"Invalid side: {side}")


def rewrite_in_mitgcm_format(v: np.ndarray) -> list[tuple[int, int]]:
    """Transforms an array of integers into MITgcm's condensed format.

    MITgcm uses a specific run-length encoding format for ob_indices to
    compress arrays containing consecutive repeated values. Each sequence of
    identical values is represented as a tuple where the first element is the
    count (number of repetitions) and the second element is the value itself.

    For example, the array [5, 5, 5, 3, 3, 7] would be encoded as
    [(3, 5), (2, 3), (1, 7)].

    Args:
        v: A one-dimensional array of integers to be encoded.

    Returns:
        A list of tuples, where each tuple contains (count, value). The count
        represents how many consecutive times the value appears in the input
        array. Returns an empty list if the input array is empty.

    Example:
        >>> rewrite_in_mitgcm_format(np.array([1, 1, 1, 2, 2]))
        [(3, 1), (2, 2)]
    """
    if len(v) == 0:
        return []

    result = []
    current_value = v[0]
    count = 1

    for i in range(1, len(v)):
        if v[i] == current_value:
            count += 1
        else:
            result.append((count, int(current_value)))
            current_value = v[i]
            count = 1

    # Append the last group
    result.append((count, int(current_value)))

    return result


def generate_ob_indices(
    river_mask: MaskWithRivers,
    sponge_extent: int = 16,
    rivers_positions: Sequence[Mapping[str, Any]] | None = None,
):
    if rivers_positions is None:
        rivers_positions = []

    ob_indices = ""
    ob_sponge = ""

    # If sponge_extent is 0, we take the first cell anyway
    cut_cells = max(sponge_extent, 1)
    for side in Side:
        LOGGER.debug("Computing ob_indices for %s", side)

        # Here we check which index will be used for the sponge and which index
        # instead chooses the cell on the domain
        if side in (Side.NORTH, Side.SOUTH):
            sponge_index = 0
        else:
            sponge_index = 1
        domain_index = 1 - sponge_index

        # Here we check which is the index that must be used to report that a
        # cell is open on the boundary. This is the beginning or the end
        # of the domain (according to the side). If it is the beginning, it is
        # enough to put a 1, but if it is the end, we need to put the number of
        # cells on the other side of the domain.
        if side in (Side.WEST, Side.SOUTH):
            open_on_boundary = 1
        else:
            # Here we add a +1 because we need to skip the depth index
            open_on_boundary = river_mask.shape[domain_index + 1]
        # This is the value that we use to report that a cell is closed on the
        # boundary
        closed_cell = 0

        LOGGER.debug("Cutting %i cells from the %s", cut_cells, side)
        current_side = cut_at_side(river_mask.get_sea_cells(), side, cut_cells)

        # We only take into account the surface
        current_side_2d = current_side[0, :, :]

        assert current_side_2d.shape[sponge_index] == cut_cells

        sponge_cells = np.full(
            (current_side_2d.shape[domain_index]), closed_cell, dtype=int
        )
        if not np.any(current_side_2d):
            LOGGER.debug("No sea cells found on %s", side)
        else:
            component_mask = ComponentMask2D(current_side_2d)
            n_components = component_mask.n_components
            LOGGER.debug("%s has %i components", side, n_components)
            for i in range(n_components):
                LOGGER.debug("Processing component %i", i)
                component_cells = component_mask.get_component(i)

                # We need to check if this component has any cell on the
                # boundary. If this is not the case, we can skip it.
                # First, we need to define a slice that will select the
                # boundary. So we need to have a slice that takes all the
                # values on the axis that represent the domain and a fixed
                # value on the axis that represents the sponge. This fixed
                # value can be 0 or -1, depending on if the boundary is the
                # maximum or the minimum value of the current axis
                boundary_index = 0 if open_on_boundary == 1 else -1
                boundary_layer = [slice(None), slice(None)]
                boundary_layer[sponge_index] = boundary_index
                boundary_layer = tuple(boundary_layer)

                if not np.any(component_cells[boundary_layer]):
                    LOGGER.debug(
                        "Component ignored because it does not have any cell "
                        "on the boundary"
                    )
                    continue

                points_to_open = component_cells.sum(axis=sponge_index).astype(
                    bool
                )
                sponge_cells[points_to_open] = open_on_boundary

        LOGGER.debug(
            "%i boundary cells will have the sponge enabled on %s",
            (sponge_cells != closed_cell).sum(),
            side,
        )

        open_bc_cells = np.copy(sponge_cells)
        for river in rivers_positions:
            if river["model"] != "stem_flux":
                continue
            if river["side"] != side.value:
                continue
            river_sources = tuple(
                cart_prod(
                    river["latitude_indices"], river["longitude_indices"]
                )
            )
            LOGGER.debug(
                "River %s on %s has the following sources: %s",
                river["name"],
                side,
                river_sources,
            )
            for cell in river_sources:
                open_bc_cells[cell[domain_index]] = cell[sponge_index] + 1

        ob_indices += f"OB_I{side}="
        for i in rewrite_in_mitgcm_format(open_bc_cells):
            ob_indices += f"{i[0]}*{i[1]},"
        ob_indices += "\n"

        ob_sponge += f"nudgOB_I{side}="
        for i in rewrite_in_mitgcm_format(sponge_cells):
            ob_sponge += f"{i[0]}*{i[1]},"
        ob_sponge += "\n"

    return ob_indices, ob_sponge
