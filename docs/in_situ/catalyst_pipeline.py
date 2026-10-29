"""Minimal ParaView 5.13 Catalyst 2 smoke-test pipeline for Nek5000."""

import os

from paraview import catalyst
from paraview.simple import SaveData


OUTPUT_DIRECTORY = os.environ.get("NEK_CATALYST_OUTPUT", "catalyst-output")
FREQUENCY = int(os.environ.get("NEK_CATALYST_FREQUENCY", "10"))


def catalyst_execute(info):
    """Write one parallel VTU data set every FREQUENCY time steps."""
    if info.cycle % FREQUENCY:
        return

    mesh = catalyst.get_producer(info, "mesh")
    filename = os.path.join(OUTPUT_DIRECTORY, f"nek_{info.cycle:08d}.pvtu")
    SaveData(filename, proxy=mesh)
