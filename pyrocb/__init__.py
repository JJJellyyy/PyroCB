"""pyroCb x FireAtlas pipeline.

Modules
-------
config      paths, constants, year resolution, the caching helper
io          readers for the four products + the geometry validity pipeline
timebins    half-day step handling and the 'snapshot before onset' rule
matching    pyroCb event -> fire (mergid)
bands       cumulative-union growth bands split at every pyroCb
active      active-burn polygon at onset, from NFP pixel footprints
variables   the per-event attribute table
"""
from . import config, io, timebins  # noqa: F401

__all__ = ["config", "io", "timebins", "matching", "bands", "variables",
           "active", "landcover", "firesim", "pipeline"]
