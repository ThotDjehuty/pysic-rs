"""
pysicrs — Python bindings for pysic-rs.

High-performance mathematical physics engine with Rust core.
"""

import pysicrs._core as _core

# Dynamically export everything from _core (raw *_py suffixed bindings)
for _name in dir(_core):
    if not _name.startswith("_"):
        globals()[_name] = getattr(_core, _name)

# Expose clean, unsuffixed aliases (gamma -> gamma_py, rk4_solve -> rk4_solve_py, ...)
# so the documented API reads naturally without the internal `_py` suffix.
for _name in dir(_core):
    if _name.endswith("_py") and not _name.startswith("_"):
        _alias = _name[:-3]
        globals().setdefault(_alias, getattr(_core, _name))

__version__ = "0.1.0"
__author__ = "HFThot Research Lab <contact@hfthot-lab.eu>"
