"""Minimal dependency-free notebook builder + executor.

Builds an executed Jupyter nbformat-4 notebook from a list of
(markdown | code) cells:

    from nblite import build
    build("05_dirac_spinors.ipynb", [
        ("md",  "# Title\\n prose"),
        ("code", "import numpy as np\\nprint(np.arange(3))"),
    ])

Code cells are executed IN ORDER in one shared namespace, exactly as a
Jupyter kernel would; stdout is captured into stream outputs and any
exception is captured as an error output and then raised so that the
build fails loudly rather than shipping a half-executed notebook.
"""
from __future__ import annotations

import contextlib
import io
import json
import platform
import traceback


def build(out_path: str, cells: list[tuple[str, str]]) -> None:
    nb = {
        "cells": [],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python",
                           "name": "python3"},
            "language_info": {"name": "python",
                              "version": platform.python_version()},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    ns: dict = {"__name__": "__main__"}
    count = 0
    for kind, src in cells:
        cell = {"cell_type": kind, "metadata": {}, "source": src}
        if kind == "code":
            count += 1
            buf = io.StringIO()
            outs = []
            try:
                with contextlib.redirect_stdout(buf):
                    exec(compile(src, "<cell %d>" % count, "exec"), ns)
            except Exception:
                cell["outputs"] = [{
                    "output_type": "error",
                    "ename": type(exc) if (exc := None) is None else "",
                    "evalue": "", "traceback": [],
                }]
                raise
            if buf.getvalue():
                outs.append({"output_type": "stream", "name": "stdout",
                             "text": buf.getvalue().splitlines(True)})
            cell["outputs"] = outs
            cell["execution_count"] = count
        nb["cells"].append(cell)
    with open(out_path, "w") as f:
        json.dump(nb, f, indent=1)
        f.write("\n")
    print(f"wrote {out_path} ({len(nb['cells'])} cells, {count} executed)")
