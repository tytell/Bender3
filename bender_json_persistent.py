"""
Strict JSON persistence for GUI autosave, templates, and debug payloads.

Converts numpy scalars/arrays to Python primitives and lists; rejects non-data
objects instead of calling ``str()`` (which caused silent corruption on restore).
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

__all__ = ['JsonPersistTypeError', 'to_json_persistent']


class JsonPersistTypeError(TypeError):
    """Value cannot be written to JSON persistence (autosave, templates, etc.)."""


def _path_join(path: str, segment: str) -> str:
    seg = str(segment)
    if not path:
        return seg
    if seg.startswith('['):
        return f'{path}{seg}'
    return f'{path}.{seg}'


def _reject(value: Any, path: str) -> None:
    raise JsonPersistTypeError(
        f'Non-JSON-persistable value at {path or "(root)"}: {type(value).__name__}'
    )


def _is_file_like(value: Any) -> bool:
    return hasattr(value, 'read') and callable(getattr(value, 'read', None))


def to_json_persistent(value: Any, *, path: str = '') -> Any:
    """
    Convert *value* to a tree of JSON primitives (null, bool, int, float, str, list, dict).

    ``np.ndarray`` → nested lists; numpy scalars → Python numbers. Raises
    :class:`JsonPersistTypeError` for callables, file handles, and other opaque objects.
    """
    if value is None:
        return None

    if isinstance(value, bool):
        return value

    if isinstance(value, str):
        return value

    if isinstance(value, int) and not isinstance(value, bool):
        return int(value)

    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return float(value)

    if isinstance(value, (np.bool_,)):
        return bool(value)

    if isinstance(value, (np.integer,)):
        return int(value.item())

    if isinstance(value, (np.floating,)):
        fv = float(value.item())
        if not math.isfinite(fv):
            return None
        return fv

    if isinstance(value, np.ndarray):
        return to_json_persistent(value.tolist(), path=path)

    if isinstance(value, (list, tuple)):
        return [to_json_persistent(x, path=_path_join(path, f'[{i}]')) for i, x in enumerate(value)]

    if isinstance(value, dict):
        out: dict[str, Any] = {}
        for kk, vv in value.items():
            out[str(kk)] = to_json_persistent(vv, path=_path_join(path, str(kk)))
        return out

    if callable(value):
        _reject(value, path)

    if _is_file_like(value):
        _reject(value, path)

    mod = getattr(type(value), '__module__', '') or ''
    if mod not in ('builtins',) and not isinstance(value, (int, float, str, bool)):
        _reject(value, path)

    _reject(value, path)
