
# Adapted from sunpy. May need to modify base on how eispac uses tags/releases
# NOTE: First try _dev.scm_version if it exists and setuptools_scm is installed
# This file is not included in wheels/tarballs, so otherwise it will
# fall back on the generated _version module.
try:
    try:
        from ._dev.scm_version import version as full_version
    except ImportError:
        from ._version import version as full_version
except Exception:
    import warnings

    warnings.warn(
        f'could not determine {__name__.split(".")[0]} package version; '
        f'this may indicate a broken installation. Defaulting to package '
        f'metadata (if available).')
    del warnings

    try:
        from importlib.metadata import version as meta_ver
        full_version = meta_ver('eispac')
    except:
        full_version = '0.0.0'

from packaging.version import parse as _parse

version = _parse(full_version).base_version # Clean version num for PyPi
_version = _parse(full_version)
major, minor, bugfix = [*_version.release, 0][:3]
release = not _version.is_devrelease
