import atexit
import hashlib
import os
import platform
import subprocess
import sys
import traceback
from datetime import datetime
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path

PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "seaborn", "plotly", "scikit-learn", "PyYAML"]


class RunLog:
    """
    Plain-text record of one pipeline run, saved as
    <output_dir>/run_log_<YYYYmmdd_HHMMSS>.txt (one file per run, never appended to).

    The log is closed automatically when Python exits: COMPLETED if ``finish()``
    was called, FAILED with the traceback or ``stop()`` message if an error
    stopped the run, INCOMPLETE otherwise (e.g. a bare ``sys.exit``).

    Parameters
    ----------
    output_dir : str or Path
        Run output directory (must exist).
    config_path : str or Path
        Config file the run was started with; copied verbatim into the log.
    scalibr_version : str
        Installed scalibr version, or 'unknown'.
    """

    def __init__(self, output_dir, config_path, scalibr_version):
        self.start = datetime.now()
        self.path = Path(output_dir) / f"run_log_{self.start:%Y%m%d_%H%M%S}.txt"
        self.config_path = Path(config_path)
        self._finished = False
        self._error = None

        self.write("sCALIBR run log")
        self.write("=" * 60)
        self.field("Started", f"{self.start:%Y-%m-%d %H:%M:%S}")
        self.field("Command", " ".join([os.path.basename(sys.executable)] + sys.argv))
        self.field("Working directory", os.getcwd())
        self.field("Config file", self.config_path.resolve())
        self.field("Output directory", Path(output_dir).resolve())

        self.section("Code")
        self.field("scalibr version", scalibr_version)
        self.field("git", _git_state(Path(__file__).resolve().parent.parent))
        self.field("Python", f"{platform.python_version()} ({sys.executable})")
        self.field("Packages", ", ".join(f"{p} {_pkg_version(p)}" for p in PACKAGES))

        previous_hook = sys.excepthook

        def _record_error(exc_type, exc, tb):
            self._error = "".join(traceback.format_exception(exc_type, exc, tb))
            previous_hook(exc_type, exc, tb)

        sys.excepthook = _record_error
        atexit.register(self._close)

    def write(self, text=""):
        with open(self.path, "a") as fh:
            fh.write(f"{text}\n")

    def section(self, title):
        self.write()
        self.write(f"[{title}]")

    def field(self, key, value):
        self.write(f"{key + ':':<28}{value}")

    def input_file(self, path, table):
        """Record the count file: location, checksum and size."""
        self.section("Input")
        self.field("Counts file", Path(path).resolve())
        self.field("MD5", _md5(path))
        self.field("Size", f"{table.shape[0]} gRNAs x {table.shape[1]} columns")
        self.field("Sample columns", ", ".join(map(str, table.columns[3:])))

    def config_text(self):
        """Copy the config file into the log exactly as written."""
        self.section(f"Config file as written ({self.config_path.name})")
        try:
            self.write(self.config_path.read_text().rstrip())
        except OSError as e:
            self.write(f"(could not read: {e})")

    def finish(self):
        self._finished = True

    def stop(self, message):
        """Record why the run is stopping, then exit with that message."""
        self._error = message
        sys.exit(message)

    def _close(self):
        end = datetime.now()
        if self._finished:
            status = "COMPLETED"
        elif self._error:
            status = "FAILED"
        else:
            status = "INCOMPLETE (stopped before the end)"
        self.section("Finished")
        self.field("Ended", f"{end:%Y-%m-%d %H:%M:%S}")
        self.field("Duration", str(end - self.start).split(".")[0])
        self.field("Status", status)
        if self._error:
            self.write(self._error.rstrip())


def _pkg_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def _git_state(repo_dir):
    """Short commit, branch and whether tracked files had uncommitted changes."""
    def git(*args):
        return subprocess.run(["git", "-C", str(repo_dir), *args], capture_output=True,
                              text=True, timeout=10, check=True).stdout.strip()
    try:
        commit = git("rev-parse", "--short", "HEAD")
        branch = git("rev-parse", "--abbrev-ref", "HEAD")
        changed = git("status", "--porcelain", "--untracked-files=no").splitlines()
    except (OSError, subprocess.SubprocessError):
        return "not available (not a git checkout)"
    dirty = f"uncommitted changes in {len(changed)} file(s)" if changed else "no uncommitted changes"
    return f"{commit} on {branch}, {dirty}"


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
