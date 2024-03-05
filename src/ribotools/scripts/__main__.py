import shlex
import subprocess
import sys
from importlib.resources import files


def get_sample_table():
    _run("get-sample-table")


def run_dea():
    _run("run-dea")


def run_tea():
    _run("run-tea")


def _get_resource(script):
    package = "ribotools"
    resource = "scripts"
    return files(package).joinpath(resource).joinpath(script)


def _run(script):
    executable = _get_resource(script)
    cmd = f"{executable} {' '.join(sys.argv[1:])}"
    safe = shlex.split(cmd)
    subprocess.run(safe, check=True, stderr=sys.stderr, stdout=sys.stdout, text=True)
