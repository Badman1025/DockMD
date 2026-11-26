import json
import os
import subprocess
import sys
from datetime import datetime

import config as cfg

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
STATE_FILE = os.path.join(PROJECT_DIR, cfg.RESULT_DIR, "state.json")


def run(cmd, cwd=None, stdin_str=None, check=True, quiet=False):
    printable = " ".join(cmd) if isinstance(cmd, list) else cmd
    print(f"\n>>> {printable}")
    try:
        return subprocess.run(
            cmd, cwd=cwd, check=check, text=True,
            input=stdin_str,
            stdout=subprocess.DEVNULL if quiet else None,
        )
    except subprocess.CalledProcessError as e:
        print(f"!! Command failed: {printable}")
        if e.returncode:
            print(f"   Exit code: {e.returncode}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"!! Executable not found: {cmd[0] if isinstance(cmd, list) else cmd}")
        print("   Please check the TOOLS paths in config.py.")
        sys.exit(1)


def load_state():
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE, encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_state(**kwargs):
    state = load_state()
    state.update(kwargs)
    state["updated_at"] = datetime.now().isoformat(timespec="seconds")
    with open(STATE_FILE, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2)
    print(f"[state] Updated {STATE_FILE}: "
          + ", ".join(f"{k}={v}" for k, v in kwargs.items()))
    return state


def clear_state_keys(keys):
    state = load_state()
    removed = [k for k in keys if k in state]
    if not removed:
        return
    for k in removed:
        state.pop(k)
    state["updated_at"] = datetime.now().isoformat(timespec="seconds")
    with open(STATE_FILE, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2)
    print(f"[state] Cleared stale keys: {', '.join(removed)}")


def get_path(state, key, fallback_finder=None, what=""):
    p = state.get(key)
    if p and os.path.isdir(p):
        return p
    if fallback_finder:
        p2 = fallback_finder()
        if p2:
            print(f"[state] {what}: no valid record in state, falling back to {p2}")
            return p2
    sys.exit(f"!! {what} not found. Run the upstream step first, "
             f"or make sure the '{key}' directory exists.")
