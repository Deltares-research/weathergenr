"""Run the portable board CLI from the activated project skill."""
import runpy
import sys
from pathlib import Path

skill = Path(__file__).resolve().parents[2] / ".agents/skills/todo-board"
if not (skill / "todoboard/__main__.py").is_file():
    skill = Path.home() / ".agents/skills/todo-board"
if not (skill / "todoboard/__main__.py").is_file():
    raise SystemExit(
        f"Board CLI not found at {skill}. Activate the Brain agent system:\n"
        '  python "$env:BRAIN_INFRA_ROOT/brain.py" refresh --agent-system --project .\n'
        "See dev/README.md#agent-activation."
    )
sys.path.insert(0, str(skill))
runpy.run_module("todoboard", run_name="__main__")
