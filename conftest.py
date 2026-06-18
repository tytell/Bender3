"""pytest configuration: add templates/ to sys.path so config modules are importable by tests."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'templates'))
