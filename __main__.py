"""Allow running plasrisk as `python -m plasrisk`."""
from plasrisk.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
