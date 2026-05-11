import os, sys

ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')
)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


from class1 import matching_diagram as match

