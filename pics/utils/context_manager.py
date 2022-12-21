import sys


class temp_path:
    """Add a temporary system path for loading or importing files or modules"""

    def __init__(self, path):
        self.path = path

    def __enter__(self):
        sys.path.insert(0, self.path)

    def __exit__(self, exec_type, exec_value, traceback):
        try:
            sys.path.remove(self.path)

        except ValueError:
            pass
