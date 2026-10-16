import re

def sample_group(name):
    """Strip a trailing _R<N> replicate tag: 'T0_F_R2' -> 'T0_F'."""
    m = re.search(r"_R\d+", name)
    return name[:m.start()] if m else name