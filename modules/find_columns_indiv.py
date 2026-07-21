import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

def find_columns_indiv(norm_tab, indiv):
    """
    Find columns containing individual identifiers like 'F' or 'M'.
    Returns column indices.
    """
    col_labels = norm_tab.columns.tolist()
    matches = [i for i, col in enumerate(col_labels) if indiv in col]
    if len(matches) < 2:
        print(f"Warning: Could not find two matching columns for '{indiv}'")
    return matches, len(matches)