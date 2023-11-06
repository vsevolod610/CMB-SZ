#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Path managment: do "import paths" in file.py
"""

# Imports: built-in
import os
import sys


# Find current dir from scr dir
current = os.path.dirname(__file__)
src = os.path.join(current, '..', '..')

# PATH is a Environment variable, add scr in PATH
sys.path.append(src)

if __name__ == '__main__':
    print(f"{current = }")
    print(f"{src = }")
    print(f"PATH = ", *sys.path, sep='\n')
