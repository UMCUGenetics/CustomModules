#! /usr/bin/env python
# Import statements, alphabetic order of main package.
from csv import writer
from errno import ENOENT as errno_ENOENT
from os import strerror as os_strerror
from pathlib import Path

def non_empty_existing_path(file_or_dir):
    input_path = Path(file_or_dir)
    if not input_path.is_file() and not input_path.is_dir():
        raise FileNotFoundError(errno_ENOENT, os_strerror(errno_ENOENT), file_or_dir)
    elif not input_path.is_dir() and not input_path.stat().st_size:
        raise OSError(f"File {file_or_dir} is empty.")
    elif input_path.is_dir() and file_or_dir[::-1][0] != '/':
        return f"{file_or_dir}/"
    return file_or_dir