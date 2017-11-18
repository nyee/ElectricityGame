"""
Run as cleariPython.py < my_notebook.ipynb > my_notebook_stripped.ipynb

"""

def strip_output(nb):
    for ws in nb.worksheets:
        for cell in ws.cells:
            if hasattr(cell, "outputs"):
                cell.outputs = []
            if hasattr(cell, "prompt_number"):
                del cell["prompt_number"]


if __name__ == "__main__":
    from sys import stdin, stdout
    # from IPython.nbformat.current import read, write
    from nbformat.current import read, write

    # stdin = "/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/sdata134k_small_polycyclic_v3_pyMaster_dbMaster.ipynb"
    # stdout = "/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/sdata134k_small_polycyclic_v3_pyMaster_dbMaster_cleared.ipynb"

    nb = read(stdin, "ipynb")
    strip_output(nb)
    write(nb, stdout, "ipynb")
    stdout.write("\n")
