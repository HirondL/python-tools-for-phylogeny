# python-tools-for-phylogeny

The two .py files contain Python tools for manipulating and curating phylogeny-related materials: sequences, alignments, and trees.
For example, one can:
- Retrieve sequences from accession numbers
- Reformat sequence names
- Filter sequences by pattern
- Check and remove redundancy
- Encode/decode full names to use .phy alignments
    This is useful when third-party computing resources require the .phy format, which trims sequence names to 30 characters.
- Sort sequences based on tree structure

The .ipynb file shows a mock workflow so one can directly copy these blocks.

Installation:
```bash
git clone https://github.com/HirondL/python-tools-for-phylogeny
cd python-tools-for-phylogeny
pip install -r requirements.txt
```
