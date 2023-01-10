from BlackLock_mod import *
import sys
from Bio import SeqIO
import gzip
import shutil
import re
import math

if len(sys.argv) != 2:
    help(BLhelp)
    exit()

f = sys.argv[1]
print(f)
Tata2(f, minlen = 0, maxlen = 'a', minPhred=0, barcode='a', totalb='a', seq_search='')
