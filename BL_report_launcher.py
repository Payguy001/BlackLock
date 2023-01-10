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

fastqFile = sys.argv[1]
"""
y = fastqFile.split('.')
ly = y[len(y)-1]
if ly != 'fastq':
    if ly != 'gz':
        help(BLhelp)
    else:
        if y[len(y)-2] == 'fastq':
            newfq = fastqFile[:-3]
        else:
            newfq = fastqFile[:-3]+'.fastq'
        with gzip.open(fastqFile, 'rb') as f_in:
            with open(newfq, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
else:
    newfq = fastqFile
"""
data_ready = file2fq(fastqFile)
print('[STATUS]: File is prepared in .fastq extension form')
gi = fq2dict(data_ready)
print('[STATUS]: Dictionary of the data is created')

#ค่า Heading Row ของ Table
Heading = ["Barcode", "Number of seq", "Average length", "Average QScore", "Density (%)"]
Row = [["001", "002"], [1000, 500], [500, 250], [100, 200], [50, 40], [0.5,0.4]]
lendict ,qsc = createdict(gi)
Row = inputtable(lendict, qsc)

def test():
    #สร้าง script
    table4script = GenTableHTML("Table4",Heading,Row)
    lenplotlist, qscplotlist = GeneratePlotObject(gi)
    Histo1script = GenerateHistogramFromObject(Name = "Histo1", Header = "Length VS count", xlabel = "Length", BarcodeList = lenplotlist)
    Histo2script = GenerateHistogramFromObject(Name = "Histo2", Header = "QScore VS count", xlabel = "QScore", BarcodeList = qscplotlist)
    filename = str(fastqFile)

    #สร้าง HTML
    Gen=GenReadableHTML("visualize4",table4script)
    re1=Gen.replace("histo1",Histo1script)
    re1=re1.replace("histo2",Histo2script)
    re1=re1.replace("filename",filename)
    print('\n=============')
    a = str(input("[Input] Please enter filename: "))
    a= a +".html"
    with open(a, "w") as output:
        output.write(re1)
    print("[STATUS]: File "+a+" is completely prepared")
    print('=============\n')
#if (__name__ == "__main__"):
test()