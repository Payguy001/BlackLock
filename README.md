# BlackLock
Visualize and Filter DNA sequence

Compatible sequence file type: `.fastq` and `.gz`
Run `BlackLock` file to start

```
python BlackLock -h
```
## Visualize DNA Sequence
Use `VisualizeSeq` to generate HTML file for visualizing data 

Example:
```
python BlackLock VisualizeSeq -f FILENAME
```
## Filter DNA Sequence
Use FilterSeq to generate filtered FASTQ file

**Example:**
```
python BlackLock FilterSeq -f FILENAME -mi 100
```
  This example generates filtered fastq file for every sequence that is at least 100 bases long

**Option argument for FilterSeq**

```
usage: BlackLock FilterSeq [-h] [-f FNAME] [-mi MINLEN] [-ma MAXLEN] [-ph MINPHRED] [-bar BARCODE] [-tb TOTALB]
                           [-seq SEQ_SEARCH]

optional arguments:
  -h, --help            show this help message and exit
  -f FNAME, --fname FNAME
                        Input Sequence filename (.fastq or .gz)
  -mi MINLEN, --minlen MINLEN
                        Provide minimum length for filtering
  -ma MAXLEN, --maxlen MAXLEN
                        Provide maximum length for filtering
  -ph MINPHRED, --minPhred MINPHRED
                        Provide minimum Phred Score for filtering
  -bar BARCODE, --barcode BARCODE
                        Provide barcode name for filtering
  -tb TOTALB, --totalb TOTALB
                        Provide total base to be read for filtering
  -seq SEQ_SEARCH, --seq_search SEQ_SEARCH
                        Provide sequence for searching for filtering
```
