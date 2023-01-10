import numpy as np
import sys
import re
from Bio import SeqIO
import gzip
import shutil
import math
import os

class BarcodePlotObject:
    def __init__(self, name = "Unnamed Barcode",x = [], y = []):
        self.name = name
        self.x = x[:]
        self.y = y[:]
        
class fastq_BL:
    """ fastq_BL class contain id, name, seq, quality, barcode """
    def __init__(self, Id='', name='', seq='', quality=list(), barcode='N/A', raw=''):
        self.id = Id
        self.name = name
        self.seq = seq
        self.len = len(self.seq)
        self.qual = quality
        self.mean = sum(self.qual)/self.len
        self.bar = barcode
        self.raw = raw

def BLhelp():
    """
    To use this, in command line type
    >>>python BlackLock_temp.py [fastqfile]

    [EXAMPLE]
    >>>python BlackLock_temp.py example.reads.fastq
    
    Please feel free to modify this! 
    ---Press q to exit---

    """

def Tata(SeqFileName, minlen = 0, maxlen = 0, minPhred=0, barcode='a', totalb='a', seq_search=''):
    saving_status = saving_dialog()

    fastqFile = SeqFileName
    data_ready = file2fq(fastqFile)
    # print('[STATUS]: File is prepared in .fastq extension form')
    gi = fq2dict(data_ready)
    print('[STATUS]: Dictionary of the data is created')

    #filter
    jj = filter_BL(gi, minlen, maxlen, minPhred, barcode, totalb, seq_search)
    print('[STATUS]: Filtered dictionary created... begin reading...')
    dict_BL_read(jj)

    #Save file
    sav = input("Would you like to save filtered data? [Y/n] : ")
    if(sav != 'n'):
        while True:
            
            ty = input("Please put in filename with extention : ")
            ex = f2ext(ty)
            if ex=='folder':
                print('[WARNING]: Please add extension (.fastq,.gz)')
            else:
                break
            
        if ex=='fastq': #.fastq
            dict_BL2fq(jj, ty[:-5]+'.fastq') #เขียน filter
        if ex=='gz':
            if f2ext(ty[:-3]) != 'fastq': #.fastq.gz
                dict_BL2fq(jj, ty[:-3]+'.fastq')
                fq2gz(ty[:-3]+'.fastq') #เปลี่ยนเป็น gz
            else: #.gz
                dict_BL2fq(jj, ty[:-3])
                fq2gz(ty[:-3])
        print("[STATUS]: File is saved as "+ty)
    else:   print("[STATUS]: File is not saved")
    
    print("\n++++++++++++++++++\nBLACKLOCK THANK YOU\n++++++++++++++++++\n")
    
def Tata2(SeqFileName, minlen = 0, maxlen = 0, minPhred=0, barcode='a', totalb='a', seq_search=''):

    data_ready = file2fq(SeqFileName)
    print('[STATUS]: File is prepared in .fastq extension form')
    s = saving_dialog()
    print(type(s))
    #print('c1')
    ult_filter(SeqFileName, minlen, maxlen, minPhred, barcode, totalb, seq_search,save_file=s)

    print("\n++++++++++++++++++\nBLACKLOCK THANK YOU\n++++++++++++++++++\n")


def ult_filter(filename, minlen=0, maxlen='a', minPhred=0, barcode='a', totalb='a', seq_search='',save_file=''):
    
    ext_d = re.compile(r'barcode=\S+',re.M)
    
    print('[FILTER]: '+'Inputed Parameters:(a=all/max)')
    print('[FILTER]: '+'min seq length\t:',minlen)
    print('[FILTER]: '+'max seq length\t:',maxlen)
    print('[FILTER]: '+'min phred score\t:',minPhred)
    print('[FILTER]: '+'barcode\t\t:',barcode)
    print('[FILTER]: '+'sequence search\t:',seq_search)
    print('[FILTER]: '+'total base pair\t:',totalb)
    print('[FILTER]: '+'Begin filtering process')
    #print('c2')
    fqfile = file2fq(filename)
    out = list()
    ttr_len= 0
    ttr_n= 0
    ttr_ex = dict()
    ttr_bar = dict()
    ttr = list()
    crec=0
    kkk=1
    for record in SeqIO.parse(fqfile,'fastq'):
        #print(record)
        #Isagi[ID] = fastq_BL(record.id,record.id,record.seq,quality,barcode)
        r_id = record.id
        #self.name = name
        r_seq = record.seq
        r_len = len(r_seq)
        r_qual = record.letter_annotations['phred_quality']
        r_mean = sum(r_qual)/r_len
        r_len = len(record.seq)
        r_dis = record.description
        for m in ext_d.finditer(r_dis):
            r_bar = m.group(0)[8:]

        #print('c3')
        if barcode != 'a' and r_bar != barcode:
            continue
    # Filter by minimum Phred score
        #print('c4')
        if minPhred > 0 and r_mean < minPhred:
            continue
        # Filter by maximum sequence length
        if maxlen != 0 and r_len > maxlen:
            continue
        # Filter by minimum sequence length
        #print('c5')
        if minlen > 0 and r_len < minlen:
            continue
        # Filter by sequence search
        if seq_search and seq_search not in r_seq:
            continue
        # Check if the total base pair count is below the specified maximum
        if totalb != 0:
            if ttr_len + r_len > totalb:
                break
        # Add the item to the output dictionary
        out.append(r_id)
        ttr_len += r_len
        ttr_n+= 1
        
        if ttr_n <6:
            ttr_ex[r_id] = fastq_BL(r_id,r_id,r_seq,r_qual,r_bar)
        if ttr_bar.get(r_bar) is not None:
            ttr_bar[r_bar]+=1
        else:
            ttr_bar[r_bar]=1
        ttr.append(record)
        #print('c3')
        if ttr_n%10000==0:
            print('[STATUS]: ',ttr_n,'records have been filtered')
        if ttr_n> 100000*kkk:
            if save_file[0]!=0:
                if save_file[0]==1:
                    print('[STATUS]: Saving... ',save_file[1][:-6]+str(kkk)+'.fastq')
                    SeqIO.write(ttr, save_file[1][:-6]+str(kkk)+'.fastq', 'fastq')
                if save_file[0]==2:
                    print('[STATUS]: Saving... ',save_file[1][:-6]+str(kkk)+'.fastq')
                    SeqIO.write(ttr, save_file[1][:-6]+str(kkk)+'.fastq', 'fastq')
                    print('[STATUS]: Saving... ',save_file[1][:-6]+str(kkk)+'.fastq.gz')
                    fq2gz(save_file[1][:-6]+str(kkk)+'.fastq')
            kkk+=1
            ttr=list()
            
    """   
    if save_file[1] == '':
        pass
    else:
        SeqIO.write(record, save_file[1], 'fastq')
"""
    if kkk>1:
        if save_file[0]!=0:
            if save_file[0]==1:
                print('[STATUS]: Saving... ',save_file[1][:-6]+str(kkk)+'.fastq')
                SeqIO.write(ttr, save_file[1][:-6]+str(kkk)+'.fastq', 'fastq')
            if save_file[0]==2:
                print('[STATUS]: Saving... ',save_file[1][:-6]+str(kkk)+'.fastq')
                SeqIO.write(ttr, save_file[1][:-6]+str(kkk)+'.fastq', 'fastq')
                print('[STATUS]: Saving... ',save_file[1][:-6]+str(kkk)+'.fastq.gz')
                fq2gz(save_file[1][:-6]+str(kkk)+'.fastq') 
    else:
        if save_file[0]!=0:
            if save_file[0]==1:
                print('[STATUS]: Saving... ',save_file[1][:-6]+'.fastq')
                SeqIO.write(ttr, save_file[1][:-6]+'.fastq', 'fastq')
            if save_file[0]==2:
                print('[STATUS]: Saving... ',save_file[1][:-6]+'.fastq')
                SeqIO.write(ttr, save_file[1][:-6]+'.fastq', 'fastq')
                print('[STATUS]: Saving... ',save_file[1][:-6]+'.fastq.gz')
                fq2gz(save_file[1][:-6]+'.fastq') 
        """
    print(ttr, save_file)
    SeqIO.write(ttr, save_file, 'fastq')
"""
            
        
    
    #print mini report
    #gi = dict(dict_BL)
    print('\n_____Start_of_mini_report_____\n')
    status_rec(ttr_n)
    l=1
    gi = ttr_ex
    for i in gi:
        print('Record ',l,' :')
        l+=1
        print('ID: ',gi[i].id)
        #print('Name: ',gi[i].name)
        se = gi[i].seq
        le = gi[i].len
        if le >=15:
            print('Sequence: ',se[:6],'...',se[-6:])
        else:
            print('Sequence: ',se)
        print('Seq Length: ',le)
        #print('Quality: ',gi[i].qual)
        print('Mean Quality: ',gi[i].mean)
        print('Barcode: ',gi[i].bar)
        print('\n=============')
        #tlist.append(gi[i].raw)
    print('\n______End_of_mini_report______')
    print('Barcode',ttr_bar)

def saving_dialog():
    out = ""
    sav = input("Would you like to save filtered data? [Y/n] : ")
    if(sav != 'n'):
        while True:
            ty = input("Please put in filename with extention : ")
            ex = f2ext(ty)
            if ex=='folder':
                print('[WARNING]: Please add extension (.fastq,.gz)')
            else:
                break
            
        if ex=='fastq': #.fastq
            out = ty[:-6]+'.fastq'
            casu = [1, out]
        if ex=='gz':
            
            if f2ext(ty[:-3]) != 'fastq':
                out = ty[:-3]+'.fastq'
            else: 
                out = ty[:-3]
            casu = [2, out]
        print("[STATUS]: File will be saved as "+out)
        return casu
    else:
        print("[STATUS]: File will not saved")
        return [0,out]
    
def Omsin(SeqFileName):
    # if len(sys.argv) != 2:
    #     help(BLhelp)
    #     exit()

    fastqFile = SeqFileName

    data_ready = file2fq(fastqFile)
    print('[STATUS]: File is prepared in .fastq extension form')
    gi = fq2dict(data_ready)
    print('[STATUS]: Dictionary of the data is created')

    #ค่า Heading Row ของ Table
    Heading = ["Barcode", "Number of seq", "Average length", "Average QScore", "Density (%)"]
    Row = [["001", "002"], [1000, 500], [500, 250], [100, 200], [50, 40], [0.5,0.4]]
    lendict ,qsc = createdict(gi)
    Row = inputtable(lendict, qsc)

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

def f2ext(f):
    ex=f.split('.')
    a = len(ex)
    if a==1:
        return('folder')
    else:
        return(ex[a-1])
    
def fq2dict(fastqFile):
    print('[STATUS]: '+str(fastqFile)+' is being read')
    fastqFile = fastqFile
    ext_d = re.compile(r'barcode=\S+',re.M)

    nx = input("Type 'a' to read all records /Type [number] to first [number] records : ")

    Isagi = dict()
    n = 0

    if nx == 'a':
        for record in SeqIO.parse(fastqFile,'fastq'):
            quality = record.letter_annotations['phred_quality']
            seq_len = len(record.seq)
            r_dis = record.description
            for m in ext_d.finditer(r_dis):
                barcode = m.group(0)[8:]
            Isagi[record.id] = fastq_BL(record.id,record.id,record.seq,quality,barcode,record)
            n+=1
            #print(n)
        #status_rec(len(Isagi))
            #hoho = input()
    else:
        nx = int(nx)
        for record in SeqIO.parse(fastqFile,'fastq'):
            if n>=nx:
                break
            else:
                quality = record.letter_annotations['phred_quality']
                seq_len = len(record.seq)
                r_dis = record.description
                ID = record.id
                for m in ext_d.finditer(r_dis):
                    barcode = m.group(0)[8:]
                
                Isagi[ID] = fastq_BL(record.id,record.id,record.seq,quality,barcode,record)
                n+=1
    status_rec(len(Isagi))
    return Isagi
    
def filter_BL(dict_BL, minlen=0, maxlen=0, minPhred=0, barcode='a', totalb='a', seq_search=''):
    
    print('[FILTER]: '+'Inputed Parameters:(a=all/max)')
    print('[FILTER]: '+'min seq length\t:',minlen)
    print('[FILTER]: '+'max seq length\t:',maxlen)
    print('[FILTER]: '+'min phred score\t:',minPhred)
    print('[FILTER]: '+'barcode\t\t:',barcode)
    print('[FILTER]: '+'sequence search\t:',seq_search)
    print('[FILTER]: '+'total base pair\t:',totalb)
    print('[FILTER]: '+'Begin filtering process')
    pdic = dict_BL
    pdic2 = dict(pdic)
    out = dict()
    bp_count = 0
    if barcode =='a':
        pass
    else:
        for i in pdic:
            if pdic[i].bar == barcode:
                pass
            else:
                x = pdic2.pop(i)
        pdic= dict(pdic2)
    #print('Checkpoint1: ')
    #print(pdic)
    
    if minPhred == 0:
        pass
    else:
        for i in pdic:
            if pdic[i].mean >= minPhred:
                pass
            else:
                x = pdic2.pop(i)
        pdic= dict(pdic2)
    #print('Checkpoint2: ')
    #print(pdic)
                
    if maxlen==0:
        pass
    else:
        for i in pdic:
            if pdic[i].len <= maxlen:
                pass
            else:
                x = pdic2.pop(i)
        pdic= dict(pdic2)
    #print('Checkpoint3: ')
    #print(pdic)
        
    if minlen==0:
        pass
    else:
        for i in pdic:
            #print(pdic[i].len)
            if pdic[i].len >= minlen:
                pass
            else:
                #print('aw')
                x = pdic2.pop(i)
        pdic= dict(pdic2)
    #print('Checkpoint4: ')
    #print(pdic)
        
    if seq_search=='':
        pass
    else:
        for i in pdic:
            if pdic[i].seq.find(seq_search)==-1:
                x = pdic2.pop(i)
            else:
                pass
        pdic= dict(pdic2)
    #print('Checkpoint5: ')
    #print(pdic)

    if totalb =='a':
        #print('Checkpoint6: ')
        #print(pdic)
        return pdic
    else:
        pdic2=dict()
        for i in pdic:
            k= pdic[i].len + bp_count
            if k > totalb:
                return pdic2
            else:
                bp_count = k
                #print('ok')
                pdic2[i] = pdic[i]
        return pdic2
    print("If this message is shown something went wrong")

def status_rec(n):
    if n == 0:
        print('[STATUS]: There is no record')
        print('=============\n')
    elif n == 1:
        print('[STATUS]: 1 record found')
        print('=============\n')
    else:
        print('[STATUS]: ',n,' records found')
        print('=============\n')
        
def dict_BL_read(dict_BL):
    
    #tlist=list()
    gi = dict(dict_BL)
    print('\n_____Start_of_mini_report_____\n')
    status_rec(len(gi))
    l=1
    
    for i in gi:
        print('Record ',l,' :')
        l+=1
        print('ID: ',gi[i].id)
        #print('Name: ',gi[i].name)
        se = gi[i].seq
        le = gi[i].len
        if le >=15:
            print('Sequence: ',se[:6],'...',se[-6:])
        else:
            print('Sequence: ',se)
        print('Seq Length: ',le)
        #print('Quality: ',gi[i].qual)
        print('Mean Quality: ',gi[i].mean)
        print('Barcode: ',gi[i].bar)
        print('\n=============')
        #tlist.append(gi[i].raw)
    print('\n______End_of_mini_report______')
    
def dict_BL2fq(dict_BL, fq_out_name):
    m=fq_out_name.split('.')
    if m[len(m)-1] != 'fastq':
        fq_out_name += '.fastq'
    l = list()
    for i in dict_BL:
        #print(dict_BL[i])
        l.append(dict_BL[i].raw)
    SeqIO.write(l, fq_out_name, "fastq")
    
def fq2gz(fqfile,filename=''):
    if filename=='':
        filename=fqfile+'.gz'
    #a = fq2dict(fqfile)
    #l = list()
    #for i in a:
    #    l.append(a[i].raw)
    with open(fqfile, 'rb') as f_in:
            with gzip.open(filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

def col_ext(dict_BL):
    id_lis=list()
    len_lis=list()
    mean_lis=list()
    bar_lis=list()
    for i in dict_BL:
        h = dict_BL[i]
        id_lis.append(h.id)
        len_lis.append(h.len)
        mean_lis.append(h.mean)
        bar_lis.append(h.bar)
    return [id_lis,len_lis,mean_lis,bar_lis]

def createdict(dict_BL):
    lendict=dict()
    qsc=dict()
    for i in dict_BL:
        h = dict_BL[i]
        if h.bar not in lendict.keys(): 
            lendict[h.bar] = []
            qsc[h.bar]=[]
        lendict[h.bar].append(h.len)
        qsc[h.bar].append(h.mean)
    return lendict, qsc

#สร้าง input ให้ table (bar,numseq,avglen,avgqsc,den)
def inputtable(lendict, qsc):
    a = dict()
    bar=[]
    count=0
    for i in lendict:
        if i not in a:
            a[i]=[]
            bar.append(i)
        a[i].append(i) #bar
        a[i].append(len(lendict[i])) #numseq
        count = count + len(lendict[i])
        avglen=sum(lendict[i])/len(lendict[i])
        a[i].append(avglen) #avglen
    for i in qsc:
        avgqsc = sum(qsc[i])/len(qsc[i])
        a[i].append(avgqsc) #avgqsc
        den = len(qsc[i])/count *100
        a[i].append(den) #den
    bar.sort()
    numseq=[]
    lenlist=[]
    qsclist=[]
    denlist=[]
    for n in range(len(bar)):
        numseq.append(a[bar[n]][1])
        lenlist.append(a[bar[n]][2])
        qsclist.append(a[bar[n]][3])
        denlist.append(a[bar[n]][4])
    return [bar,numseq,lenlist,qsclist,denlist]
    
def GenerateBinnedList(Xlist = [], nbins = 10):
    if(len(Xlist)>1):
        if(len(Xlist)>= nbins):
            minx = min(Xlist)
            maxx = max(Xlist)
            width = (maxx - minx)/nbins
            out = dict()
            out_x = []
            out_y = []
            for i in Xlist:
                part = (i-minx)//width
                if(i == maxx): part -=1
                if part not in out.keys(): out[part] = 1
                else: out[part]+=1 
            #print(out)
            for part in sorted(out):
                out_x.append(minx + (part+0.5)*width)
                out_y.append(out[part])
            return out_x, out_y
        else: GenerateBinnedList(Xlist, len(Xlist))
    else: return Xlist, [1]

def GeneratePlotObject(dict_BL):
    lendict, qsc = createdict(dict_BL)
    lenplotlist = []
    qscplotlist = []
    #print(lendict)
    for Barcodename, valuelist in lendict.items():
        #print(Barcodename, len(valuelist))
        nbins = math.floor(math.sqrt(len(valuelist)))+1
        x, y = GenerateBinnedList(valuelist, nbins)
        lenplotlist.append(BarcodePlotObject(Barcodename, x, y))
    for Barcodename, valuelist in qsc.items():
        nbins = math.floor(math.sqrt(len(valuelist)))+1
        x, y = GenerateBinnedList(valuelist, nbins)
        qscplotlist.append(BarcodePlotObject(Barcodename, x, y))
    return lenplotlist, qscplotlist     

def file2fq(fastqFile):
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
    return newfq

def GenReadableHTML(name,script):
    #สร้าง full HTML (replace)
    out=""
    with open("Report_template.html", "r") as template:
        temp=template.read()
        out = temp.replace(name, script)
    return out

def GenTableHTML(Name,Heading,Row):
    #สร้าง HTML ของ Table เท่านั้น
    tabletem=""
    with open("Table_Template.html", "r") as table:
        tabletem=table.read()
        tabletem=tabletem.replace("TableID",Name)
        location = re.finditer("\"values\": \[\]", tabletem)
        m1 = next(location)
        m2 = next(location)
        newscript = tabletem[:m1.end()-2]+ str(Row) + tabletem[m1.end():m2.end()-2] + str(Heading) + tabletem[m2.end():]
    return newscript

#สร้าง HTML ของ Histo เท่านั้น
def GenerateHistogramFromObject(Name,Header = "Header", xlabel = "xlabel", ylabel = "Count", BarcodeList = []):
    oldscript = ""
    with open("Histo_template.html", "r") as temp:
        oldscript = temp.read()
        oldscript = oldscript.replace("HistogramID",Name)
        pointer = re.search("data = \[\];", oldscript)
        added_text = ""
        for i in range(len(BarcodeList)):
            added_text += "{\"fill\": \"tozeroy\", \"mode\": \"none\", \"name\": " + "\"" + BarcodeList[i].name + "\", \"showlegend\": true, \"type\": \"scatter\", \"x\": " + str(BarcodeList[i].x)+", \"y\": " + str(BarcodeList[i].y) +"},"
        added_text = added_text.rstrip(",")
        added_text += "]"
        script_updated = oldscript[:pointer.end()-2]+ added_text + oldscript[pointer.end():]
        location = re.finditer("\"text\": \"\"", script_updated)
        m1 = next(location)
        m2 = next(location)
        m3 = next(location)
        newscript = script_updated[:m1.end()-2] + "\"" + Header + "\"" + script_updated[m1.end():m2.end()-2] + "\"" + xlabel + "\"" + script_updated[m2.end():m3.end()-2] + "\"" + ylabel + "\"" + script_updated[m3.end():]
    return newscript
