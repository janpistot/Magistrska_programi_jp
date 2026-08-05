
import os
import shutil
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from primer3 import bindings

#wdir
wdir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir"

#obdelava_introni
obd_int = os.path.join(wdir,"obdelava_intronskih")

obd_seq = os.path.join(obd_int, "z_obdel_seq")

#Biodeli
Biodeli = os.path.join(wdir,"Biodeli")

res_encm = os.path.join(Biodeli,"res_encm")

#restrikcija
res_rez = os.path.join(wdir,"restrikcija")


programi_dir = r"C:\Users\HP\Desktop\working"

#------------------------------------------------------------------------------

from Bio import Restriction
from collections import Counter

def purge_folder(folder_path):
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isfile(item_path) or os.path.islink(item_path):
            os.remove(item_path)

purge_folder(res_rez)

last = []
def rez (inpt, encimi, outpt):
    with open (os.path.join(inpt),"r") as rd, open (os.path.join(encimi), "r") as encm, open(os.path.join(outpt),"a") as wr:
        enzymes = []
        for line in encm:
            enzymes.append(line.strip())
        lin = ""   
        for record in SeqIO.parse(rd, "fasta"):
            seq = record.seq
            for enzyme in enzymes:
                enz = getattr(Restriction, enzyme) 
                cutting_sites = enz.search(seq)
                if not cutting_sites == []:
                    last.append(enzyme)
                    try:    
                        a = record.id
                        b = a.split("|")[1]
                    except:
                        b = record.id
                        
                    if not lin == b:
                        wr.write("\n" + b + "-----" + enzyme + "_:_" + str(cutting_sites))
                        lin = b
                    else:
                        wr.write("-----" + enzyme + "_:_" + str(cutting_sites))
                       
#aaa = os.path.join(zdru_neint, "zdru_neint.fasta")
bbb = os.path.join(res_encm, "rezalni ecimi v pRG MCS.txt")
ccc = os.path.join(res_rez, "zareze.txt")
#rez(aaa, bbb, ccc)

aaa = os.path.join (obd_seq, "zdr_u_seq.fasta")
rez(aaa, bbb, ccc)


for file in os.listdir(Biodeli):
    if file in("A2", "Kozak", "promotor", "terminator"):
        for fil in os.listdir (os.path.join(Biodeli,file)):
            aaa = os.path.join(Biodeli,file,fil)
            rez(aaa, bbb, ccc)

count = Counter(last)
print(count) 

subprocess.run(["python", os.path.join(programi_dir, "5_2_dolocitev_unikatnih_rezalnih_promotor.py")])
