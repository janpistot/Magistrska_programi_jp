import os
import shutil
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from primer3 import bindings
from Bio.Restriction import *

#wdir
wdir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir"

#Biodeli
Biodeli = os.path.join(wdir,"Biodeli")
res_encm = os.path.join(Biodeli,"res_encm")

#procesirane sekvence
proc_seq = os.path.join(wdir,"procesirane_sekvence")

#restrikcija
res_rez = os.path.join(wdir,"restrikcija")

#Gibson
gibson = os.path.join(wdir,"Gibson")

asembly = os.path.join(gibson,"assembly")
rezrez = os.path.join(gibson, "razrez")


programi_dir = r"C:\Users\HP\Desktop\working"

#------------------------------------------------------------------------------
from Bio import Restriction
os.remove(os.path.join(rezrez,"gibson_razrez.txt"))
os.remove(os.path.join(rezrez,"promotor_razrez.txt"))

def rez (inpt, encimi, outpt):
    with open (os.path.join(inpt),"r") as rd, open (os.path.join(encimi), "r") as encm, open(os.path.join(outpt),"a") as wr:
        enzymes = []
        for line in encm:
            enzymes.append(line.strip())
        for record in SeqIO.parse(rd, "fasta"):
            seq = record.seq
            wr.write("\n"+">"+ record.description + "\n")
            for enzyme in enzymes:
                enz = getattr(Restriction, enzyme) 
                cutting_sites = enz.search(seq)
                if cutting_sites == []:
                    wr.write(enzyme + str(cutting_sites)+ "__")   

                       
aaa = os.path.join(asembly, "asembly.fasta")
bbb = os.path.join(res_encm, "rezalni ecimi v pRG MCS.txt")
ccc = os.path.join(rezrez, "gibson_razrez.txt")
rez(aaa, bbb, ccc)


def rezprom (inpt, encimi, outpt):
    with open (os.path.join(inpt),"r") as rd, open (os.path.join(encimi), "r") as encm, open(os.path.join(outpt),"a") as wr:
        enzymes = []
        for line in encm:
            enzymes.append(line.strip())
        for record in SeqIO.parse(rd, "fasta"):
            seq = record.seq
            bb = record.id.split("__")[1]
            wr.write("\n"+">"+ record.description + "\n")
            for enzyme in enzymes:
                if enzyme == bb:
                    enz = getattr(Restriction, enzyme) 
                    cutting_sites = enz.search(seq)
                    if len(cutting_sites) == 1:
                        wr.write(enzyme + str(cutting_sites)+ "__")
                    else:
                        print(record.id)
                        
                    
bbb = os.path.join(res_rez, "rezanje_vektorjev.txt")
ccc = os.path.join(rezrez, "promotor_razrez.txt")
rezprom(aaa, bbb, ccc)


subprocess.run(["python", os.path.join(programi_dir, "7_nova.py")])