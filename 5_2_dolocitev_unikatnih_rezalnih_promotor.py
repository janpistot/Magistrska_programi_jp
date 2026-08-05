import os
import shutil
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from primer3 import bindings
from Bio.Restriction import *

#wdir
wdir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir"

#obdelava_introni
obd_int = os.path.join(wdir,"obdelava_intronskih")

#Biodeli
Biodeli = os.path.join(wdir,"Biodeli")

#procesirane sekvence
proc_seq = os.path.join(wdir,"procesirane_sekvence")
vektorji_fasta = os.path.join(proc_seq,"vektorji_fasta")

#restrikcija
res_rez = os.path.join(wdir,"restrikcija")

programi_dir = r"C:\Users\HP\Desktop\working"

#------------------------------------------------------------------------------
for file in os.listdir(res_rez):
    if file in ("rezanje_vektorjev.fasta", "rezanje_vektorjev.txt"):
        os.remove(os.path.join(res_rez,file))
        
#------------------------------------------------------------------------------
from Bio import Restriction
from Bio.Restriction import (CommOnly)
ll = {}
nek = []

def rez2 (inpt):
    
    with open (os.path.join(inpt),"r") as rd:
        enzymes = list(CommOnly)            
                
        for record in SeqIO.parse(rd, "fasta"):
            nerez = []
            ll[record.id] = []
            seq = record.seq

            for enzyme in enzymes:
                envim = str(enzyme)
                enz = getattr(Restriction, envim) 
                cutting_sites = enz.search(seq)
                
                if cutting_sites == []:
                    if envim not in nerez:
                        nerez.append(envim)
                        ll[record.id].extend(nerez)                    
                    
for file in os.listdir(vektorji_fasta):
    vvv = os.path.join(vektorji_fasta, file)
    rez2 (vvv)

for file in os.listdir(Biodeli):
    if file in ("A2", "Kozak", "Promotor", "Terminator"):
        for fil in os.listdir(os.path.join(Biodeli, file)):
            vvv = os.path.join(Biodeli, file, fil)
            rez2 (vvv)

for file in os.listdir(obd_int):
    if file in ("z_obdel_seq", "zdru_neintron"):
        for fil in os.listdir(os.path.join(obd_int, file)):
            vvv = os.path.join(obd_int, file, fil)
            rez2 (vvv)

common_items = set.intersection(*(set(value) for value in ll.values()))

outpt = os.path.join(res_rez, "rezanje_vektorjev.fasta" )
enzyme_batch = RestrictionBatch(common_items)

ambiguous_bases = set("WRYSMKHBVDN")

with open(os.path.join(outpt),"a") as wr:
    for enzyme in enzyme_batch:
        sss = set(enzyme.site)
        c = 0
        for item in ambiguous_bases:
            if item in sss:
                c = 1
        if enzyme and not enzyme.is_blunt() and c == 0:
            print(f"{enzyme}: {enzyme.site}")

            wr.write(">" + str(enzyme) + "\n") 
            wr.write(str(enzyme.site) + "\n")

with open(os.path.join(res_rez,"rezanje_vektorjev.txt"),"a") as rd2:
    for enzyme in common_items:
        rd2.write(enzyme + "\n")
        
subprocess.run(["python", os.path.join(programi_dir, "6_1_Nova.py")])
