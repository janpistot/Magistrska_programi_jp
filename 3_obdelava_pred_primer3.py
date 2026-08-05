import os
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq

#wdir
wdir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir"

#obdelava_introni
obd_int = os.path.join(wdir,"obdelava_intronskih")

obd_seq = os.path.join(obd_int, "z_obdel_seq")

#Primerji
Primerji = os.path.join(wdir,"Primerji")
primer_input = os.path.join(Primerji,"primer3input")

#Biodeli
Biodeli = os.path.join(wdir,"Biodeli")
dvaA = os.path.join(Biodeli,"A2")

zapi = os.path.join(Biodeli,"zaporedje","zaporedje.txt")
promotor = os.path.join(Biodeli,"promotor")

programi_dir = r"C:\Users\HP\Desktop\working"

#------------------------------------------------------------------------------
def purge_folder(folder_path):
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isfile(item_path) or os.path.islink(item_path):
            os.remove(item_path)

for file in os.listdir(Primerji):
    if file in ("primer3input"):
        purge_folder(os.path.join(Primerji, file))
        
#------------------------------------------------------------------------------

STOP_kodoni = {Seq("TAG"), Seq("TAA"), Seq("TGA")}

zaporedje = []
with open (os.path.join(zapi), "r") as rd:
    for line in rd:
        if not line.startswith("rez"):
            if line.startswith("A2"):
                uu = line.strip()
                bb = uu.split("__")[0]
                zaporedje.append(bb)
            else:
                zaporedje.append(line.strip())
            
print(zaporedje)


lag = []
for file in os.listdir(promotor):
    if "promotor" in zaporedje:
        with open(os.path.join(promotor, file),"r") as pr:
            for recr in SeqIO.parse(pr,"fasta"):
                v = recr.id
                lag.append(v)
    else:
        break

def obdelava_primer_encimi (inpt, outpt):
    with open(os.path.join(inpt),"r") as neint, open(os.path.join(outpt),"a") as noSTOP:
        
        for record in SeqIO.parse(neint, "fasta"):
                
            if record.id in lag and not str(record.id) == "":
                noSTOP.write(">" + record.id + "__" + "promotor" +"\n")
                noSTOP.write(str(record.seq) + "\n")
                
            elif not str(record.id) == "":
                
                c = record.id.split("|")[0]
                vrsta_enc = c.split("_")[1]
                if vrsta_enc in zaporedje:
                    
                    zap = zaporedje.index(vrsta_enc)
                    
                    if zaporedje[zap - 1].startswith("kozak"):
                        if vrsta_enc == zaporedje[-1]:
                            noSTOP.write(">" + record.id + "__"+ vrsta_enc +"\n")
                            noSTOP.write(str(record.seq[3:]) + "\n")
                        else:
                            noSTOP.write(">" + record.id + "__"+ vrsta_enc +"\n")
                            sekv = record.seq[3:]
                            while sekv[-3:] in STOP_kodoni:
                                sekv = sekv[:-3]
                            noSTOP.write(str(sekv) + "\n")
                    else:
                        if vrsta_enc == zaporedje[-1]:
                            noSTOP.write(">" + record.id + "__"+ vrsta_enc +"\n")
                            noSTOP.write(str(record.seq) + "\n")
                        else:
                            sekv = record.seq[3:]
                            while sekv[-3:] in STOP_kodoni:
                                sekv = sekv[:-3]
                                print(sekv)
                            noSTOP.write(">" + record.id + "__"+ vrsta_enc +"\n")
                            noSTOP.write(str(sekv) + "\n")

                              
outpt = os.path.join(primer_input, "primer3inpt.fasta")
with open(os.path.join(outpt),"a") as noSTOP:
    noSTOP.write(">DaNeCrkne" + "\n")
    
inpt = os.path.join(obd_seq,"zdr_u_seq.fasta")
obdelava_primer_encimi(inpt, outpt)

for file in os.listdir(promotor):
    if "promotor" in zaporedje:
        inpak = os.path.join(promotor, file)
        obdelava_primer_encimi(inpak, outpt)
    else:
        break


subprocess.run(["python", os.path.join(programi_dir, "4_primer_design1.py")])