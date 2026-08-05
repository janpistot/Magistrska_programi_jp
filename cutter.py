
import os
import shutil
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from primer3 import bindings

#wdir
wdir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir"

#obdelava_introni
obd_int = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\obdelava_intronskih"

introni = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\obdelava_intronskih\introni"
neintron = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\obdelava_intronskih\neintron"
wfol = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\obdelava_intronskih\working_folder"
zdru_neint = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\obdelava_intronskih\zdru_neintron"
zdr_splc = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\obdelava_intronskih\zdruzene_splc"
obd_seq = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\obdelava_intronskih\z_obdel_seq"

#BlastDB
blastdb_dir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Blastdb"

primerBDB = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Blastdb\primer_blastDB"

#Primerji
Primerji = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Primerji"

primer_input = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Primerji\primer3input"
primerBLAST = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Primerji\Primer_blast"
pr_inpt = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Primerji\Primer3input"
pr_out = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Primerji\primer3output"

#Biodeli
Biodeli = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Biodeli"

gen_dir =r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Biodeli\genbank_genomi"
encimi_dir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Biodeli\Encimi"
vektorjigb = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Biodeli\Vektorjigb"
res_encm = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Biodeli\res_encm"

#procesirane sekvence
proc_seq = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\procesirane_sekvence"

out = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\procesirane_sekvence\genbanksplit"
in_dir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\procesirane_sekvence\Genomi"
file_path = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\procesirane_sekvence\concat"
vektorji_fasta = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\procesirane_sekvence\vektorji_fasta"

#Blastoutput-i
B_out = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Blastoutput"

Blastoutput =  r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Blastoutput\Blastout"
E_out = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Blastoutput\E_out"
u_blast = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\Blastoutput\uporabniBLAST"

#uporabni_encimi
u_encm = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\uporabni_encimi"

dvojne = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\uporabni_encimi\sekvence"
encimi_list = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\uporabni_encimi\encimi_list"

#restrikcija
res_rez = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir\restrikcija"



evalue = str(10e-5)
#------------------------------------------------------------------------------

from Bio import Restriction
from collections import Counter


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
                       
aaa = os.path.join(zdru_neint, "zdru_neint.fasta")
bbb = os.path.join(res_encm, "rezalni ecimi v pRG MCS.txt")
ccc = os.path.join(res_rez, "zareze.txt")
rez(aaa, bbb, ccc)

aaa = os.path.join (obd_seq, "zdr_u_seq.fasta")
rez(aaa, bbb, ccc)


for file in os.listdir(Biodeli):
    if file in("2A", "Kozak", "promotor", "terminator"):
        for fil in os.listdir (os.path.join(Biodeli,file)):
            aaa = os.path.join(Biodeli,file,fil)
            rez(aaa, bbb, ccc)

count = Counter(last)
print(count) 
