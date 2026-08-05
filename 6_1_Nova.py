import os
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from primer3 import bindings
from Bio.Restriction import *

#wdir
wdir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir"

#Primerji
Primerji = os.path.join(wdir,"Primerji")

pr_inpt =  os.path.join(Primerji,"Primer3input")
pr_out = os.path.join(Primerji, "primer3output")

#Biodeli
Biodeli = os.path.join(wdir,"Biodeli")

zapi = os.path.join(Biodeli,"zaporedje","zaporedje.txt")

promotor = os.path.join(Biodeli,"promotor")
dvaA = os.path.join(Biodeli,"A2")
kozak= os.path.join(Biodeli, "Kozak")

#procesirane sekvence
proc_seq = os.path.join(wdir,"procesirane_sekvence")

#restrikcija
res_rez = os.path.join(wdir,"restrikcija")

#Gibson
gibson = os.path.join(wdir,"Gibson")

asembly = os.path.join(gibson,"assembly")

#Encmi_skupaj
encimi_skupaj = os.path.join(wdir,"Encimi_skupaj")

programi_dir = r"C:\Users\HP\Desktop\working"

#------------------------------------------------------------------------------
lista_assembly = {}
pro_list =[]

u_pro = "FUS1promotor-BBa_K2111002"



def purge_folder(folder_path):
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isfile(item_path) or os.path.islink(item_path):
            os.remove(item_path)
for file in os.listdir(encimi_skupaj):
    di = os.path.join(encimi_skupaj, file)
    purge_folder(di)

for file in os.listdir(encimi_skupaj):
    os.rmdir(os.path.join(encimi_skupaj,file))

for file in os.listdir(asembly):
    os.remove(os.path.join(asembly,file))
    
necleni = ["promotor","kozak","A2"]
zaporedje = []
zapredje = []
with open (os.path.join(zapi), "r") as rd:
    for line in rd:
        if line.startswith("A2") or line.startswith("kozak"):
            uu = line.strip()
            bb = uu.split("__")[0]
            zaporedje.append(bb)
        else:
            zaporedje.append(line.strip())
        if not line.startswith(("rez", "A2","kozak")) and not line.strip() in necleni:
            zapredje.append(line.strip())

with open(os.path.join(pr_inpt,"primer3inpt.fasta"),"r") as inp:
    for record in SeqIO.parse(inp,"fasta"):
        vrsta = record.id.split("__")[-1]
        
        if vrsta in zapredje:
            ime = record.id.split("|")[1]
            
            if not os.path.isdir(os.path.join(encimi_skupaj,vrsta)):
                os.mkdir(os.path.join(encimi_skupaj,vrsta))
                
            with open (os.path.join(encimi_skupaj,vrsta,vrsta +".fasta"),"a") as wr:
                wr.write(">" + record.id + "\n")
                wr.write (str(record.seq) + "\n")

    
for file in os.listdir(pr_out):
    aa = os.path.join(pr_out, file)
    if not os.path.getsize(aa) == 0 and file.endswith("promotor.fasta"):
       bb = file.split("__")[1]
       cd = bb.replace(".fasta","")
       pro_list.append(cd)
       
lista_assembly["promotor"]=[]
pro = []
for file in os.listdir(promotor):
   cc = file.split(".fasta")[0]
   if cc in pro_list:
    pro.append(cc)

lista_assembly["promotor"].extend(pro)


lista_assembly["kozak"]=[]
with open(os.path.join(kozak,"vsi_kozak.fasta"),"r") as koz:
    if not record.id == "":
        for record in SeqIO.parse(koz,"fasta"):
           lista_assembly["kozak"].append(record.id)


plo = 0
lista_assembly["rez__pro"]=[]
por = []
with open(os.path.join(res_rez,"rezanje_vektorjev.fasta"),"r") as rd9:
    for record in SeqIO.parse(rd9,"fasta"):
        if plo == 0:
            por.append(record.id)
            plo = plo + 1
        else:
             break
         
lista_assembly["rez__pro"].extend(por)

lista_assembly["A2"] = []
with open(os.path.join(dvaA, "A2.fasta"),"r") as A2:
    for record in SeqIO.parse(A2,"fasta"):
        if not record.id == "":
            lista_assembly["A2"].append(record.id)


def ekstrakt_list(ipt):
    inp = []
    lista_assembly[ipt]=[]
    with open (os.path.join(encimi_skupaj,ipt,ipt+ ".fasta"),"r") as rd:
        for record in SeqIO.parse(rd,"fasta"):
            inp.append(record.id.split("|")[1])           
    lista_assembly[ipt].extend(inp)
      
for file in os.listdir(encimi_skupaj):
    ekstrakt_list(file)

n = 0
z = 0
A2_count = 0
kozak_count = 0

bar = []

zap = 0

print(lista_assembly)
for i in zaporedje:
    zap = zap + 1
    if i not in ("rez__zac", "rez__kon"):
        
        n = n + 1
        globals()["__" + str(n)] = {}
        bar.append(n)
        
        z = 0
        if i.startswith("A2"):
            u = "A2"
            A2_count =  A2_count + 1
            
        elif i.startswith("kozak"):
            u = "kozak"
            kozak_count = kozak_count + 1
            
        else:
            u = i
        
        for value in lista_assembly[u]:
            if not n == 1:
                if u in("kozak", "A2"):
                    
                    for valu in globals()["__" + str(n - 1)].values():
                        
                        if value == lista_assembly[u][globals()[u + "_count"]-1]:
                            
                            print(value)
                            
                            z = z + 1
                            a = []
                            a.extend(valu)
                            globals()["__" + str(n)]["__" + str(z)] = a
                            globals()["__" + str(n)]["__" + str(z)].append(value)
                
                else:
                    for valu in globals()["__" + str(n - 1)].values():
                        
                        z = z + 1
                        
                        a = []
                        a.extend(valu)
                        globals()["__" + str(n)]["__" + str(z)] = a
                        globals()["__" + str(n)]["__" + str(z)].append(value)
                   
     
            else:
                z = 0
                for value in lista_assembly[i]:
                    if u in("kozak", "A2"):
                        if value == lista_assembly[u][globals()[u + "_count"]-1]:
                            globals()["__" + str(n)]["__" + str(z)] = []
                            globals()["__" + str(n)]["__" + str(z)].append(value)
                            
                    else:       
                        z = z + 1
                        globals()["__" + str(n)]["__" + str(z)] = []
                        globals()["__" + str(n)]["__" + str(z)].append(value)
                    
combinations = (globals()["__" + str(bar[-1])])

with open(os.path.join(asembly,"asembly.fasta"),"a") as wr:
    for st in combinations.values():
        
        
        cc = "__".join(st)
        print(cc)
        
        wr.write ("\n" + ">" + cc + "\n")
        for i in st:
            en = ""
            dA = ""
            if i in lista_assembly["promotor"]:
                with open(os.path.join(promotor, i + ".fasta"),"r") as rd1:
                    for rcr in SeqIO.parse(rd1, "fasta"):
                        if rcr.id == u_pro:
                            wr.write(str(rcr.seq))
           
            elif i in lista_assembly["rez__pro"]:
                with open(os.path.join(res_rez,"rezanje_vektorjev.fasta"),"r") as rd5:
                          for r in SeqIO.parse(rd5,"fasta"):
                              if r.id in lista_assembly["rez__pro"]:
                                  print(r.id)
                                  wr.write(str(r.seq))
           
            elif i in lista_assembly["kozak"]:
                with open(os.path.join(kozak,"vsi_kozak.fasta"),"r") as rd4:
                    for rca in SeqIO.parse(rd4, "fasta"):
                        if rca.id == i:
                            wr.write(str(rca.seq))
            
            elif i in lista_assembly["A2"]:
                with open(os.path.join(dvaA,"A2.fasta"),"r") as rd3:
                    for rec in SeqIO.parse(rd3, "fasta"):
                        if rec.id == i:
                            wr.write(str(rec.seq))
            else:
                nek = 0
                for key, value in lista_assembly.items():
                    if nek == 0:
                        for z in value:
                            if i == z:
                                en = key
                                nek = 1
                                break
                path = os.path.join(encimi_skupaj, en, en + ".fasta")
                
                with open(path,"r") as rd2:
                    for recrd in SeqIO.parse(rd2,"fasta"):
                        if recrd.id.split("|")[1] == i:
                            wr.write(str(recrd.seq))

        
subprocess.run(["python", os.path.join(programi_dir, "6_2_Gibson_razrez.py")])
        
