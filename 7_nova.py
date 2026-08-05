import os
from Bio import SeqIO
#import subprocess
from Bio.Seq import Seq
from Bio.Restriction import *
from Bio.Restriction import RestrictionBatch, AllEnzymes
import itertools

from primer3 import bindings
from primer3 import calc_tm
from primer3 import calc_hairpin
from primer3 import calc_homodimer
from primer3 import calc_heterodimer

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import csv

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Promotor_v_vektorju = True
Terminator_v_vektorju = True
podaljski = 18 #maksimalna dolžina podaljškov primerjev
koncni_primerji_z_2A_podaljski = ["r_ICS","f_SAMT"]
st_pr = 1

zac_rez = Seq("AA")
kon_rez =Seq("TATA")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#wdir
wdir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir"

#Primerji
Primerji = os.path.join(wdir,"Primerji")
pr_out = os.path.join(Primerji, "primer3output")
primer_input = os.path.join(Primerji,"primer3input")


#Biodeli
Biodeli = os.path.join(wdir,"Biodeli")

dvaA = os.path.join(Biodeli,"A2")
zapi = os.path.join(Biodeli,"zaporedje","zaporedje.txt")
res_encm = os.path.join(Biodeli,"res_encm")
kozak = os.path.join(Biodeli,"Kozak")

#restrikcija
res_rez = os.path.join(wdir,"restrikcija")

#Gibson
gibson = os.path.join(wdir,"Gibson")

asembly = os.path.join(gibson,"assembly")
rezrez = os.path.join(gibson, "razrez")

#rezultati
Rezultati = os.path.join(wdir,"Rezultati")

Genbank_rezultati = os.path.join(Rezultati, "anotacije")
Skrajsan_GB = os.path.join(Rezultati,"Kratke_sekvence")

#obdelava_introni
obd_int = os.path.join(wdir,"obdelava_intronskih")

obd_seq = os.path.join(obd_int, "z_obdel_seq")

#------------------------------------------------------------------------------

if os.path.exists(os.path.join(pr_out,"CDS.fasta")):
    os.remove(os.path.join(pr_out,"CDS.fasta"))
if os.path.exists(os.path.join(pr_out,"reverz.fasta")):
    os.remove(os.path.join(pr_out,"reverz.fasta"))

if os.path.exists(os.path.join(Rezultati,"primerji.fasta")):
    os.remove(os.path.join(Rezultati,"primerji.fasta"))
if os.path.exists(os.path.join(Rezultati,"asembly.fasta")):
    os.remove(os.path.join(Rezultati,"asembly.fasta"))
    
if os.path.exists(os.path.join(Rezultati,"asembly_koncni.fasta")):
    os.remove(os.path.join(Rezultati,"asembly_koncni.fasta"))

for file in os.listdir(os.path.join(Rezultati, "Kratke_sekvence")):
    os.remove(os.path.join(Rezultati,"Kratke_sekvence",file))

if os.path.exists(os.path.join(Rezultati, "A2_primerji.fasta")):
    os.remove(os.path.join(Rezultati, "A2_primerji.fasta"))

for file in os.listdir(Genbank_rezultati):
    os.remove(os.path.join(Genbank_rezultati,file))

if os.path.exists(os.path.join(Rezultati, "primerji.csv")):
    os.remove(os.path.join(Rezultati, "primerji.csv"))

if os.path.exists(os.path.join(Rezultati, "drugi_pomnozek.fasta")):
    os.remove(os.path.join(Rezultati, "drugi_pomnozek.fasta"))
if os.path.exists(os.path.join(Rezultati, "prvi_pomnozek.fasta")):
    os.remove(os.path.join(Rezultati, "prvi_pomnozek.fasta"))
#------------------------------------------------------------------------------
with open(os.path.join(pr_out,"CDS.fasta"),"a") as wr: #brezveze, popravi v prejšnjem programu, vse napiše v eno, da ni treba združevati
    for file in os.listdir(pr_out):
        if file not in ("CDS.fasta"):
            with open(os.path.join(pr_out,file),"r") as rd:
                for record in SeqIO.parse(rd,"fasta"):
                    wr.write("\n" + ">" + record.id)
                    wr.write("\n" + str(record.seq))

#--------------preberi zaporedje, inicializiraj vse knjižnice------------------
kok = []
for k in koncni_primerji_z_2A_podaljski:
    op = k.split("_")
    kok.append(op[1])
    
zaporedje = []
zapiredje = []
with open (os.path.join(zapi), "r") as rd:
    for line in rd:
        zapiredje.append(line.strip())
        if line.startswith("A2") or line.startswith("kozak"):
            cc = line.strip()
            dd = cc.split("__")[0]
            zaporedje.append(dd)
            globals()[dd + "_dic"] = {}
            
        else:
            zaporedje.append(line.strip())
            variable_name = line.strip()
            globals()[variable_name + "_dic"] = {} #inicializacija knjižnice vsakega dela (ICS, IPL, SAMT, kozak, A2, rezalna mesta)


#-------------------------konstrukcija knjižnic--------------------------------
CDS_primer = {}# iz fila da kodo v knižnico -> popraavi, brezveze da najprej naredi file, pol pa iz njega bere, združi s prejšnjo fazo, brezveze, da se file združuje, samu v enega bi se moglo pisat.
koz1 = 0
for i in zaporedje:
    if not i.startswith("rez") and not (i in ("kozak","A2")):
        with open(os.path.join(pr_out,"CDS.fasta"),"r") as cd:
            for record in SeqIO.parse(cd,"fasta"):
                try:
                    prot = record.id.split("|")[1]
                except:
                    prot = record.id.split("__")[2]
                
                CDS_primer[record.id] = record.seq
                    
                if record.id.split("__")[-2] ==  i:
                    if prot not in globals()[i + "_dic"]:
                        globals()[i + "_dic"][prot] = [record.seq]
                    else:
                        globals()[i + "_dic"][prot].append(record.seq)
                       
        
    elif i == "kozak": #če v sekvenci je kozak, 
        if koz1 == 0:# preberi vse kozake v bazi le enkrat, daj v knjižnico kozakov
            koz1 = koz1 + 1 
            with open(os.path.join(kozak,"vsi_kozak.fasta"),"r") as rd8:
                for record in SeqIO.parse(rd8,"fasta"):
                        if record.id not in globals()[i + "_dic"]:
                            globals()[i + "_dic"][record.id] = [record.seq]
                            globals()[i + "_dic"][record.id].append(record.seq.reverse_complement())
                        else:
                            globals()[i + "_dic"][record.id].append(record.seq)


    elif i ==  "rez__pro": #če v sekvenci interno rezalno mesto
        
        plo = 0
        with open(os.path.join(res_rez,"rezanje_vektorjev.fasta"),"r") as rd9:
            for record in SeqIO.parse(rd9,"fasta"):
                if plo == 0: #preberi potencialna rezalna mesta le enkrat, daj v knjižnico 
                    if record.id not in globals()[i + "_dic"]:
                        globals()[i + "_dic"][record.id] = [record.seq]
                    else:
                        globals()[i + "_dic"][record.id].append(record.seq)
                    plo = plo + 1
                else:

                     break      
A2_dic = {} #inicializacija A2 knjižnice, ne vem zakaj sem to dal posebej, praviloma naj bi program bil za združevanje z A2, ampak ni le to
with open(os.path.join(dvaA,"A2.fasta"),"r") as cd:
    for record in SeqIO.parse(cd,"fasta"):
        if not record.id == "":
            if record.id not in A2_dic:
               A2_dic[record.id] = [record.seq]     #je pa koda veliko preprostejša in sproži incializacijo le enkrat
            else:
               A2_dic[record.id].append(record.seq)

primerji = {}
rezalna_mesta = {}
vse_seq = []
with open(os.path.join(rezrez,"gibson_razrez.txt")) as rd: #preberi rezalna mesta za insercijo, tiste ki ne režejo znotraj inserta in režejo v vektorju
    for line in rd: #veliko bolje bi bilo če bi tudi ta del združil s prejšnjim programom, da ni treba brt iz filov (program za določanje restrikcijskih mest)
        if line.startswith(">"):
            cde = line.split(">")[1]# ni fasta file, ampak sem ohranil fasta headerje, tako da ga ne berem kot sekvenco, ima le info o imenih encimov v viabilnih parih
            cfg = cde.strip()
            rezalna_mesta[cfg] = []
            vse_seq.append(cfg)
        else:
            enc = line.split("__")
            for i in enc[0:-1]:
                e = i.split("[]")[0].strip()
                rezalna_mesta[cfg].append(e)

#-Kombiniranje vseh viabilnih kombinacij ustreznih restrikcijskih encimov za insercijo---
r_zap = []
#imam zapisano zaporedje v MCS vektorja, glede na to ali sta terminator in promotor v vektorju že prisotna je orientacija pomembna ali pa ni.
#ne moremo okoli obračati vektorja, če je karkoli od tega že prisotno, v nasprotnem primeru, pa je orienatcija nepomembna
#hkrati določa pozicija terminatorja ali promotorja to kako bo insert obrnjen

with open(os.path.join(res_encm, "rezalni ecimi v pRG MCS.txt")) as rd:
    for line in rd:
        r_zap.append(line.strip()) #preberi vse rezalne encime v vektorju in jih daj na list
      
com = list(itertools.combinations(r_zap, 2)) #naredi vse kombinacije dveh encimov v MCS vektorja
# list setov ("[('XbaI', 'ClaI'), ('XbaI', 'SalI'), ('XbaI', 'XhoI'), ('ClaI', 'SalI'), ('ClaI', 'XhoI'), ('SalI', 'XhoI')]")

u_kombinacije = [] #inicializacija liste uporabnih kombinacij

for i in com:
    set(i)
    if Promotor_v_vektorju is False and Terminator_v_vektorju is False: #orientacija je nepomembna, ni ne promtorja ne terminatorja, največ nmožnih kombinacij
        b = [] 
        b.extend(r_zap) #b je svoja spremenljivka, ker potem brišem člene
        indx = []
        
        for z in i:
            pos = r_zap.index(z) # slaba koda, naredi listo iz indeksov: kateri encim v setu dveh encimov se nahaja kje na vektorju
            indx.append(pos)
        
        a = int(indx[0])
        c = int(indx[1]) + 1
        
        del b [a:c] #brišem člene v b, da vidim kaj vse po razrezu ostane na vektorju (razrez od enega rezalnega mesta pa do drugega za vsako kombinacijo rezalnih mest)

        cam = list(itertools.combinations(b, 2)) #od ostalih rezalnih mest kakšne so kombinacije -> namen naj bi bil, da ostanejo še rezalna mesta za promotor in terminator: od izrezanih mest ali sta na voljo dve kombinaciji za promotor in dve za terminator.
        
        
        for d in cam: #enako kot prej, preverjam kombinacije za promotor in terminator
            set(d)
            mah = []
            for k in d:
                g = r_zap.index(k)
                mah.append(g)
                
            km = (int(mah[1]) - int(mah[0])) #če obstaja najmanj ena kombinacija rezalnih encimov, ki po izrezu ostanejo in so po indexu zaporedja eden zraven drugega-> kombinacija encimov viabilna
            # vedno tudi če jih je veče viabilnih eden zraven drugega, bo check pozitiven
        
        
            if km == 1 and i not in u_kombinacije:
                u_kombinacije.append(i)
    elif Promotor_v_vektorju is False and Terminator_v_vektorju is True: # terminator je, promotorja pa ni, orientacijo določa promotor, rezalni encimi morajo ostati prosti zadaj
        for x in i:
            if not x in r_zap[0:1]:
                u_kombinacije.append(i)
    elif Promotor_v_vektorju is True and Terminator_v_vektorju is False: #orientacijo določa promotor, rezalni encimi morajo ostati prosti spredaj
         for x in i:
             if not x in r_zap[-2:-1]:
                 u_kombinacije.append(i)
    elif Promotor_v_vektorju is True and Terminator_v_vektorju is True: #orientacijo določa oboje, ne rabimo ekstra rezalnih mest za insert promotorja in terminatorja.
        u_kombinacije.append(i)
dolocanje_kom = {}

for kom in u_kombinacije:
    kim = kom[0] + "__"+ kom [1]
    dolocanje_kom[kim] = [] #vse uporabne kombinacije zberem v listo, kjer je vsaka uporabna kombinacija člen z encimi ločenimi z "__" (slaba koda, potem ločujem, bolje bi bilo, če bi to bil set, popravi)
    
for key,value in rezalna_mesta.items():
    combinations = list(itertools.combinations(value, 2)) #med vsemi ustreznimi kombinacijami preverjam vse kombinavije encimov, ki so pomembne nam.
    
    for njeh in combinations:
        ab = njeh[0] + "__" + njeh[1]
        for ky,val in dolocanje_kom.items():
            if ky == ab:
                dolocanje_kom[ky].append(key) # lista ustreznih kombinacij restrikcijskih encimov za določen set proteinov(MenF, EntC,AtICS1...) v razredu (ICS, IPL, SAMT), kot jih želimo izražati
                
kak = False
for  i in range(1,3):
    if kak is False:
        aa = list(dolocanje_kom.keys())
        kombinacije =list(itertools.combinations(aa,i))
    else:
        break
    for h in kombinacije:
        if kak is False:
            u = []
            for f in h:
                for key,value in dolocanje_kom.items():
                    if key == f:
                        u.extend(value)
            ga = set(u)
            meha = set(vse_seq)
            if ga == meha:
                kak = True
                lista_uporabnih = h
                break
        else:
            break

sekvence = []

for k in r_zap: #pridobi sekvence za vse restrikcijske encime v MCS (bolje kot da pridobivam vsakič znova)
    enzyme = getattr(Restriction, k, None)
    sekvence.append(enzyme.site)

dejanska_u = {} 

for i in lista_uporabnih: #pridobi imena in sekvence danih encimov
    if not i == "":
        for key,value in dolocanje_kom.items():
            if key == i:
                zeh = key.split("__")
                gah = []
                for o in zeh:
                    gah.append(o)
                for cb in value:
                    dejanska_u[cb] = []
                    dejanska_u[cb].extend(gah)
zp = {}

for i in lista_uporabnih: #če je usmerjenost pomembna (terminator/promotor že notri) => le ena pravilna orientacija
    smr = True
    splat = i.split("__")
    zp[i] = []
    
    if Promotor_v_vektorju is False and Terminator_v_vektorju is False: # če ni ne enega ne drugega, potem lahko na katerikoli konec damo katerokoli od restrikcijskih zaporedij
        if splat[1] in r_zap[-2:]:
            smr = False
    else:
        smr = True
        
    for p in splat:
        k = r_zap.index(p)
        m = sekvence[k]
        
        if smr is True:
            zp[i].append(m)
        elif smr is False:
            zp[i].insert(0, m)

primer_elementi = []
zapzap = []
for i in zapiredje:
    if not i.startswith(("rez","kozak","A2")):
        zapzap.append(i)

#______________________________________________________________________________
##########################_KATERE ELEMENTE POTREBUJE PRIMER_###################
delcki = {}
for z in zapzap:
    zapo_ind = zapiredje.index(z)
    zap_ind = zapzap.index(z)
    if zap_ind == 0:
        zapo_ind_nas = zapiredje.index(zapzap[zap_ind + 1])        
        delcki[z] =  zapiredje[0:zapo_ind_nas + 1]
    elif z == zapzap[-1]:
        zapo_ind_prej = zapiredje.index(zapzap[zap_ind - 1])
        delcki[z] = zapiredje[zapo_ind_prej:]
    else:
        zapo_ind_nas = zapiredje.index(zapzap[zap_ind + 1])
        zapo_ind_prej = zapiredje.index(zapzap[zap_ind - 1])
        delcki[z] = zapiredje[zapo_ind_prej:zapo_ind_nas  + 1]
        
#______________________________________________________________________________
rez__zac_dic = {}
rez__kon_dic = {}

for key,value in zp.items():
    a = key.split("__")[0]
    b = key.split("__")[1]
    rez__zac_dic[a] = [(Seq(value[0])), (Seq(value[0]))]
    rez__kon_dic[b] = [(Seq(value[1])), (Seq(value[1]))]
#______________________________________________________________________________

try:
    A2_list =[]
    for key in A2_dic.keys():
        A2_list.append(key)
    A2_zap = []
    for i in zapiredje:
        if i.startswith("A2"):
            A2_zap.append(i)
except:
    print("ni A2")

try:
    kozak_list =[]
    for key in kozak_dic.keys():
        kozak_list.append(key)
    kozak_zap = []
    for i in zapiredje:
        if i.startswith("kozak"):
            kozak_zap.append(i)
except:
    print("ni kozak")


forward__primer = {} # setavljanje primerjev:prejšnji primer že imam, zdaj želim sestaviti skupaj primer iz vseh delov, ki jih potrebujem
# primerje sem in bom dobil iz primer 3 programa po sekvencah posameznega proteina
# alternativa sestavljanju bi bilo iskanje po znani združeni sekvenci
for x in zapzap:
    for key, value in delcki.items():
        if key == x:
            y = value.index(x)
            forward_primer = value[:y + 1] #lista prejšnjega, zdajšnjega in sledečega biodela za konstrukcijo primerja, za forward primer rabim le prejšnji del, za reverse primer pa sledečega
            num = 0    
            for z in forward_primer:
                globals()["f_" + str(num)] = {}
                
                if not num == 0:
                    for ka, va in globals()["f_" + str(num-1)].items():
                        if z.startswith("A2"):
                            u = "A2"
                        elif z.startswith("kozak"):
                            u = "kozak"
                        else:
                            u = z
                            
                        for key,value in globals()[u + "_dic"].items():
                            if z.startswith("A2"):
                                ko = A2_list[A2_zap.index(z)]
                                if key == ko:
                                    ki = ka + "__" + key
                                    seqe = va[0] + value[0]
                                    globals()["f_" + str(num)][ki] = [seqe]
                                    
                            elif z.startswith("kozak"):
                                ko = kozak_list[kozak_zap.index(z)]
                                if key == ko:
                                    ki = ka + "__" + key
                                    seqe = va[0] + value[0]
                                    globals()["f_" + str(num)][ki] = [seqe]
                            
                            else:    
                                ki = ka + "__" + key
                                seqe = va[0] + value[0]
                                globals()["f_" + str(num)][ki] = [seqe]
                           
                elif num == 0: # najprej rabim vse incializirati
                    if z.startswith("A2"):
                        u = "A2"
                    elif z.startswith("kozak"):
                        u = "kozak"
                    else:
                        u = z

                    for key,value in globals()[u + "_dic"].items():
                        if z.startswith("A2"):
                            ko = A2_list[A2_zap.index(z)]
                            if key == ko:
                                globals()["f_" + str(num)][key] = [value[1].reverse_complement()]
                                
                        elif z.startswith("kozak"):
                            ko = kozak_list[kozak_zap.index(z)]
                            if key == ko:
                                globals()["f_" + str(num)][key] = [value[1].reverse_complement()]
                                
                        else:        
                            globals()["f_" + str(num)][key] = [value[1].reverse_complement()]
                num = num + 1
            forward__primer.update(globals()["f_" + str(num-1)])


izlocen_reverse = delcki[zapzap[0]] # delcki prvega clena, tistega, ki ga izločim
rpi = izlocen_reverse.index(zapzap[0]) # v bistvu pridobim celotne še ne skrajšane sekvence biodela in tistega kki mu je predhoden za forward primer,
#posledično je večina reverse primerjev enakih reverznim transkriptom forward primerjev. S to metodo moramo le izločiti reverse primer prvega forward primerja in dodati končni reverse primer, ki nima forward komplementa.

irev = []
for i in izlocen_reverse: #izločim nepotreben šum (številke na koncu A2 in kozakov) (Slava koda?)
    if i.startswith("kozak"):
        irev.append("kozak")
    elif i.startswith("A2"):
        irev.append("A2")
    else:
        irev.append(i)

izrev = irev[:rpi+1] # izločamo prvi forward primer za pridobivanje reverse primerjev z reverznimi komplementi

rev = {}

for i in izrev: #naredi njižnico z vsemi kombinacijami uporabljenih delov za neuporabni komplement (vmes se koda konča, slaba koda, verjetno bi to lahko rešil že na začetku določanja primerjev)
    for key in globals()[i+"_dic"].keys():
        if not i in rev:
            rev[i] =[key]
        else:
            rev[i].append(key)


neuporabni_revpr = []

for key in forward__primer.keys():
    check = True
    count = 0
    if check is True:
        lista = key.split("__") # če ime ustreza zaporedju, ki ga je določil neuporabni => izloči
        for i in izrev:
            o = izrev.index(i)
            try:
                if lista[o] in rev[i]:
                    check = True
                    count = count + 1
                else:
                    check = False
            except:
                check = False
                
    if count == len(izrev):
        neuporabni_revpr.append(key)
        
reverse__primer = {}
for key, value in forward__primer.items():
    if not key in neuporabni_revpr:
        rki = key.split("__")
        rki.reverse()
        reverse_key = "__".join(rki)
        reverse_seq = value[0].reverse_complement()
        reverse__primer[reverse_key] = [reverse_seq]

dodan_reverse = delcki[zapzap[-1]]
drp = dodan_reverse.index(zapzap[-1])
dorev = []
for i in dodan_reverse:
    if i.startswith("kozak"):
        dorev.append("kozak")
    elif i.startswith("A2"):
        dorev.append("A2")
    else:
        dorev.append(i)
drev = dorev[drp:]

u_rev = {}
for i in drev:
    for key in globals()[i+"_dic"].keys():
        if not i in u_rev:
            u_rev[i] =[key]
        else:
            u_rev[i].append(key)

num = 0
drev.reverse()

for z in drev: #dodaj še zadnji reverzni primer, slaba koda? verjetno vseeno hitreje reverzni komplement kot tole še enkrat
    globals()["r_" + str(num)] = {}
    
    if not num == 0:
        for ka, va in globals()["r_" + str(num-1)].items():
            if z.startswith("A2"):
                u = "A2"
            elif z.startswith("kozak"):
                u = "kozak"
            else:
                u = z            
            for key,value in globals()[u + "_dic"].items():
                if z.startswith("A2"):
                    ko = A2_list[A2_zap.index(z)]
                    if key == ko:
                        ki = ka + "__" + key
                        seqe = va[0] + value[1]
                        globals()["f_" + str(num)][ki] = [seqe]
                elif z.startswith("kozak"):
                    ko = kozak_list[kozak_zap.index(z)]
                    if key == ko:
                        ki = ka + "__" + key
                        seqe = va[0] + value[1]
                        globals()["f_" + str(num)][ki] = [seqe]
                else:
                    ki = ka + "__" + key
                    seqe = va[0] + value[1]
                
                globals()["r_" + str(num)][ki] = [seqe]
                
    elif num == 0:
        if z.startswith("A2"):
            u = "A2"
        elif z.startswith("kozak"):
            u = "kozak"
        else:
            u = z
        
        for key,value in globals()[u + "_dic"].items():
            if z.startswith("A2"):
                ko = A2_list[A2_zap.index(z)]
                if key == ko:
                    globals()["f_" + str(num)][key] = [value[0].reverse_complement()]
            elif z.startswith("kozak"):
                ko = kozak_list[kozak_zap.index(z)]
                if key == ko:
                    globals()["f_" + str(num)][key] = [value[0].reverse_complement()]
            else:
                globals()["r_" + str(num)][key] = [value[0].reverse_complement()]
    num = num + 1
reverse__primer.update(globals()["r_" + str(num-1)])
#______________________________________________________________________________                  


assembly_dic = {}
koncni_assembly = [] # prepiši cel assembly + dodaj ustrezne restrikcijske encime na koncih -> napiši, tako da bo naerdilo genebank sekvenco z anotacijami
with open (os.path.join(asembly,"asembly.fasta"),"r") as rd, open(os.path.join(Rezultati,"asembly.fasta"),"a") as wr:
    for record in SeqIO.parse(rd,"fasta"):
        for key, value in dejanska_u.items():
            if record.id == key:
               aa = "__".join(value)
               for k, val in zp.items():
                   if aa == k:
                       res_encim = k.split("__")
                       ime = res_encim[0] + "__" + record.id + "__" + res_encim[1]
                       wr.write(">" + ime + "\n")
                       seq = Seq(zac_rez) + Seq(val[0]) + record.seq + Seq(val[1]) + Seq(kon_rez).reverse_complement()
                       wr.write(str(seq)+ "\n")
                       koncni_assembly.append(ime)
                       assembly_dic[ime] = seq


#______________________________________________________________________________
#primerji 2A

vrste =[]
with open(os.path.join(primer_input,"primer3inpt.fasta"),"r") as rd:
    for record in SeqIO.parse(rd,"fasta"):
        if record.id.startswith(("sp","tr")):
            vrstat = record.id.split("|")[0]
            vrsta = vrstat.split("_")[1]
            if not vrsta in vrste:
                vrste.append(vrsta)
                globals()[vrsta + "__CDS"] = {}
           
with open(os.path.join(primer_input,"primer3inpt.fasta"),"r") as rd:
    for record in SeqIO.parse(rd,"fasta"):
        if record.id.startswith(("sp","tr")): #potrebno ker Primer3 bere prvi heading kot komentar, zato vključen še en heading (">DaNeCrkne")
            vrstat = record.id.split("|")[0]
            vrsta = vrstat.split("_")[1]
            ime  = record.id.split("|")[1]
            globals()[vrsta + "__CDS"][ime] = record.seq

fseq_list = {}
rseq_list = {}
A2_primerji_celi = {}
for clen in koncni_primerji_z_2A_podaljski: 
    c_list = clen.split("_")
    smer = c_list[0] # 2A podaljški na forward primerjih na zadnjih členih, na reverse primerjih na sprednjem členu
    vrsta = c_list[1]
    
    print(clen)

    for key, value in reverse__primer.items():
        kljuc = key.split("__")[-1]
        try:
            if kljuc in globals()[vrsta + "_dic"]:    
                for ky, vale in forward__primer.items():
                    klju = ky.split("__")[-1]
                    
                    if klju == kljuc:
                        for k in koncni_assembly:
                            zapor = k.split("__")
                            
                            try:
                                a = zapor.index(kljuc)
                                if (smer == "r" and zapor[a + 1] in A2_list) or (smer == "f" and ((zapor[a - 2] in A2_list) or (zapor[a - 1] in A2_list))):
                                    for kay, velue in globals()[vrsta + "__CDS"].items(): #če je 2A sekvenca 2 enoti za/pred dano enoto, potem iz tega skonstruiraj primer
                                        for kiy, vula in globals()[vrsta + "_dic"].items():
                                            if kiy == kljuc and kay == kljuc:

                                                dolzina_f_pr = len(vula[0]) 
                                                dolzina_r_pr = len(vula[1])

                                                sekvenca = velue[dolzina_f_pr: -dolzina_r_pr]
                                                sekvenca = vale[0] + sekvenca + value[0].reverse_complement() #celoten primer bi lahko pokrival dolžino od začetka ene enote pa do konca druge

                                                ind = zapzap.index(vrsta)
                                                                                                
                                                if smer == "r":   
                                                    inde = zapzap[ind + 1]
                                                    for kz,vz in globals()[inde + "_dic"].items():
                                                        if kz in zapor:
                                                            sekvenca = sekvenca[:-len(vz[0])]
                                                            rseq_list[kljuc] = [Seq(sekvenca)]
                                                            rseq_list[kljuc].extend([Seq(vula[0])])
                                                            rseq_list[kljuc].extend([Seq(vula[1])])
                                                                                                            
                                                
                                                elif smer == "f":
                                                    inde = zapzap[ind - 1]
                                                    for kz,vz in globals()[inde + "_dic"].items():
                                                        if kz in zapor:
                                                            sekvenca = sekvenca[len(vz[1]):]
                                                            fseq_list[kljuc] = [Seq(sekvenca)]
                                                            fseq_list[kljuc].extend([Seq(vula[0])])
                                                            fseq_list[kljuc].extend([Seq(vula[1])])
                                                                                                            
                                                
                                                
                                                            
                                                            
                                                            
                                                            
                            except:
                                continue
        except:
            continue
        

f_za_pri_2 = {}
r_za_pri_2 = {}

f2A = {}
r2A = {}
fr2A = {}
rf2A = {}

primerji_2A_f = {}
primerji_2A_r = {}
with open(os.path.join(Rezultati, "A2_primerji.fasta"),"a") as primer: 
    #najdi mi primerje, ki podaljšajo sekvenco tako, da vežejo obstoječe konce in uvedejo zadosti velik overlap za Gibsona
    for key,value in rseq_list.items():
        rpr_v_okv = value[2].reverse_complement()
        reverse_primer_start = value[0].find(rpr_v_okv) #najdi binding primerja 1 na  reverse primerju (celotna sekvenca prejšnji CDS + 2A + Kozak + iskani CDS)
        
        fwd_primer_start = value[0].find(value[1]) #najdi binding primerja1 na forward primerju(celotna sekvenca pri CDS + 2A + Kozak + naslednji CDS)
        #forward primer prejšnjega predstavlja kje se reverse primer naslednjega konča, KODA VERJETNO ZARADI TEGA DELA NE BO DELAL; ČE BOSTA VEČ kot 2 2A PRIMERJA (ne privzema, da lahko 2A primerja zrasteta it obeh strani)
        pet_pozicija_rpr = reverse_primer_start + 59 #na začetku smo se ravnali, kot da je maksimalna velikost oligov 60 nt.
        # bolje bi bilo če bi lokacijo primerjev izpisal, ko jih naredim, namesto da jih potem iščem.
        #maksimum dolžina drugega nukleotida je bila 60 nt, podaljšamo ga lahko za max 60 nt, minimalni overlap je 18 nt (nastavljen zgoraj)
        
        tri_pozicija_rpr = len(value[0]) - 61
        print("#####")
        seq_args = {
            'SEQUENCE_ID': "r__A2__" + key,
            'SEQUENCE_TEMPLATE': value[0],
            'SEQUENCE_INCLUDED_REGION': [0, (pet_pozicija_rpr)],
            'SEQUENCE_FORCE_RIGHT_END': [tri_pozicija_rpr + 1],
            'SEQUENCE_FORCE_LEFT_END': [fwd_primer_start + len(value[1]) - 1],
        }
        
        primer_params = {        
        'PRIMER_OPT_SIZE': 18,
        'PRIMER_MIN_SIZE': 15,
        'PRIMER_MAX_SIZE': (pet_pozicija_rpr) - (tri_pozicija_rpr) + 1,
        'PRIMER_MIN_TM': 54.5,
        'PRIMER_OPT_TM': 60,
        'PRIMER_MAX_TM': 75,
        'PRIMER_NUM_RETURN': st_pr,
        'PRIMER_MIN_GC': 0,
        'PRIMER_MAX_GC': 100,                        
        'PRIMER_PRODUCT_SIZE_RANGE': [tri_pozicija_rpr - (fwd_primer_start + len(value[1])), tri_pozicija_rpr],
        'PRIMER_EXPLAIN_FLAG': 1,
        
        'PRIMER_DNA_CONC': 500.0,
        'PRIMER_DNTP_CONC': 0.8
        }
                            
        result = bindings.design_primers(primer_params,seq_args)
        
        try:
            print(result["PRIMER_LEFT_0_SEQUENCE"])
            print(result["PRIMER_RIGHT_0_SEQUENCE"])
             
        except:
            print(".................................")
            print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
            print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
            print(".................................")
            continue   
        
        
        for i in range(st_pr):
            try:
                a = Seq(result[f'PRIMER_LEFT_{i}_SEQUENCE']) 
                b = Seq(result[f'PRIMER_RIGHT_{i}_SEQUENCE'])
                fTm = str(result[f'PRIMER_LEFT_{i}_TM'])
                rTm = str(result[f'PRIMER_RIGHT_{i}_TM'])
                
                smak = key.strip()
                
                primer.write(f">f_{i}_"+  smak  +"__"+ fTm + "\n")
                primer.write(str(a)+ "\n" )
                primer.write(f">r_{i}_"+  smak  +"__"+ rTm + "\n")
                primer.write(str(b)+ "\n" )
                
                
                r_za_pri_2[smak +"__"+ str(i)] = len(value[0])-reverse_primer_start - 60 + len(Seq(b))   #(len(value[0][reverse_primer_start:-60]) + len(Seq(b))) 


                r2A[smak + "__" + rTm] = b
                rf2A[smak+ "__" + fTm] = a
                
                primerji_2A_r[key] = [a,b,fTm,rTm]
            except:
                continue
# izpiše primer z maksimalno dolžino (60 nt)           
#______________________________________________________________________________  
    for key,value in fseq_list.items():
        rpr_v_okv = value[2].reverse_complement()
        reverse_primer_start = value[0].find(rpr_v_okv)
        
        fwd_primer_start = value[0].find(value[1])
        
        pet_pozicija_fpr = fwd_primer_start - (60 - len(value[1])) - 1
        tri_pozicija_fpr = 60
        max_size = (tri_pozicija_fpr) - (pet_pozicija_fpr) + 1

        seq_args = {
            'SEQUENCE_ID': "r__A2__" + key,
            'SEQUENCE_TEMPLATE': value[0],
            'SEQUENCE_INCLUDED_REGION': [pet_pozicija_fpr, len(value[0]) - pet_pozicija_fpr],
            'SEQUENCE_FORCE_RIGHT_END': [reverse_primer_start],
            'SEQUENCE_FORCE_LEFT_END': [tri_pozicija_fpr-1]
        }
        
        primer_params = {        
        'PRIMER_OPT_SIZE': 18,
        'PRIMER_MIN_SIZE': 15,
        'PRIMER_MAX_SIZE': max_size,
        'PRIMER_MIN_TM': 50,
        'PRIMER_OPT_TM': 60,
        'PRIMER_MAX_TM': 75,
        'PRIMER_NUM_RETURN': st_pr,
        'PRIMER_MIN_GC': 0,
        'PRIMER_MAX_GC': 100,                        
        'PRIMER_PRODUCT_SIZE_RANGE': [reverse_primer_start - tri_pozicija_fpr, len(value[0]) - pet_pozicija_fpr],
        'PRIMER_PAIR_MAX_DIFF_TM' : 5.0,
        'PRIMER_EXPLAIN_FLAG': 1,
        
        'PRIMER_DNA_CONC': 500.0,
        'PRIMER_DNTP_CONC': 0.8
        
        }
                            
        result = bindings.design_primers(primer_params,seq_args)
        
        try:
            print(result["PRIMER_LEFT_0_SEQUENCE"])
            print(result["PRIMER_RIGHT_0_SEQUENCE"])
             
        except:
            print(".................................")
            print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
            print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
            print(".................................")
            continue   
        
        
        for i in range(st_pr):
            try:
                a = Seq(result[f'PRIMER_LEFT_{i}_SEQUENCE']) 
                b = Seq(result[f'PRIMER_RIGHT_{i}_SEQUENCE'])
                fTm = str(result[f'PRIMER_LEFT_{i}_TM'])
                rTm = str(result[f'PRIMER_RIGHT_{i}_TM'])
                
                smak = key.strip()
                
                primer.write(f">f_{i}_"+  smak  +"__"+ fTm + "\n")
                primer.write(str(a)+ "\n" )
                primer.write(f">r_{i}_"+  smak  +"__"+ rTm + "\n")
                primer.write(str(b)+ "\n" )
                                
                f_za_pri_2[smak +"__"+ str(i)] = len(Seq(a)) #+ (fwd_primer_start - 60)

                f2A[smak + "__" + fTm] = a
                fr2A[smak + "__" + rTm] = b
                
                primerji_2A_f[key] = [a,b,fTm,rTm]
                
            except:
                continue

# izpiše primer z maksimalno dolžino (60 nt) 

#______________________________________________________________________________
Da_se_seq_ne_ponavlja = [] # ker sem določal za forward in reverse primerje posebej, moram zdaj izložiti tiste, ki so isti, SLABA KODA, bolje bi bilo samo narediti isto kot prej, 
#predhodno v reverse primerje postaviti le en primer, tisti, ki ga forward primerji ne zajemajo
forw = {}

for o in koncni_primerji_z_2A_podaljski: # filtriram le tiste na katere daje, 2A podaljške, Slaba koda, to bi moral narediti že prej, zdaj primer 3 za brezveze gleda za primerje zadev, ki jih ne rabijo
    if o.startswith("f"):
        a = o.split("_")[1]
        for z in zapzap:
            if z == a:
                for key in globals()[z + "_dic"].keys(): # daj mi na listo vse komponente razreda (ICS, IPL,SAMT), ki jih imajo koncni f2A primerji
                    if z not in forw:
                        forw[z] = [key]
                    else:
                        forw[z].append(key)

rew = {}
for o in koncni_primerji_z_2A_podaljski:
    if o.startswith("r"):
        a = o.split("_")[1]
        for z in zapzap:
            if z == a:
                for key in globals()[z + "_dic"].keys():# daj mi na listo vse komponente razreda (ICS, IPL,SAMT), ki jih imajo koncni r2A primerji
                    if z not in rew:
                        rew[z] = [key]
                    else:
                        rew[z].append(key)

koncni_sestavljeni_f_primerji = {}
koncni_sestavljeni_r_primerji = {}
primer_pari = {}

with open(os.path.join(Rezultati,"primerji.fasta"),"a") as wr:
    for key, value in forward__primer.items():        
        mezam = key.split("__")[-1]
        for ki, valu in CDS_primer.items():
            ide = ki.split("|")[1]

            if mezam == ide and ki.startswith("forward"):
                Tmf = ki.split("__")[1]
                
                val = value[0]
                
                flag = 0
                
                for ka, va in forw.items():
                    if flag == 0:
                        if ide in va:
                            flag = 1
                
                if flag == 1:
                    for k, v in f_za_pri_2.items():
                        if k.split("__")[0] == ide:
                            for muk,vuk in f2A.items():
                                njeheh = muk.split("__")[0]
                                if njeheh == ide:
                                    
                                    lokacija = val.find(vuk)
                                    
                                    #velikost =(v+len(valu))
                                    
                                    f = val[lokacija:]
                                    if f not in Da_se_seq_ne_ponavlja:
                                        wr.write(">" + "f_" + mezam + "_____" + key + "__"+ Tmf + "\n")
                                        wr.write(str(f) + "\n")
                                        koncni_sestavljeni_f_primerji[mezam + "__" + Tmf + "__" + key] = f
                                        Da_se_seq_ne_ponavlja.append(f)
                
                                        primer_pari[mezam] = [f]
                                
                else:
                    for ko, vo in zp.items():

                        for i in vo:
                            if i in val:
                                f = zac_rez + val
                                #f = val
                                break   
                            else:
                                
                                f = val[-(podaljski+len(valu)):]
                                
                    
                    if f not in Da_se_seq_ne_ponavlja:
                        wr.write(">" + "f_" + mezam + "_____" + key + "__"+ Tmf + "\n")
                        wr.write(str(f) + "\n")
                        if "XbaI" in key.split("__")[0]:
                            pod_f = zac_rez + f
                            koncni_sestavljeni_f_primerji[mezam + "__" + Tmf + "__" + key] = pod_f
                            Da_se_seq_ne_ponavlja.append(pod_f)
                            primer_pari[mezam] = [pod_f]
                        else:
                            
                            koncni_sestavljeni_f_primerji[mezam + "__" + Tmf + "__" + key] = f
                            Da_se_seq_ne_ponavlja.append(f)
                            primer_pari[mezam] = [f]
        
        for ky, vale in reverse__primer.items():
            sezam = ky.split("__")[-1]
            if mezam == sezam:
                for ki, valu in CDS_primer.items():
                    ide = ki.split("|")[1]
                    if mezam == ide and ki.startswith("reverse"):

                        Tmr = ki.split("__")[1]
                        
                        val = vale[0]
                        fleg = 0

                        
                        for ka, va in rew.items():
                            if fleg == 0:
                                if ide in va:
                                    fleg = 1
                        
                        if fleg == 1:
                            for k, v in r_za_pri_2.items():
                                if k.split("__")[0] == ide:
                                    
                                    velikost = len(vale)-v -1
                                    r = val[velikost:]
                                    
                                    if r not in Da_se_seq_ne_ponavlja:
                                        wr.write(">" + "r_" + mezam + "_____" +key + "__"+ Tmr + "\n")
                                        wr.write(str(r) + "\n")
                                        koncni_sestavljeni_r_primerji[mezam + "__" + Tmr + "__" + ky] = r
                                        Da_se_seq_ne_ponavlja.append(r)
                                        primer_pari[mezam].append(r)
                                    
                                    
                        else:
                            for ku, vu in zp.items():
                                for e in vu:
                                    if e in val:
                                        r = kon_rez + val
                                        #r= val
                                        break
                                    else:
                                        # zakaj pa imam to samo na reverse primerjih?
                                        r = val[-(podaljski+len(valu)):]
                                        
                                        rewers = r[:podaljski]
                                        Tm_Gibson = 0
                                        
                                        Tm_Gibson = calc_tm(str(rewers),dntp_conc = 0.8)
                                        podaljsaj = podaljski
                                        
                                        while Tm_Gibson <= 49:
                                            podaljsaj = podaljsaj + 1 #podaljšuj dokler ne prideš do zadostnega Gibsona
                                            r = val[-(podaljsaj + len(valu)):]
                                            rewers = r[:podaljsaj]
                                            
                                            Tm_Gibson = calc_tm(str(rewers), dna_conc = 25.0, dntp_conc = 0.8)
                                        

                            if r not in Da_se_seq_ne_ponavlja:
                                wr.write(">" + "r_" + mezam + "_____" + ky + "__" + Tmr + "\n")
                                wr.write(str(r) + "\n")
                                
                                if "XhoI" in key.split("__")[-1]:
                                    pod_r = r + kon_rez
                                    koncni_sestavljeni_r_primerji[mezam + "__" + Tmr + "__" + ky] = pod_r
                                    Da_se_seq_ne_ponavlja.append(pod_r)
                                    primer_pari[mezam].append(pod_r)
                                else:
                                    koncni_sestavljeni_r_primerji[mezam + "__" + Tmr + "__" + ky] = r
                                    Da_se_seq_ne_ponavlja.append(r)
                                    primer_pari[mezam].append(r)
                                    

#___________________________Genbank_assembly___________________________________        

STOP_kodoni = {Seq("TAG"), Seq("TAA"), Seq("TGA")}

with open(os.path.join(obd_seq, "zdr_u_seq.fasta"),"r") as rd:
    for record in SeqIO.parse(rd, "fasta"):
        a = record.id.split("|")
        vrsta = (a[0]).split("_")[1]
        ime = a[1]
        if vrsta in zapzap:
            try:
                globals()[vrsta + "_sekvence"][ime] = record.seq
            except:
                globals()[vrsta + "_sekvence"] = {}
                globals()[vrsta + "_sekvence"][ime] = record.seq
                


ne2a_dic = {}

forward_podaljšek = "AAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACCCCGGAT"
reverse_podaljšek = "TCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACAT"
f_Gibson_overalap =Seq(forward_podaljšek + "T")
r_Gibson_overlap = Seq ("G" + reverse_podaljšek).reverse_complement()


cel_csv_input = [[],["ICS_GFP","fPrimer","rPrimer","fTm","rTm","fHp","dh","ds","rHp","dh","ds","f_Dimer","dh","ds","r_Dimer","dh","ds","HeteroDimer","dh","ds"]]


def uporabne_tm_info (f_pm,r_pm):
                fHp = calc_hairpin(str(f_pm),
                                       dntp_conc = 0.2,
                                       output_structure=True,
                                       dna_conc = 500.0,
                                       mv_conc= 50,
                                       dv_conc= 1.5)
                
                rHp = calc_hairpin(str(r_pm),
                                       dntp_conc = 0.2,
                                       output_structure=True,
                                       dna_conc = 500.0,
                                       mv_conc= 50,
                                       dv_conc= 1.5)
            
                dimer_f = calc_homodimer(str(f_pm),
                                       dntp_conc = 0.2,
                                       dna_conc = 500.0,
                                       mv_conc= 50,
                                       dv_conc= 1.5)
                
                dimer_r =calc_homodimer(str(r_pm),
                                       output_structure=True,
                                       dntp_conc = 0.2,
                                       dna_conc = 500.0,
                                       mv_conc= 50,
                                       dv_conc= 1.5)
                
                Hetero = calc_heterodimer(str(f_pm),str(r_pm),
                                          output_structure=True,
                                          dntp_conc = 0.2,
                                          dna_conc = 500.0,
                                          mv_conc= 50,
                                          dv_conc= 1.5)
                
                
                print("forward_hairpin")
                print(fHp)
                print("reverse_hairpin")
                print(rHp)
                print("forward_dimer")
                print(dimer_f)
                print("reverse_dimer")
                print(dimer_r)
                
                for line in dimer_r.ascii_structure_lines:
                        art = line.split("\t")[1]
                        print(art)
                        
                for line in Hetero.ascii_structure_lines:
                        art = line.split("\t")[1]
                        print(art)
                        
                print("heterodimer")
                print(Hetero)
                
                globals()["csv_primer_data"] = [str(fHp.tm),str(fHp.dh),str(fHp.ds),str(rHp.tm),str(rHp.dh),str(rHp.ds),str(dimer_f.tm),str(dimer_f.dh),str(dimer_f.ds),str(dimer_r.tm),str(dimer_r.dh),str(dimer_r.ds),str(Hetero.tm),str(Hetero.dh),str(Hetero.ds)]
                    
                    
                    
A2_primerji_lista = {}
for key, value in assembly_dic.items():
    with open (os.path.join(Genbank_rezultati,key + ".gb"),"a") as wr:
        rez = []

        
        #modifikacija_za_gibsona_______________________________________________
        
        #celotni primerji:
            
        print("______________________CELOTNI_PRIMERJI__________________________")
        
        with open(os.path.join(Rezultati,"primerji.fasta"),"a") as primer:
            
            seq_args = {
                'SEQUENCE_ID': "r__A2__" + "X82644yopt",
                'SEQUENCE_TEMPLATE': value,
                'PRIMER_TASK':  'pick_pcr_primers',
                
                'SEQUENCE_FORCE_LEFT_START': [0],
                #'SEQUENCE_PRIMER': "TCTAGAATAACCATGTCTCAATCACT"
                #'SEQUENCE_FORCE_RIGHT_START': [len(value)]
            }
            
            primer_params = {
            #'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 0,
            'PRIMER_OPT_SIZE': 18,
            'PRIMER_MIN_SIZE': 15,
            #'PRIMER_MAX_SIZE': max_size,
            'PRIMER_MIN_TM': 50,
            'PRIMER_OPT_TM': 60,
            'PRIMER_MAX_TM': 72,
            'PRIMER_NUM_RETURN': st_pr,
            'PRIMER_MIN_GC': 0,
            'PRIMER_MAX_GC': 100,     
            'PRIMER_PRODUCT_SIZE_RANGE':[len(value), len(value)],          
            'PRIMER_PAIR_MAX_DIFF_TM' : 5.0,
            'PRIMER_EXPLAIN_FLAG': 1,
            
            #'PRIMER_MAX_HAIRPIN_TH ': 1000000,
            #'PRIMER_MAX_SELF_ANY_TH':   10000000,
            #'PRIMER_MAX_SELF_END_TH': 10000000,
            #'PRIMER_MAX_HAIRPIN_TH': 10000000,
            'PRIMER_DNA_CONC': 500.0,
            'PRIMER_DNTP_CONC': 0.8
            
            }
                                
            result = bindings.design_primers(primer_params,seq_args)
            print("########")
            try:
                print(result["PRIMER_LEFT_0_SEQUENCE"])
                print(result["PRIMER_RIGHT_0_SEQUENCE"])
                 
            except:
                print(".................................")
                print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
                print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
                print(".................................")
                continue   
            
            
            for i in range(st_pr):
                try:
                    cel_f = Seq(result[f'PRIMER_LEFT_{i}_SEQUENCE'])
                    cel_r = Seq(result[f'PRIMER_RIGHT_{i}_SEQUENCE'])
                    fTm = str(result[f'PRIMER_LEFT_{i}_TM'])
                    rTm = str(result[f'PRIMER_RIGHT_{i}_TM'])
                                        
                    primer.write(f">f_{i}_"+  "celotni_primer"  +"__"+ fTm + "\n")
                    primer.write(str(cel_f)+ "\n" )
                    primer.write(f">r_{i}_"+  "Celotni_primer"  +"__"+ rTm + "\n")
                    primer.write(str(cel_r)+ "\n" )
            
                except:
                    continue
            
            
            print("________cel_______")
            
            str_cel_f = result[f'PRIMER_LEFT_{i}_SEQUENCE']
            str_cel_r = result[f'PRIMER_RIGHT_{i}_SEQUENCE']
            uporabne_tm_info(str_cel_f,str_cel_r)
            
            celotni_csv_podatki = ["cel_primer_brez_podaljškov_in_brez_his",str_cel_f, str_cel_r, fTm, rTm]
            celotni_csv_podatki.extend(csv_primer_data)
            cel_csv_input.append(celotni_csv_podatki)
        
        print()
        
        Gibson_f_Tm = 100.0
        while Gibson_f_Tm >=  62.0:
            Gibson_f_Tm  = calc_tm(str(f_Gibson_overalap),
                                   dntp_conc = 0.2,
                                   dna_conc = 500.0,
                                   mv_conc= 50,
                                   dv_conc= 1.5)
            f_Gibson_overalap = f_Gibson_overalap[1:]
        
        print("______________________OVERLAP_TM___________________________")
        
        
        print(f_Gibson_overalap)
        print(Gibson_f_Tm)
        
        
        Gibson_r_Tm = 100.0
        while Gibson_r_Tm >=  62.0:
            Gibson_r_Tm  = calc_tm(str(r_Gibson_overlap),
                                   dntp_conc = 0.2,
                                   dna_conc = 500.0,
                                   mv_conc= 50,
                                   dv_conc= 1.5)
            r_Gibson_overlap = r_Gibson_overlap[1:]
        print(r_Gibson_overlap)
        print(Gibson_r_Tm)
        print()
        
        print("______________________sekvence_cel_primer__________________________")
        print(cel_f)
        print(cel_r)
        
        podaljšan_cel_f =f_Gibson_overalap[:-1] + cel_f
        podaljšan_cel_r = r_Gibson_overlap[:-1] + cel_r 
        
        str_pod_cel_f = str(podaljšan_cel_f)
        str_podaljšan_cel_r = str(podaljšan_cel_r)
        
        print(podaljšan_cel_f)
        print(podaljšan_cel_r)
        
        uporabne_tm_info(str_pod_cel_f,str_podaljšan_cel_r)
        
        celotni_csv_podatki = ["cel_primer_z_vektor_podaljški_brez_his_taga",str_pod_cel_f, str_podaljšan_cel_r, fTm, rTm]
        celotni_csv_podatki.extend(csv_primer_data)
        
        cel_csv_input.append(celotni_csv_podatki)
        
        podaljški_F = podaljšan_cel_f[:podaljšan_cel_f.find(cel_f)]
        podaljški_R = podaljšan_cel_r[:podaljšan_cel_r.find(cel_r)]
        
        value = podaljški_F + value + podaljški_R.reverse_complement()

        konstrukt = SeqRecord(
            value,
            id = "Popolni_konstrukt",
            name="Popolni_konstrukt",
            description="koncni konstrukt, ki ga bo dalo Gibsonovo kloniranje",
            annotations={"molecule_type": "DNA"}
        )
        
        #konstrukt.features = []
        f_vektor = SeqFeature(
            FeatureLocation(value.find(podaljški_F),value.find(podaljški_F) + len(podaljški_F)),
            type = "primer_bind",
            qualifiers = {"label" : ["f_primer_podaljsek_za_celo_sekvenco(GIBSON_VEKTOR_OVERLAP)"]})
        r_vektor = SeqFeature(
            FeatureLocation(value.find(podaljški_R.reverse_complement()),value.find(podaljški_R.reverse_complement()) + len(podaljški_R),strand = -1),
            type = "primer_bind",
            qualifiers = {"label" : ["f_primer_podaljsek_za_celo_sekvenco(GIBSON_VEKTOR_OVERLAP)"]})
        
        f_ful_primer = SeqFeature(
            FeatureLocation(value.find(podaljšan_cel_f),value.find(podaljšan_cel_f) + len(podaljšan_cel_f)),
            type = "primer_bind",
            qualifiers = {"label" : ["f_cel_primer"]})
        r_ful_primer = SeqFeature(
            FeatureLocation(value.find(podaljšan_cel_r.reverse_complement()),value.find(podaljšan_cel_r.reverse_complement()) + len(podaljšan_cel_r),strand = -1),
            type = "primer_bind",
            qualifiers = {"label" : ["r_cel_primer"]})
        
        konstrukt.features = [f_vektor,r_vektor,f_ful_primer,r_ful_primer]
        
        aba = key.split("__")
        rez.append(aba[0])
        rez.append(aba[-1])

        for a in zapzap:
            for k, v in globals()[a + "_sekvence"].items():
                if k in key.split("__"):
                    
                    n = 0                    
                    sekv = v
                    while sekv[-3:] in STOP_kodoni:
                        sekv = sekv[:-3]
                        n = n + 1
                    
                    globals()[k] = SeqRecord(
                    v,
                    id = k,
                    name= a +"__" + k,
                    description="osnovna sekvenca " + a + ": " + k,
                    annotations={"molecule_type": "DNA"}
                    )
                    if v[0:3] in (Seq("ATG"), Seq("GTG")):
                        start_kodon = SeqFeature(
                            FeatureLocation(0,3),
                            type = "misc_feature",
                            qualifiers = {"label" : ["START_kodon"], "note": ["START kodon na osnovni sekvenci"]}
                            )
                    
                    STOP_kodon = SeqFeature(
                        FeatureLocation(len(v)-n*3, len(v)),
                        type = "misc_feature",
                        qualifiers = {"label" : ["STOP_kodoni"], "note": ["STOP kodoni na osnovni sekvenci"]}
                        )                    

                    for ko, vo in CDS_primer.items():
                        ab = ko.split("|")[1]
                        if ab == k:
                            if ko.startswith("f"):
                                f_mesto = v.find(vo)
                                f_kon_mesto = value.find(vo)
                                
                                f_kon_bind = SeqFeature(
                                    FeatureLocation(f_kon_mesto,f_kon_mesto + len(vo)),
                                    type = "primer_bind",
                                    qualifiers = {"label" : ["forward_primer_" + k], 
                                                  "note": ["vezavno mesto za forward primer " + ab +": " + k +" podaljski niso oznaceni"],
                                                  "Tm": [ko.split("__")[1] ]})
                                
                                f_primer_bind = SeqFeature(
                                    FeatureLocation(f_mesto,f_mesto + len(vo)),
                                    type = "primer_bind",
                                    qualifiers = {"label" : ["forward_primer_" + k], 
                                                  "note": ["vezavno mesto za forward primer " + ab +": " + k +" podaljski niso oznaceni"],
                                                  "Tm": [ko.split("__")[1] ]}
                                    )
                                konstrukt.features.append(f_kon_bind)

                            elif ko.startswith("r"):
                                r_mesto = v.find(vo.reverse_complement())
                                r_primer_bind = SeqFeature(
                                    FeatureLocation(r_mesto,r_mesto + len(vo),  strand = -1),
                                    type = "primer_bind",
                                    qualifiers = {"label" : ["reverse_primer_" + k],
                                                  "note": ["vezavno mesto za reverse primer " + ab +": " + k +" podaljski niso oznaceni"],
                                                  "Tm": [ko.split("__")[1] ]}
                                    )
                                r_kon_mesto = value.find(vo.reverse_complement())
                                r_kon_bind = SeqFeature(
                                    FeatureLocation(r_kon_mesto,r_kon_mesto + len(vo),  strand = -1),
                                    type = "primer_bind",
                                    qualifiers = {"label" : ["reverse_primer_" + k],
                                                  "note": ["vezavno mesto za reverse primer " + ab +": " + k +" podaljski niso oznaceni"],
                                                  "Tm": [ko.split("__")[1] ]}
                                    )
                                
                                CDS = SeqFeature(
                                    FeatureLocation(f_mesto,r_mesto + len (vo)),
                                    type = "CDS",
                                    qualifiers = {"label" : ["pomnozek_" + k], "note": ["del gena, ki ga bomo dobili po PCR reakciji " + ab +": " + k +" brez podaljskov"]}
                                    )
                                konstrukt.features.append(r_kon_bind)

                            


                    globals()[k].features = [start_kodon, STOP_kodon, CDS, f_primer_bind, r_primer_bind]
                    
                    SeqIO.write(globals()[k], wr, "genbank")
                    
                    
#____________________________sestavljene sekvence______________________________
                    for u, vu in koncni_sestavljeni_f_primerji.items():
                        
                        tat = u.split("__")
                        tut = tat[0]

                        if tut == k:
                            for e, ve in koncni_sestavljeni_r_primerji.items():
                                tit = e.split("__")
                                tet = tit[0]
                                if tet == k:
                                    for z, m in CDS_primer.items():
                                        if z.split("|")[1] == tut:
                                            if z.startswith("f") :
                                                dolzina_podaljska_f = len(vu) - len(m)
                                            if z.startswith("r"):
                                                dolzina_podaljska_r = len(ve) - len(m)
                                                
                                    with open(os.path.join(Rezultati,"primerji.csv"),"a", newline='') as neki:
                                        
                                        nin = csv.writer(neki)
                                        
                                        vstopni_podatki = ["tarca","reakcija", "forward primer sekvenca", "forward primer Tm", "reverse primer sekvenca", "reverse primer Tm", "pomnozek", "matrica"]
                                        data = []
                                        
                                        with open(os.path.join(Rezultati,"primerji.csv"),"r",newline='') as red:
                                            n = 0
                                            for line in red:
                                                n = n + 1
                                            if n == 0:
                                                data.append(vstopni_podatki)                              
                                        
                                        if a not in zapzap[-1]:
                                            sekvenca_po_prvem_pomnozevanju = vu[:dolzina_podaljska_f] + sekv[3:] + ve.reverse_complement()[len(ve)-dolzina_podaljska_r:]
                                        else:
                                            sekvenca_po_prvem_pomnozevanju = vu[:dolzina_podaljska_f] + v[3:] + ve.reverse_complement()[len(ve)-dolzina_podaljska_r:]

                                        with open(os.path.join(Rezultati,"prvi_pomnozek.fasta"),"a") as wr99: 
                                            wr99.write(">" + tut + "\n")
                                            wr99.write(str(sekvenca_po_prvem_pomnozevanju) + "\n")


                                        sekvenca_prvo_pomnozevanje = SeqRecord(
                                            sekvenca_po_prvem_pomnozevanju,
                                            id = tut,
                                            name= tut +"_podaljsan",
                                            description="produkt pomnozevanja " + a + ": " + tut  ,
                                            annotations={"molecule_type": "DNA"}
                                            )
                                        
                                        if a not in kok:#
                                            ne2a_dic[tut] = [vu]#
                                            ne2a_dic[tut].append(ve)#
                                            
                                        
                                        f_kon_pod_loc = value.find(vu)
                                        f_kon_podaljsek = SeqFeature(
                                            FeatureLocation(f_kon_pod_loc, f_kon_pod_loc+dolzina_podaljska_f),
                                            type = "misc_feature",
                                            qualifiers = {"label" : ["f_primer_podaljsek__" + tut], "note": ["forward podaljsek, ki smo ga uvedli s primerjem"], "Tm" : [u.split("__")[1]]}
                                            )
                                        konstrukt.features.append(f_kon_podaljsek)
                                        
                                        f_primer_podaljsek = SeqFeature(
                                            FeatureLocation(0, dolzina_podaljska_f),
                                            type = "misc_feature",
                                            qualifiers = {"label" : ["f_primer_podaljsek__" + tut], "note": ["forward podaljsek, ki smo ga uvedli s primerjem"], "Tm" : [u.split("__")[1]]}
                                            )
                                        
                                        r_primer_podaljsek = SeqFeature(
                                            FeatureLocation(len(sekvenca_po_prvem_pomnozevanju)-dolzina_podaljska_r, len(sekvenca_po_prvem_pomnozevanju), strand=-1),
                                            type = "misc_feature",
                                            qualifiers = {"label" : ["r_primer_podaljsek__" + tet], "note": ["reverse podaljsek, ki smo ga uvedli s primerjem"], "Tm" : [e.split("__")[1]]}
                                            )


                                        r_kon_pod_loc = value.find(ve.reverse_complement())
                                        
                                        r_kon_podaljsek = SeqFeature(
                                            FeatureLocation( r_kon_pod_loc + len(ve)-dolzina_podaljska_r, r_kon_pod_loc + len(ve),strand=-1),
                                            type = "misc_feature",
                                            qualifiers = {"label" : ["r_primer_podaljsek__" + tet], "note": ["reverse podaljsek, ki smo ga uvedli s primerjem"], "Tm" : [e.split("__")[1]]}
                                            )
                                        konstrukt.features.append(r_kon_podaljsek)

                                        
                                        f_primr_bind = SeqFeature(
                                            FeatureLocation(0,len(vu)),
                                            type = "primer_bind",
                                            qualifiers = {"label" : ["f_primer__" + tut], "note": ["forward_primer " + tut], "Tm" : [u.split("__")[1]]}
                                            )
                                        
                                        f_kon_bind = SeqFeature(
                                            FeatureLocation(f_kon_pod_loc,f_kon_pod_loc + len(vu)),
                                            type = "primer_bind",
                                            qualifiers = {"label" : ["f_primer__" + tut], "note": ["forward_primer " + tut], "Tm" : [u.split("__")[1]]}
                                            )
                                        
                                        konstrukt.features.append(f_kon_bind)
                                        
                                        
                                        r_primr_bind = SeqFeature(
                                            FeatureLocation(len(sekvenca_po_prvem_pomnozevanju)-len(ve),len(sekvenca_po_prvem_pomnozevanju),strand=-1),
                                            type = "primer_bind",
                                            qualifiers = {"label" : ["r_primer__" + tet], "note": ["reverse_primer " + tet], "Tm" : [e.split("__")[1]]}
                                            )
                                        
                                        r_kon_bind = SeqFeature(
                                            FeatureLocation(r_kon_pod_loc,r_kon_pod_loc + len(ve),strand=-1),
                                            type = "primer_bind",
                                            qualifiers = {"label" : ["r_primer__" + tet], "note": ["reverse_primer " + tet], "Tm" : [e.split("__")[1]]}
                                            )
                                        konstrukt.features.append(r_kon_bind)
                                        
                                        sekvenca_prvo_pomnozevanje.features=[f_primer_podaljsek, r_primer_podaljsek,f_primr_bind, r_primr_bind]
                                        
                                        f_data_pri = sekvenca_po_prvem_pomnozevanju[0:len(vu)]
                                        r_data_pri = sekvenca_po_prvem_pomnozevanju[len(sekvenca_po_prvem_pomnozevanju)-len(ve):].reverse_complement()
                                        CDS_original = v
                                        f_Tm = u.split("__")[1]
                                        r_Tm = e.split("__")[1]
                                        
                                        prvo_pom_data = [k,"prvo_pomnozevanje",f_data_pri , f_Tm , r_data_pri, r_Tm, sekvenca_po_prvem_pomnozevanju, CDS_original]
                                        data.append(prvo_pom_data)
                                        
                                        fg = True
                                        tot = []
                                        for c in tat[2:]:
                                            if not c in tot:
                                                tot.append(c)
                                        for c in reversed(tit[2:]):
                                            if not c in tot:
                                                tot.append(c)
                                        print(tot)        
                                        zadnji_2A =""
                                        if "2A" in tot[-3]:
                                            zadnji_2A = tot[-3]
                                            
                                        for m in tot:
                                            if not "2A" in m and fg is True:
                                                if m in r_zap:
                                                    
                                                    k = r_zap.index(m)
                                                    seq_poisci = sekvence[k]
                                                    najdi = sekvenca_po_prvem_pomnozevanju.find(seq_poisci)
                                                    Res_mesto = SeqFeature(
                                                        FeatureLocation(najdi, najdi + len(seq_poisci)),
                                                        type = "misc_feature",
                                                        qualifiers = {"label" : [m  + "__restrikcjisko mesto"], "note": ["restrikcijsko mesto za insert v vektor"]})
                                                    sekvenca_prvo_pomnozevanje.features.append(Res_mesto)
                                                    
                                                elif "kozak" in m:
                                                    with open(os.path.join(kozak, "vsi_kozak.fasta"),"r") as rd1:
                                                        for rf in SeqIO.parse(rd1,"fasta"):
                                                            if rf.id == m:
                                                                paz = sekvenca_po_prvem_pomnozevanju.find(rf.seq[5:])
    
                                                                if paz -5 >=0:
                                                                    globals()[rf.id] = SeqFeature(
                                                                        FeatureLocation(paz-5, paz-5+len(rf.seq)),
                                                                        type = "misc_feature",
                                                                        qualifiers = {"label" : [rf.id], "note": ["kozak zaporedje"]})
                                                                    sekvenca_prvo_pomnozevanje.features.append(globals()[rf.id])
                                                                else:
                                                                    globals()[rf.id] = SeqFeature(
                                                                        FeatureLocation(0, paz + len(rf.seq)-5),
                                                                        type = "misc_feature",
                                                                        qualifiers = {"label" : [rf.id], "note": ["kozak zaporedje"]})
                                                                    sekvenca_prvo_pomnozevanje.features.append(globals()[rf.id])
                                                                    
                                                                
    
                                                
                                                elif not m in tot[0]:
                                                    with open(os.path.join(primer_input,"primer3inpt.fasta"),"r") as rd3:
                                                        for re in SeqIO.parse(rd3,"fasta"):
                                                            if not re.seq == "":
                                                                if re.id.split("|")[1] == m:
                                                                    vrstica = (re.id.split("|")[0]).split("_")[1]
                                                                    piz = sekvenca_po_prvem_pomnozevanju.find(re.seq[0:30])
                                                                    globals()[re.id] = SeqFeature(
                                                                        FeatureLocation(piz, piz+len(re.seq)),
                                                                        type = "CDS",
                                                                        qualifiers = {"label" : [re.id], "note": [vrstica]})
                                                                    sekvenca_prvo_pomnozevanje.features.append(globals()[re.id])
                                                                    break
                                            
                                            else:
                                                with open(os.path.join(dvaA,"A2.fasta"),"r") as rd2:
                                                    for rc in SeqIO.parse(rd2,"fasta"):
                                                        if rc.id == m:
                                                            if m in tot[1]:
                                                                with open(os.path.join(kozak, "vsi_kozak.fasta"),"r") as rd1:
                                                                    for rf in SeqIO.parse(rd1,"fasta"):
                                                                        if rf.id == tot[2]:
                                                                            paz = sekvenca_po_prvem_pomnozevanju.find(rf.seq[-podaljski:])
             
                                                                poz = paz #po 2A vedno kozak
                                                                if paz >=2:
                                                                    globals()[rc.id] = SeqFeature(
                                                                         FeatureLocation(0, poz),
                                                                         type = "misc_feature",
                                                                         qualifiers = {"label" : [rc.id], "note": ["2A zaporedje od " + str(len(rc.seq) -poz) + " do konca"]})
                                                                    sekvenca_prvo_pomnozevanje.features.append(globals()[rc.id])
                                                            else:
                                                                
                                                                poz = piz +len(re.seq) #pred vedno CDS
                                                                globals()[rc.id] = SeqFeature(
                                                                     FeatureLocation(poz, len(sekvenca_po_prvem_pomnozevanju)),
                                                                     type = "misc_feature",
                                                                     qualifiers = {"label" : [rc.id], "note": ["2A zaporedje od zacetka do " + str(len(sekvenca_po_prvem_pomnozevanju) - poz)]})
                                                                sekvenca_prvo_pomnozevanje.features.append(globals()[rc.id])
                                                            
                                                            
                                                            
                                                            if m == zadnji_2A:
                                                                pozicija = tot.index(m)
                                                                for ik in koncni_primerji_z_2A_podaljski:
                                                                    pop = ik.split("_")
                                                                    if pop[1] in zapzap:
                                                                        for klj in globals()[pop[1]+ "_dic"].keys():
                                                                            
                                                                            if klj not in A2_primerji_celi.keys():
          
                                                                                
                                                                                A2_primerji_celi[klj] = []
                                                                            if klj == tot[pozicija-1] and pop[0] == "r":
                                                                               for kw, vw in r2A.items():
                                                                                   if kw.split("__")[0] == tot[pozicija-1]:
                                                                                       mest = rc.seq.find(vw.reverse_complement())
                                                                                       for kp, vp in kozak_dic.items():
                                                                                           if kp == tot[pozicija+1]:
                                                                                               pravi_kozak = vp
                                                                                               break
                                                                                       
                                                                                       print(A2_primerji_celi[klj])
                                                                                       A2_pri_seq = vw.reverse_complement() + rc.seq[mest + len(vw):] + pravi_kozak[0]
                                                                                       
                                                                                       
                                                                                       if len(A2_pri_seq) >=61:
                                                                                           A2_pri_seq= A2_pri_seq[:60]
                                                                                       
                                                                                       
                                                                                       print("______________________2A_primerji__________________________")
                                                                                       print(klj + "__r")
                                                                                       
                                                                                       
                                                                                       if "0" in A2_primerji_celi[klj]:
                                                                                           A2_primerji_celi[klj][1] = A2_pri_seq
                                                                                       else:
                                                                                           A2_primerji_celi[klj] = ["0", A2_pri_seq]
                                                                                       
                                                                                      
                                                                                       
                                                                                       for kj, vlu in rseq_list.items():
                                                                                           print(kj)
                                                                                           for kk, vv in primerji_2A_r.items():
                                                                                               if kj == klj and kk == klj:
                                                                                                   fpr = vlu[0][:vlu[0].find(vlu[1]) + len(vlu[1])]
                                                                                                   
                                                                                                   fTM = vv[2]
                                                                                                   rTM = vv[3]
                                                                                                   
                                                                                       uporabne_tm_info(fpr,A2_pri_seq.reverse_complement())
                                                                                       celotni_csv_podatki = [klj + "__r" +"2A_primer in pripadajoči forward primer",str(fpr), str(A2_pri_seq.reverse_complement()), fTM, rTM]
                                                                                       celotni_csv_podatki.extend(csv_primer_data)
                                                                                       cel_csv_input.append(celotni_csv_podatki)
                                                                                       
                                                                                       print("_________Tm__________")
                                                                                       print(fTM)
                                                                                       print(rTM)
                                                                                       
                                                                                       
                                                                                       A2_sekvenca = SeqRecord(
                                                                                            A2_pri_seq,
                                                                                            id = rc.id,
                                                                                            name= rc.id + "_2A_primer",
                                                                                            description= "2A primer za drugo PCR reakcijo"  ,
                                                                                            annotations={"molecule_type": "DNA"}
                                                                                            )
                                                                                       
                                                                                       r_lok_A2 = value.find(A2_pri_seq)
                                                                                        
                                                                                       konstrukt_r2A = SeqFeature(
                                                                                           FeatureLocation(r_lok_A2,r_lok_A2 + len(A2_pri_seq), strand = -1),
                                                                                           type = "primer_bind",
                                                                                           qualifiers = {"label" : [rc.id + "_r_2"], "note": ["reverse primer druge PCR reakcije"], "Tm" : [kw.split("__")[1]]}
                                                                                           )
                                                                                       
                                                                                       konstrukt.features.append(konstrukt_r2A)
                                                                                       
                                                                                       A2_primer = SeqFeature(
                                                                                            FeatureLocation(0,len(vw), strand = -1),
                                                                                            type = "primer_bind",
                                                                                            qualifiers = {"label" : [rc.id + "_r_2"], "note": ["reverse primer druge PCR reakcije"], "Tm" : [kw.split("__")[1]]})
                                                                                       
                                                                                       
                                                                                       drugi_pomnozek = sekvenca_po_prvem_pomnozevanju[0: -len(vw)] + A2_pri_seq
                                                                                       for kz in rf2A.keys():
                                                                                           if kz.split("__")[0] == tut:
                                                                                               TM = kz.split("__")[1]
                                                                                               break
                                                                                               
                                                                                       drugo_pom_dat = [tut, "drugo_pomnozevanje", f_data_pri, TM, A2_pri_seq.reverse_complement(), kw.split("__")[1], drugi_pomnozek, sekvenca_po_prvem_pomnozevanju]
                                                                                       data.append(drugo_pom_dat)
                                                                            
                                                                            elif klj == tot[pozicija+2] and pop[0] == "f":
                                                                                for kw, vw in f2A.items():
                                                                                    if kw.split("__")[0] == tot[pozicija+2]:
                                                                                        mest = rc.seq.find(vw[:15])
                                                                                        
                                                                                        A2_pri_seq = rc.seq[:mest] + vw
                                                                                        
                                                                                        #____________________________________
                                                                                        if len(A2_pri_seq) >=61:
                                                                                            A2_pri_seq= A2_pri_seq[-60:]
                                                                                        #____________________________________
                                                                                        
                                                                                        print("______________________2A_primerji__________________________")

                                                                                        print(klj + "__f")
                                                                                        
                                                                                        if "0" in A2_primerji_celi[klj]: 
                                                                                            A2_primerji_celi[klj][0] = A2_pri_seq
                                                                                        else:
                                                                                            A2_primerji_celi[klj] = [A2_pri_seq,"0"]
                                                                                        
                                                                                       
                                                                                        for kj, vlu in fseq_list.items():
                                                                                            for kk, vv in primerji_2A_f.items():
                                                                                                if kj == klj and kk == klj:
                                                                                                    rpr = vlu[0][vlu[0].find(vlu[2].reverse_complement()):]
                                                                                                    print(rpr)
                                                                                                    fTM = vv[2]
                                                                                                    rTM = vv[3]

                                                                                        
                                                                                        uporabne_tm_info(A2_pri_seq,rpr)
                                                                                        celotni_csv_podatki = [klj + "__f" +"2A_primer in pripadajoči reverse primer",str(A2_pri_seq), str(rpr), fTM, rTM]
                                                                                        celotni_csv_podatki.extend(csv_primer_data)
                                                                                        cel_csv_input.append(celotni_csv_podatki)
                                                                                        print("_________Tm__________")
                                                                                        print(fTM)
                                                                                        print(rTM)
                                                                                        
                                                                                        
                                                                                        A2_sekvenca = SeqRecord(
                                                                                            A2_pri_seq,
                                                                                            id = rc.id,
                                                                                            name= rc.id + "_2A_primer",
                                                                                            description= "2A primer za drugo PCR reakcijo"  ,
                                                                                            annotations={"molecule_type": "DNA"}
                                                                                            )
                                                                                        
                                                                                        f_lok_A2 = value.find(A2_pri_seq)
                                                                                         
                                                                                        konstrukt_f2A = SeqFeature(
                                                                                            FeatureLocation(f_lok_A2,f_lok_A2 + len(A2_pri_seq)),
                                                                                            type = "primer_bind",
                                                                                            qualifiers = {"label" : [rc.id + "_f_2"], "note": ["forward primer druge PCR reakcije"], "Tm": [kw.split("__")[1]]}
                                                                                            )
                                                                                        
                                                                                        konstrukt.features.append(konstrukt_f2A)
                                                                                        
                                                                                        A2_primer = SeqFeature(
                                                                                            FeatureLocation(mest, len(A2_pri_seq)),
                                                                                            type = "primer_bind",
                                                                                            qualifiers = {"label" : [rc.id + "_f_2"], "note": ["forward primer druge PCR reakcije"], "Tm": [kw.split("__")[1]]})
                                                                                        
                                                                                        

                                                                        A2_sekvenca.features = []
                                                                        A2_sekvenca.features.append(A2_primer)
                                                                        
                                                                        
                                                               
                                                                
                                                                SeqIO.write(A2_sekvenca,wr,"genbank")
                                                                fg = False
                                                                break

                                        ja = False
                                        inde = zaporedje.index(vrstica)
                                        for g in koncni_primerji_z_2A_podaljski:
                                            if vrstica in g:
                                                ja = True
                                        
                                        if "A2" in zaporedje[inde - 2] and ja is True:
                                            
                                            for kz in fr2A.keys():
                                                if kz.split("__")[0] == tut:
                                                    tm = kz.split("__")[1]
                                                    break
                                            
                                            drugi_pomnozek = A2_pri_seq + sekvenca_po_prvem_pomnozevanju[len(vw):]
                                            with open(os.path.join(Rezultati, "drugi_pomnozek.fasta"),"a") as wr90:
                                                wr90.write(">" + tut + "\n")
                                                wr90.write(str(drugi_pomnozek) + "\n")
                                                
                                            drugo_pom_dat = [tut,"drugo_pomnozevanje", A2_pri_seq, kw.split("__")[1], r_data_pri, tm, drugi_pomnozek, sekvenca_po_prvem_pomnozevanju]
                                            data.append(drugo_pom_dat)
                                        
                                        nin.writerows(data)
                                        SeqIO.write(sekvenca_prvo_pomnozevanje, wr, "genbank")
#______________________________________________________________________________
        
       
        for i in r_zap:
            if i in rez:
                k = r_zap.index(i)
                seq_poisci = sekvence[k]
                najdi = value.find(seq_poisci)
                
                if i == rez[0]:
                    f_Res_mesto = SeqFeature(
                        FeatureLocation(najdi, najdi + len(seq_poisci)),
                        type = "misc_feature",
                        qualifiers = {"label" : [i  + "__restrikcjisko mesto"], "note": ["5' restrikcijsko mesto za insert v vektor"]})
                if i == rez[1]:
                    r_Res_mesto = SeqFeature(
                        FeatureLocation(najdi, najdi + len(seq_poisci)),
                        type = "misc_feature",
                        qualifiers = {"label" : [i  + "__restrikcjisko mesto"], "note": ["3' restrikcijsko mesto za insert v vektor"]})
    
        
        konstrukt.features.append(f_Res_mesto)
        konstrukt.features.append(r_Res_mesto)
        kok = 0
        sok = 0
        for p in aba [1:-1]:
            if "kozak" in p:
                with open(os.path.join(kozak, "vsi_kozak.fasta"),"r") as rd1:
                    for rf in SeqIO.parse(rd1,"fasta"):
                        if rf.id == p:
                            poz = value.find(rf.seq)                            
                            globals()[rf.id] = SeqFeature(
                                FeatureLocation(poz, poz+len(rf.seq)),
                                type = "misc_feature",
                                qualifiers = {"label" : [rf.id], "note": ["kozak zaporedje"]})
                            konstrukt.features.append(globals()[rf.id])
                            if sok == 0:
                                sok = sok + 1
                                START = rf.seq.find("ATG") + poz
                                
                                
            elif "2A" in p:
                with open(os.path.join(dvaA,"A2.fasta"),"r") as rd2:
                    for rc in SeqIO.parse(rd2,"fasta"):
                        if rc.id == p:
                            poz = value.find(rc.seq)
                            globals()[rc.id] = SeqFeature(
                                FeatureLocation(poz, poz+len(rc.seq)),
                                type = "misc_feature",
                                qualifiers = {"label" : [rc.id], "note": ["2A zaporedje"]})
                            konstrukt.features.append(globals()[rc.id])
            
            else:
                with open(os.path.join(primer_input,"primer3inpt.fasta"),"r") as rd3:
                    kok = kok + 1
                    for re in SeqIO.parse(rd3,"fasta"):
                        if not re.seq == "":
                            if re.id.split("|")[1] == p:
                                poz = value.find(re.seq[0:30])
                                if zapzap[-1] in re.id.split("|")[0]:
                                    STOP = poz + len(re.seq)
                                globals()[re.id] = SeqFeature(
                                    FeatureLocation(poz, poz+len(re.seq)),
                                    type = "CDS",
                                    qualifiers = {"label" : [re.id]}) #"note": [zapzap[kok-1]]
                                konstrukt.features.append(globals()[re.id])
                                
        CDS_cel = SeqFeature(
            FeatureLocation(START, STOP),
            type = "CDS",
            qualifiers = {"label" : ["nelocen_CDS"], "note": ["celoten nelocen CDS"]})
        konstrukt.features.append(CDS_cel)
        
        with open(os.path.join(Rezultati, "asembly_koncni.fasta"),"w") as as_k:
            as_k.write(">" + key +"\n")
            as_k.write(str(value) + "\n")
            
#IPL_primer2 design____________________________________________________________

        try:
            print("________________IPL__________________________")
            with open(os.path.join(Rezultati,"primerji.fasta"),"a") as primer:

                IPL_p_2A_seq = A2_dic["ERBV_2A"][0] + kozak_dic["K503kozak_XuLetal"][0] + IPL_sekvence["X82644yopt"][3:-3] + A2_dic["Opbc18_2A"][0]
                print(IPL_p_2A_seq)
                forward_start1 = IPL_p_2A_seq.find(ne2a_dic["X82644yopt"][0])
                print(forward_start1)
                reverse_end1 = IPL_p_2A_seq.find(ne2a_dic["X82644yopt"][1].reverse_complement()) + len(ne2a_dic["X82644yopt"][1])
                print(reverse_end1)
                IPL1_pri_seq = IPL_p_2A_seq[forward_start1:reverse_end1]
                aaaa = IPL1_pri_seq.find(IPL_sekvence["X82644yopt"][3:-3])
                
                print(aaaa)
                pozicija_IPL = IPL1_pri_seq.find(IPL_sekvence["X82644yopt"][3:-3])

                start_left = 0
                len_left = pozicija_IPL
                start_right = pozicija_IPL + len(IPL_sekvence["X82644yopt"][3:-3])
                len_right = len(IPL1_pri_seq) - start_right

                ok_region_str = str(start_left) + "," + str(len_left) + "," + str(start_right) +"," + str(len_right) + ";"
                EX_START = IPL1_pri_seq.find(IPL_sekvence["X82644yopt"][3:-3])
                ex_end = len(IPL_sekvence["X82644yopt"][3:-3])

                print("IPL_PRI_SEQ")
                print(IPL1_pri_seq)
                
                seq_args = {
                    'SEQUENCE_ID': "r__A2__" + "X82644yopt",
                    'SEQUENCE_TEMPLATE': IPL1_pri_seq,
                    'PRIMER_TASK':  'pick_pcr_primers',
                    #'SEQUENCE_INCLUDED_REGION': [0,len(IPL1_pri_seq)],
                    'SEQUENCE_EXCLUDED_REGION': [str(EX_START) + "," + str(ex_end)],
                    
                    'SEQUENCE_FORCE_LEFT_START': [0]
                    #'SEQUENCE_FORCE_RIGHT_START': [len(IPL1_pri_seq)]
                }
                
                primer_params = {
                #'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 0,
                'PRIMER_OPT_SIZE': 17,
                'PRIMER_MIN_SIZE': 12,
                'PRIMER_MAX_SIZE': max_size,
                'PRIMER_MIN_TM': 50,
                'PRIMER_OPT_TM': 60,
                'PRIMER_MAX_TM': 72,
                'PRIMER_NUM_RETURN': 10,
                'PRIMER_MIN_GC': 0,
                'PRIMER_MAX_GC': 100,     
                'PRIMER_PRODUCT_SIZE_RANGE':[300, len(IPL1_pri_seq)],          
                'PRIMER_PAIR_MAX_DIFF_TM' : 5.0,
                'PRIMER_EXPLAIN_FLAG': 1,
                
                'PRIMER_MAX_HAIRPIN_TH ': 1000000,
                'PRIMER_MAX_SELF_ANY_TH':   10000000,
                'PRIMER_MAX_SELF_END_TH': 10000000,
                'PRIMER_MAX_HAIRPIN_TH': 10000000,
                'PRIMER_DNA_CONC': 500.0,
                'PRIMER_DNTP_CONC': 0.8
                
                }
                                    
                result = bindings.design_primers(primer_params,seq_args)
                
                try:
                    print(result["PRIMER_LEFT_0_SEQUENCE"])
                    print(result["PRIMER_RIGHT_0_SEQUENCE"])
                     
                except:
                    print(".................................")
                    print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
                    print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
                    print(".................................")
                    continue   
                
                
                for i in range(st_pr):
                    try:
                        IPL_pri_forw = Seq(result[f'PRIMER_LEFT_{i}_SEQUENCE']) 
                        IPL_pri_rev = Seq(result[f'PRIMER_RIGHT_{i}_SEQUENCE'])
                        fTm = str(result[f'PRIMER_LEFT_{i}_TM'])
                        rTm = str(result[f'PRIMER_RIGHT_{i}_TM'])
                        
                        smak = key.strip()
                        
                        primer.write(f">f_{i}_"+  "X82644yopt_IPL"  +"__"+ fTm + "\n")
                        primer.write(str(IPL_pri_forw)+ "\n" )
                        primer.write(f">r_{i}_"+  "X82644yopt_IPL"  +"__"+ rTm + "\n")
                        primer.write(str(IPL_pri_rev)+ "\n" )
        
                        #primerji_2A_f[key] = [a,b,fTm,rTm]
                
                    except:
                        continue
    
                find_fow = IPL_p_2A_seq.find(IPL_pri_forw)
                find_rev =IPL_p_2A_seq.find(IPL_pri_rev.reverse_complement())
                nov_IPL_f = IPL_p_2A_seq[find_fow-13:find_fow+len(IPL_pri_forw)]
                nov_IPL_r =IPL_p_2A_seq[find_rev:find_rev+ len(IPL_pri_rev)+15].reverse_complement()
                
                
                print("________________IPL_Polni_primerji__________________")
                print(nov_IPL_f)
                print(nov_IPL_r)
                
                lokacija_f_na_celotni = value.find(nov_IPL_f)
                lokacija_r_na_celotni = value.find(nov_IPL_r.reverse_complement())
                F_IPL =SeqFeature(
                    FeatureLocation(lokacija_f_na_celotni,lokacija_f_na_celotni+len(nov_IPL_f)),
                    type = "primer_bind",
                    qualifiers = {"label" : ["forward_primer_IPL2"]}
                    )
                konstrukt.features.append(F_IPL)
                
                
                R_IPL = SeqFeature(
                    FeatureLocation(lokacija_r_na_celotni,lokacija_r_na_celotni +len(nov_IPL_r),strand = -1),
                    type = "primer_bind",
                    qualifiers = {"label" : ["reverse_primer_IPL2"]}
                    )
                konstrukt.features.append(R_IPL)
            

            uporabne_tm_info(nov_IPL_f,nov_IPL_r)
            
            celotni_csv_podatki = ["novi_IPL_primerji",str(nov_IPL_f) , str(nov_IPL_r) , fTm , rTm]
            celotni_csv_podatki.extend(csv_primer_data)
            cel_csv_input.append(celotni_csv_podatki)

        except:
            print("ni IPL")
        SeqIO.write(konstrukt, wr, "genbank")
        
        #xhoI obdži le G (Gibson razgrajuje 3->5, razgradi štrleči konec) 5 ((ctcga)(g)) 3
        #XbaI obdrži samo T 5 ((t)(ctaga)) 3
        #inf_farme podaljški nemudoma nad rezalnima mestoma

with open(os.path.join(Rezultati,"primerji2.csv"),"a", newline="") as neki:
    nin = csv.writer(neki)
    nin.writerows(cel_csv_input)

with open (os.path.join(Skrajsan_GB,key + ".gb"),"a") as wr2:
    SeqIO.write(konstrukt,wr2,"genbank")     


