
'Tagging program'


import os
from Bio import SeqIO
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

#obdelava_introni
obd_int = os.path.join(wdir,"obdelava_intronskih")

obd_seq = os.path.join(obd_int, "z_obdel_seq")


#primerji
Primerji = os.path.join(wdir, "Primerji")

pr_out = os.path.join(Primerji, "primer3output")

#______________________________________________________________________________
if os.path.exists(os.path.join(Rezultati,"Full_seq.gb")):
    os.remove(os.path.join(Rezultati,"Full_seq.gb"))
if os.path.exists(os.path.join(Rezultati,"OL661280__His.gb")):
    os.remove(os.path.join(Rezultati,"OL661280__His.gb"))

cel_csv_input = [[],["ICS_IPL_SAMT_HIS_podaljški","fPrimer","rPrimer","fTm","rTm","fHp","dh","ds","rHp","dh","ds","f_Dimer","dh","ds","r_Dimer","dh","ds","HeteroDimer","dh","ds"]]


f_IPL_ful = Seq("TGAATTGAATCCAGGCCCACCAACCATGT")
r_IPL_ful = Seq("CACGTCCCCGCATAACTTCAGCAAATCATAATTTGATTT")
f_IPL =("GGCCCACCAACCATG")
r_IPL = Seq("CTTCAGCAAATCATAATTTGATTT")

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
                
                print()
                print("forward_hairpin")
                print(fHp)
                print("reverse_hairpin")

                        
                print(rHp)
                print()
                try:
                    for line in rHp.ascii_structure_lines:
                            art = line.split("\t")[1]
                            print(art)
                except:
                    print("ni_hp")
                    
                print()
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


Tagi = {"His":Seq("CATCACCATCACCATCAC"),}
        #"HA":Seq("TATCCATACGATGTTCCAGATTACGCT"),
        #"FLAGx3":Seq("GATTATAAGGATCACGACGGTGACTACAAGGATCATGACGGTGATTATAAGGATCATGATGGT")} #Vsi so optimizirani za kvasovko

STOP_kodoni = {Seq("TAG"), 
               Seq("TAA"), 
               Seq("TGA")}

forward_podaljšek = Seq("CAAGGAGAAAAAACCCCGGAT")
reverse_podaljšek = Seq("TCATGTAATTAGTTATGTCACGCTTACA").reverse_complement()


Dodatek_STOP_kodon = Seq("TAA")

vrsta_tag = "SAMT"

uporabljeni_tagi = ("His")

rezalno_mesto_koncna = "XhoI"

dodan_encim = getattr(Restriction, rezalno_mesto_koncna, None)
dodan_encim_seq = dodan_encim.site

ime_vrste = {}

primerji = {}
with open (os.path.join(obd_seq, "zdr_u_seq.fasta"),"r") as r:
    for rz in SeqIO.parse(r,"fasta"):
        neki = rz.id.split("|")
        ime = neki[1]
        vrsta = neki[0].split("_")[1]
        try:
            ime_vrste[vrsta].append(ime)
        except:
            ime_vrste[vrsta] = [ime]
zaporedje = []
     
with open(os.path.join(zapi),"r") as zapo:
    for line in zapo:
        if not line.startswith(("A2","kozak","rez")):
            zaporedje.append((line.split("__")[0]).strip())

print("zaporedje")
print(zaporedje)

with open(os.path.join(Rezultati, "prvi_pomnozek.fasta"),"r") as prva_dat, open (os.path.join(Rezultati,"asembly.fasta"),"r") as assembly, open(os.path.join(Rezultati,"A2_primerji.fasta"),"r") as primer1, open(os.path.join(Rezultati,"tag_primerji1.fasta"),"a") as tag_primer_1:
    for record in SeqIO.parse(prva_dat,"fasta"):
        n = 0
        for ki, val in ime_vrste.items():
            if record.id in val:
                z = ki
                break
        
        
        if vrsta_tag == z:
            Sekvenca_taganega = record.seq
                
            for rec in SeqIO.parse(assembly,"fasta"):
                
                if record.id in rec.id:
                    
                    prva_zadnja = rec.id.split("__")
                    
                    zadnja_res = prva_zadnja[-1]
                    encim = getattr(Restriction, zadnja_res, None)
                    
                    Sekvenca_taganega = Sekvenca_taganega[:-len(encim)]
                    
                    while Sekvenca_taganega[-3:] in STOP_kodoni:
                        Sekvenca_taganega = Sekvenca_taganega[:-3]
                        n= n + 1
                        
                    
                    Cela_sekvenca = record.seq[:-(n*3+len(encim.site)+4)]                
                    for tag, seq in Tagi.items():
                        if tag in uporabljeni_tagi:

                            globals()["Cela_sekvenca_" + tag] = Cela_sekvenca + seq + Dodatek_STOP_kodon + dodan_encim_seq + Seq("TATA").reverse_complement()
                            
                            Taggane_sekvence = SeqRecord(
                                globals()["Cela_sekvenca_" + tag],
                                id = record.id + "__" + tag,
                                name= record.id +"__" + tag + "__produkt",
                                description="produkt pomnozevanja " + record.id + "s tag-animi primerji"  ,
                                annotations={"molecule_type": "DNA"}
                                )
                            
                            for re in SeqIO.parse(primer1,"fasta"):
                                ide = re.id.split("_")[2]
                                for k, v in ime_vrste.items():
                                    if ide in v and k in vrsta_tag:
                                        if re.id.startswith("f"):
                                            forward_primer = re.seq
                                            začetek_f = globals()["Cela_sekvenca_" + tag].find(forward_primer)
                                            f_primer = SeqFeature(
                                                FeatureLocation( začetek_f,začetek_f + len(forward_primer)),
                                                type = "primer_bind",
                                                qualifiers = {"label" : ["f_primer__" + ide], "note": ["forward_primer " + ide], "Tm" : [re.id.split("__")[-1]]}
                                                )
                                            seq_args = {
                                                
                                                'SEQUENCE_ID': tag,
                                                'SEQUENCE_PRIMER': forward_primer,
                                                'SEQUENCE_TEMPLATE': Sekvenca_taganega,
                                            }
                                            
                                            primer_params = {        
                                            'PRIMER_OPT_SIZE': 18,
                                            'PRIMER_MIN_SIZE': 15,
                                            'PRIMER_MIN_TM': 50,
                                            'PRIMER_OPT_TM': 60,
                                            'PRIMER_MAX_TM': 72,
                                            'PRIMER_MIN_GC': 0,
                                            'PRIMER_MAX_GC': 100,                        
                                            'PRIMER_PRODUCT_SIZE_RANGE': [len(Sekvenca_taganega),len(Sekvenca_taganega)],
                                            'PRIMER_EXPLAIN_FLAG': 1,
                                            
                                            'PRIMER_PAIR_MAX_DIFF_TM' : 5,
                                            'PRIMER_DNA_CONC': 500.0,
                                            'PRIMER_DNTP_CONC': 0.8
                                            }
                                                                
                                            result = bindings.design_primers(primer_params,seq_args)
                                            
                                            
                                            try:
                                                print(result["PRIMER_LEFT_0_SEQUENCE"])
                                                print(result["PRIMER_RIGHT_0_SEQUENCE"])

                                                primerji[ide] = [result["PRIMER_RIGHT_0_TM"],result["PRIMER_RIGHT_0_SEQUENCE"]]

                                            except:
                                                print(".................................")
                                                print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
                                                print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
                                                print(".................................")
                                                continue
                                            
                                            print(primerji)
                                            
                                            r_pri = Seq(result["PRIMER_RIGHT_0_SEQUENCE"]).reverse_complement()                                            
                                            začetek_r = Sekvenca_taganega.find(r_pri)
                                            
                                            r_podalj = globals()["Cela_sekvenca_" + tag][začetek_r + len(r_pri):]
                                            
                                            
                                            r_primer = SeqFeature(
                                                FeatureLocation(začetek_r, začetek_r + len(result["PRIMER_RIGHT_0_SEQUENCE"]), strand=-1),
                                                type = "primer_bind",
                                                qualifiers = {"label" : ["r_primer__" + ide], "note": ["reverse_primer " + ide], "Tm" : [result["PRIMER_RIGHT_0_TM"]]}
                                                )
                                            
                                            
                                            Taggane_sekvence.features.append(f_primer)
                                            Taggane_sekvence.features.append(r_primer)
                                            
                                            
                                            with open(os.path.join(Rezultati, ide + "__" + tag + ".gb"),"a") as gb:
                                                SeqIO.write(Taggane_sekvence, gb, "genbank")
                                            
                                            for ima, pri_sekv in primerji.items():
                                                tag_primer_1.write(">" + ima + "__"+ str(pri_sekv[0]) + "\n")
                                                tag_primer_1.write(str(pri_sekv[1]) + "\n")
                                                                                                
                                                polna_seq = Seq(dodan_encim_seq).reverse_complement()+ Dodatek_STOP_kodon.reverse_complement() + seq.reverse_complement() + pri_sekv[1]
                                                tag_primer_1.write(str(polna_seq) + "\n")




primerji_in_matrice = {}


print()
print("__________GIBSON_PRIMERJI_________________")

print(zaporedje)
flag_r = False
flag_f = False
with open(os.path.join(Rezultati, "prvi_pomnozek.fasta"),"r") as a:
    for record in SeqIO.parse(a,"fasta"):
        print(record.id)
        for key,value in ime_vrste.items():

            if record.id in value and key in zaporedje[0]:
                sekvenca_za_primanje = record.seq[2:]
                with open(os.path.join(Rezultati,"A2_primerji.fasta"),"r") as primer1:
                    for rec in SeqIO.parse(primer1, "fasta"):
                        if record.id in rec.id.split("_") and rec.id.startswith("r"):
                            primer = rec.seq
                            primerji_in_matrice[record.id] = [sekvenca_za_primanje,primer]
                            break
                        elif rec.id.startswith("r"):
                            primer = r_IPL
                            flag_r = True

                seq_args = {
                    
                    'SEQUENCE_ID': record.id + "__GIBSON__",
                    'SEQUENCE_PRIMER_REVCOMP': primer,
                    'SEQUENCE_TEMPLATE': sekvenca_za_primanje,
                    'SEQUENCE_FORCE_LEFT_START': 0,
                }
                
                primer_params = {        
                'PRIMER_OPT_SIZE': 18,
                'PRIMER_MIN_SIZE': 15,
                'PRIMER_MIN_TM': 50,
                'PRIMER_OPT_TM': 60,
                'PRIMER_MAX_TM': 72,
                'PRIMER_MIN_GC': 0,
                'PRIMER_MAX_GC': 100,                        
                'PRIMER_PRODUCT_SIZE_RANGE': [len(sekvenca_za_primanje)-100,len(sekvenca_za_primanje)],
                'PRIMER_EXPLAIN_FLAG': 1,
                'PRIMER_MAX_HAIRPIN_TH': 1000000,
                'PRIMER_MAX_SELF_ANY_TH':   10000000,
                'PRIMER_MAX_SELF_END_TH': 10000000,
                'PRIMER_PAIR_MAX_DIFF_TM' : 5,
                'PRIMER_DNA_CONC': 500.0,
                'PRIMER_DNTP_CONC': 0.8
                }
                                    
                result = bindings.design_primers(primer_params,seq_args)
                
                
                try:
                    print(result["PRIMER_LEFT_0_SEQUENCE"])
                    print(result["PRIMER_RIGHT_0_SEQUENCE"])

                    f_gibson_primerji_Tm = [result["PRIMER_LEFT_0_TM"],result["PRIMER_RIGHT_0_TM"]]

                except:
                    print(".................................")
                    print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
                    print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
                    print(".................................")
                    continue
                
                print(primerji)
                        
                for_pr_gibson = forward_podaljšek + Seq(result["PRIMER_LEFT_0_SEQUENCE"])
                
                with open(os.path.join(Rezultati,"primerji.csv"),"r", newline="") as neki2: 
                    reader = csv.reader(neki2, delimiter=',')
                    for row in reader:
                        if "drugo_pomnozevanje" in row and record.id in row and flag_r is False:
                            forr_pr_gibson = row[4]
                        elif flag_r is True:
                            forr_pr_gibson = r_IPL_ful
                        
                        
            elif record.id in value and key in zaporedje[-1]:
                sekvenca_za_primanje = globals()["Cela_sekvenca_" + tag]
                with open(os.path.join(Rezultati,"A2_primerji.fasta"),"r") as primer1:
                    for rec in SeqIO.parse(primer1, "fasta"):
                        if record.id in rec.id.split("_") and rec.id.startswith("f"):
                            primer = rec.seq
                            primerji_in_matrice[record.id] = [sekvenca_za_primanje,primer]
                            break
                        elif rec.id.startswith("f"):
                            primer = f_IPL
                            flag_f = True
                
                seq_args = {
                    
                    'SEQUENCE_ID': record.id + "__GIBSON__",
                    'SEQUENCE_PRIMER': primer,
                    'SEQUENCE_TEMPLATE': sekvenca_za_primanje,
                    'SEQUENCE_FORCE_RIGHT_START': len(sekvenca_za_primanje)-1,
                }
                
                primer_params = {        
                'PRIMER_OPT_SIZE': 18,
                'PRIMER_MIN_SIZE': 15,
                'PRIMER_MIN_TM': 50,
                'PRIMER_OPT_TM': 60,
                'PRIMER_MAX_TM': 72,
                'PRIMER_MIN_GC': 0,
                'PRIMER_MAX_GC': 100,                        
                'PRIMER_PRODUCT_SIZE_RANGE': [len(sekvenca_za_primanje)-100,len(sekvenca_za_primanje)],
                'PRIMER_EXPLAIN_FLAG': 1,
                'PRIMER_MAX_HAIRPIN_TH': 1000000,
                'PRIMER_MAX_SELF_ANY_TH':   10000000,
                'PRIMER_MAX_SELF_END_TH': 10000000,
                
                'PRIMER_PAIR_MAX_DIFF_TM' : 5,
                'PRIMER_DNA_CONC': 500.0,
                'PRIMER_DNTP_CONC': 0.8
                }
                                    
                result = bindings.design_primers(primer_params,seq_args)
                
                
                try:
                    print(result["PRIMER_LEFT_0_SEQUENCE"])
                    print(result["PRIMER_RIGHT_0_SEQUENCE"])

                    r_gibson_primerji_Tm = [result["PRIMER_LEFT_0_TM"],result["PRIMER_RIGHT_0_TM"]]

                except:
                    print(".................................")
                    print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
                    print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
                    print(".................................")
                    continue
                
                print(primerji)
                with open(os.path.join(Rezultati,"primerji.csv"),"r", newline="") as neki2: 
                    reader = csv.reader(neki2, delimiter=',')
                    for row in reader:
                        if "drugo_pomnozevanje" in row and record.id in row and flag_f is False:
                            rewf_pr_gibson = row[2]
                        elif flag_f is True:
                            rewf_pr_gibson = f_IPL_ful

                
                rew_pr_gibson =  reverse_podaljšek + Seq(result["PRIMER_RIGHT_0_SEQUENCE"])
                print(rew_pr_gibson)
                
    
    
print()


        

try:
    print("forward")
    print(for_pr_gibson)
    print(forr_pr_gibson)
    uporabne_tm_info (for_pr_gibson,forr_pr_gibson)
    celotni_csv_podatki = ["Gibson_podaljški_forward",for_pr_gibson,forr_pr_gibson]
    celotni_csv_podatki.extend(f_gibson_primerji_Tm)
    celotni_csv_podatki.extend(csv_primer_data)
    cel_csv_input.append(celotni_csv_podatki)
except:
    print("ni forward")

try:   
    print("reverse")
    print(rew_pr_gibson)
    print(rewf_pr_gibson)
    uporabne_tm_info (rewf_pr_gibson,rew_pr_gibson)
    celotni_csv_podatki = ["Gibson_podaljški_reverse",rewf_pr_gibson,rew_pr_gibson]
    celotni_csv_podatki.extend(r_gibson_primerji_Tm)
    celotni_csv_podatki.extend(csv_primer_data)
    cel_csv_input.append(celotni_csv_podatki)
except:
    print("ni reverse")


primerj ={}


with open(os.path.join(Rezultati,"asembly_koncni.fasta"),"r") as as_k:                                              
    for rec in SeqIO.parse(as_k,"fasta"):
        XhoIrez =rec.seq.find(dodan_encim_seq)

        for i in STOP_kodoni:
            if i in rec.seq[XhoIrez-3:XhoIrez]:
                a = rec.seq[XhoIrez-3:XhoIrez].find(i)
        
                seq = rec.seq[:XhoIrez-(a+3)]+ Seq("CATCACCATCACCATCAC")+Seq("TAACTCGAGTATA")  + Seq(reverse_podaljšek).reverse_complement()
                
                Full_seq = SeqRecord(
                    seq,
                    id = record.id + "__" + tag,
                    name= "Full_seq",
                    description="Polna_seq",
                    annotations={"molecule_type": "DNA"}
                    )
                His_TAG = seq.find(Seq("CATCACCATCACCATCAC"))
                
                His_ano = SeqFeature(
                    FeatureLocation(His_TAG,His_TAG +len("CATCACCATCACCATCAC")),
                    type = "CDS",
                    qualifiers = {"label" : ["His_Tag"]}
                    )
                #seqfind = seq.find(Seq("CTCGAG"))+2
                print("________FUL_PRIMERJI___________")
                seq_args = {
                    
                    'SEQUENCE_ID': "Ful_primerji",
                    'SEQUENCE_TEMPLATE': seq,
                    'SEQUENCE_FORCE_RIGHT_START': [seq.find(reverse_podaljšek.reverse_complement())-1],
                    'SEQUENCE_FORCE_LEFT_START': [len(forward_podaljšek)]
                }
                
                primer_params = {        
                'PRIMER_OPT_SIZE': 18,
                'PRIMER_MIN_SIZE': 15,
                'PRIMER_MIN_TM': 50,
                'PRIMER_OPT_TM': 60,
                'PRIMER_MAX_TM': 72,
                'PRIMER_MIN_GC': 0,
                'PRIMER_MAX_GC': 100,                        
                'PRIMER_EXPLAIN_FLAG': 1,
                'PRIMER_PRODUCT_SIZE_RANGE': [1800,3000],
                'PRIMER_PAIR_MAX_DIFF_TM' : 5,
                
                #'PRIMER_MAX_HAIRPIN_TH ': 1000000,
                #'PRIMER_MAX_SELF_ANY_TH':   10000000,
                #'PRIMER_MAX_SELF_END_TH': 10000000,
                #'PRIMER_MAX_HAIRPIN_TH': 10000000,
                
                'PRIMER_DNA_CONC': 500.0,
                'PRIMER_DNTP_CONC': 0.8
                }
                                    
                result = bindings.design_primers(primer_params,seq_args)
                
                
                try:
                    print(result["PRIMER_LEFT_0_SEQUENCE"])
                    print(result["PRIMER_RIGHT_0_SEQUENCE"])
                    rPr = result["PRIMER_RIGHT_0_SEQUENCE"]
                    fPr = result["PRIMER_LEFT_0_SEQUENCE"]
                
                    primerj[ide] = [result["PRIMER_LEFT_0_TM"],fPr,result["PRIMER_RIGHT_0_TM"],rPr]
                
                except:
                    print(".................................")
                    print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
                    print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
                    print(".................................")
                    
                    #continue
                
                print(primerj)
                
                f_loc = seq.find(fPr)
                r_loc = seq.find(Seq(rPr).reverse_complement())
                
                f_loc_ano = SeqFeature(
                    FeatureLocation(f_loc, f_loc+ len(fPr)),
                    type = "primer_bind",
                    qualifiers = {"label" : ["forward_celoten_primer"]})
                
                r_loc_ano = SeqFeature(
                    FeatureLocation(r_loc,r_loc + len(rPr),strand = -1),
                    type = "primer_bind",
                    qualifiers = {"label" : ["reverse_celoten_primer"]})

                print()
                
                fTm = result["PRIMER_LEFT_0_TM"]
                rTm =result["PRIMER_RIGHT_0_TM"]
                
                print("nepodaljšani")
                print(fPr)
                print(rPr)
                uporabne_tm_info(fPr,rPr)
                print()
                
                celotni_csv_podatki = ["nepodaljšani_za_OE",fPr,rPr, fTm, rTm]
                celotni_csv_podatki.extend(csv_primer_data)
                cel_csv_input.append(celotni_csv_podatki)
                
                with open(os.path.join(Rezultati, "Full_seq" + ".gb"),"a") as gb:
                    Full_seq.features = [His_ano,f_loc_ano,r_loc_ano]#,f_pod_loc,r_pod_loc]
                    SeqIO.write(Full_seq, gb, "genbank")
            
print(cel_csv_input)

with open(os.path.join(Rezultati,"primerji2.csv"),"a", newline="") as neki:
    nin = csv.writer(neki)
    nin.writerows(cel_csv_input)
    
                
                
                                            
                                            
                            
                