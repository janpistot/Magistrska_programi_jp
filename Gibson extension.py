# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 15:16:49 2025

@author: HP
"""
import subprocess, tempfile

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



f_IPL_ful = Seq("TGAATTGAATCCAGGCCCACCAACCATGT")
r_IPL_ful = Seq("CACGTCCCCGCATAACTTCAGCAAATCATAATTTGATTT")
f_IPL =("GGCCCACCAACCATG")
r_IPL = Seq("CTTCAGCAAATCATAATTTGATTT")

cel_csv_input =[]

potrebujejo_2a_primerje_za_GIBSON = ["YFP"]
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
                    
                    

    
print()
print("__________GIBSON_PRIMERJI_________________")

forward_podaljšek = Seq("CAAGGAGAAAAAACCCCGGAT")
reverse_podaljšek = Seq("TCATGTAATTAGTTATGTCACGCTTACA").reverse_complement()

primerji_in_matrice = {}

zaporedje =[]
zappredje = []
with open(os.path.join(zapi),"r") as zap:
    for line in zap:
        zaporedje.append(line.strip())


for z in zaporedje:
    if not z.startswith(("rez","A2","kozak")):
        zappredje.append(z)

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
            
            
with open(os.path.join(Rezultati, "prvi_pomnozek.fasta"),"r") as a:
    for record in SeqIO.parse(a,"fasta"):
        print(record.id)
        for i in potrebujejo_2a_primerje_za_GIBSON:
            if i in record.id:
                for key,value in ime_vrste.items():
        
                    if record.id in value and key in zappredje[0]:
                        sekvenca_za_primanje = record.seq

                     
                        primerji_in_matrice[record.id] = [sekvenca_za_primanje,r_IPL_ful]
                        print("forward")
                        seq_args = {
                            
                            'SEQUENCE_ID': record.id + "__GIBSON__",
                            'SEQUENCE_PRIMER_REVCOMP': r_IPL,
                            'SEQUENCE_TEMPLATE': sekvenca_za_primanje,
                            'SEQUENCE_FORCE_LEFT_START': 0,
                        }
                        
                        primer_params = {        
                        'PRIMER_OPT_SIZE': 15,
                        'PRIMER_MIN_SIZE': 12,
                        'PRIMER_MIN_TM': 50,
                        'PRIMER_OPT_TM': 60,
                        'PRIMER_MAX_TM': 72,
                        'PRIMER_MIN_GC': 0,
                        'PRIMER_MAX_GC': 100,
                        
                        'PRIMER_PAIR_MAX_DIFF_TM' : 5.0,
                        'PRIMER_EXPLAIN_FLAG': 1,
                        
                        'PRIMER_MAX_HAIRPIN_TH': 1000000,
                        'PRIMER_MAX_SELF_ANY_TH':   10000000,
                        'PRIMER_MAX_SELF_END_TH': 10000000,
                        'PRIMER_DNA_CONC': 500.0,
                        'PRIMER_DNTP_CONC': 0.8,
                        
                        
                        
                        'PRIMER_PRODUCT_SIZE_RANGE': [len(sekvenca_za_primanje)-20,len(sekvenca_za_primanje)]
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
                        
                        
                        forr_pr_gibson = r_IPL_ful
                                
                        print("forward")
                        print(for_pr_gibson)
                        print(forr_pr_gibson)
                        uporabne_tm_info (for_pr_gibson,forr_pr_gibson)
                        celotni_csv_podatki = ["Gibson_podaljški_forward",for_pr_gibson,forr_pr_gibson]
                        celotni_csv_podatki.extend(f_gibson_primerji_Tm)
                        celotni_csv_podatki.extend(csv_primer_data)
                        cel_csv_input.append(celotni_csv_podatki)
                    
                    
                    elif record.id in value and key in zappredje[-1]:
                        sekvenca_za_primanje = record.seq
                        primerji_in_matrice[record.id] = [sekvenca_za_primanje,f_IPL_ful]
                        print(sekvenca_za_primanje)
                        print(f_IPL_ful)
                        
                        print("reverse")                   
                        seq_args = {
                            
                            'SEQUENCE_ID': record.id + "__GIBSON__",
                            'SEQUENCE_PRIMER': f_IPL,
                            'SEQUENCE_TEMPLATE': sekvenca_za_primanje,
                            'SEQUENCE_FORCE_RIGHT_START' : len(sekvenca_za_primanje)-1,
                        }
                        
                        primer_params = {        
                        'PRIMER_OPT_SIZE': 15,
                        'PRIMER_MIN_SIZE': 12,
                        'PRIMER_MIN_TM': 50,
                        'PRIMER_OPT_TM': 60,
                        'PRIMER_MAX_TM': 72,
                        'PRIMER_MIN_GC': 0,
                        'PRIMER_MAX_GC': 100,
                        
                        'PRIMER_PAIR_MAX_DIFF_TM' : 5.0,
                        'PRIMER_EXPLAIN_FLAG': 1,
                        
                        
                        'PRIMER_MAX_HAIRPIN_TH ': 1000000,
                        'PRIMER_MAX_SELF_ANY_TH':   10000000,
                        'PRIMER_MAX_SELF_END_TH': 10000000,
                        'PRIMER_MAX_HAIRPIN_TH': 10000000,
                        'PRIMER_DNA_CONC': 500.0,
                        'PRIMER_DNTP_CONC': 0.8,
                        
                        
                        'PRIMER_PRODUCT_SIZE_RANGE': [len(sekvenca_za_primanje)-20,len(sekvenca_za_primanje)]
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
                        
                        rewf_pr_gibson = f_IPL_ful
                                   
                                    
                        rew_pr_gibson =  reverse_podaljšek + Seq(result["PRIMER_RIGHT_0_SEQUENCE"])
                        print("reverse")
                        print(rew_pr_gibson)
                        print(rewf_pr_gibson)
                        uporabne_tm_info (rewf_pr_gibson,rew_pr_gibson)
                        celotni_csv_podatki = ["Gibson_podaljški_reverse",rewf_pr_gibson,rew_pr_gibson]
                        print(celotni_csv_podatki)
                        celotni_csv_podatki.extend(r_gibson_primerji_Tm)
                        celotni_csv_podatki.extend(csv_primer_data)
                        cel_csv_input.append(celotni_csv_podatki)

    
print()


with open(os.path.join(Rezultati,"primerji2.csv"),"a", newline="") as neki:
    nin = csv.writer(neki)
    nin.writerows(cel_csv_input)