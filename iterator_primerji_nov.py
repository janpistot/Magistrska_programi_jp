
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import *

from primer3 import bindings
from primer3 import calc_tm
from primer3 import calc_hairpin
from primer3 import calc_homodimer
from primer3 import calc_heterodimer

import math

import pickle
import os

import pandas as pd
import openpyxl


#__________________________________POPRAVI_ZA_SVOJ_RAČUNALNIK__________________
#pickle file vsebuje podatke iz Urejanje_UI_3.01.py, brez pravilnega vnosa program ne bo deloval

file_path = r"C:\Users\HP\Desktop\Python_proj\UI testing\Test_fasta\MeSA_SK_rez.pkl"

#_v to lokacijo program_napiše_končne datoteke excel. POPRAVI!!!_______________
test_program_lok = r"C:\Users\HP\Desktop\test_program"
#______________________________________________________________________________





with open(file_path,"rb") as file_s_sekvencami:
    sk_ime = pickle.load(file_s_sekvencami)
    
    globals()["skupine_popolni_record"] = pickle.load(file_s_sekvencami)
    globals()["lista_sestavljenih_sekvenc"] = pickle.load(file_s_sekvencami)

#______________________________________________________________________________

koncentracije_primerjev_v_M = 0.5*(1e-6)
zeljen_delez = 0.1
    
zelen_delez_homodimerov = (koncentracije_primerjev_v_M * zeljen_delez)/2
delez_nehomodimeriziranih = koncentracije_primerjev_v_M - koncentracije_primerjev_v_M*(1-zeljen_delez)
   
razmerje_za_homodimere = zelen_delez_homodimerov /(pow(delez_nehomodimeriziranih,2))

hp_zelen_lnKxR = -8.314*math.log((1-zeljen_delez)/(zeljen_delez))
homodimer_zelen_lnKxR = -8.314 * math.log(razmerje_za_homodimere)

#hetero_A = pow(2*koncentracije_primerjev_v_M*zeljen_delez-2*koncentracije_primerjev_v_M,2)
#hetero_B = 4*koncentracije_primerjev_v_M
#hetero_C = -1

A = (1-zeljen_delez)*koncentracije_primerjev_v_M
AB = zeljen_delez * koncentracije_primerjev_v_M

heterodimer_K1 = AB / (A * A)#(-hetero_B+math.sqrt(pow(hetero_B, 2)-4*hetero_A*hetero_C))/(2*hetero_A)
hetero_zelen_lnKxR = math.log(heterodimer_K1) * -8.314

imena_elementov = globals()["skupine_popolni_record"]["zaporedje"].copy()

#______________________________________________________________________________

stevilo_izdelanih_primerjev = 10

#______________________________________________________________________________

def test_heterodimer (f_pr, r_pr, tm):
    primerjalni_hetero_dG = hetero_zelen_lnKxR * tm
    heterodimer_TD = calc_heterodimer(f_pr, r_pr)
    heterodimer_dG = (heterodimer_TD.dh - tm*heterodimer_TD.ds)* 4.184
    
    if heterodimer_dG > primerjalni_hetero_dG:
        return True
    else:
        return False
    
def iterator_za_konstrukcijo_optimalnih_podaljškov(podaljšan_primer, prvotni_primer,tm):

    hp_zelen_dG = hp_zelen_lnKxR * tm
    homodimer_zelen_dG = homodimer_zelen_lnKxR*tm
    
    if len(podaljšan_primer)>=60:
        
        primer_iterator = podaljšan_primer[-59:]
    else:
        primer_iterator = podaljšan_primer
    
    podaljšan_primer_hp_dG = -10e9
    podaljšan_primer_homo_dG = -10e9
    
    check_za_celo = 0
    
    print("primer_1")
    
    while True:
        
        primer_iterator = primer_iterator[check_za_celo:]
        
        podaljšan_primer_hp_TD = calc_hairpin(str(primer_iterator))
        podaljšan_primer_hp_dG = (podaljšan_primer_hp_TD.dh-podaljšan_primer_hp_TD.ds*tm)* 4.184             
        
        podaljšan_primer_homo_TD = calc_homodimer(str(primer_iterator))
        podaljšan_primer_homo_dG = (podaljšan_primer_homo_TD.dh - podaljšan_primer_homo_TD.ds*tm)* 4.184
        
        check_za_celo = 1
        
        print()
        print(primer_iterator)
        #print("--__--")
        print(str(podaljšan_primer_homo_dG) + " : " + str(homodimer_zelen_dG))
        print(str(podaljšan_primer_hp_dG) + " : " + str(hp_zelen_dG))
        print("#####")

        if (len(primer_iterator) < len(prvotni_primer)):
            
            print("__Fail__")
            return(Seq(""))
        
        
        
        if (((podaljšan_primer_hp_dG >= hp_zelen_dG) or podaljšan_primer_hp_dG == 0.0)
            and 
            (podaljšan_primer_homo_dG >= homodimer_zelen_dG) or podaljšan_primer_homo_dG == 0.0):
        
            return primer_iterator


def optimizacija_konstrukcije_oligov(neobdelan_podaljsan_smerni_primer, 
                                     prvotni_smerni_primer,
                                     neobdelan_podaljsan_protismerni_primer, 
                                     prvotni_protismerni_primer,
                                     tm):
    
    podaljsan_smerni_primer = iterator_za_konstrukcijo_optimalnih_podaljškov(neobdelan_podaljsan_smerni_primer,prvotni_smerni_primer,tm)
    print("F")
    print(podaljsan_smerni_primer)
    
    if not podaljsan_smerni_primer == Seq(""):
        
        podaljsan_protismerni_primer = iterator_za_konstrukcijo_optimalnih_podaljškov(neobdelan_podaljsan_protismerni_primer,prvotni_protismerni_primer,tm)
        
        print("R")
        print(podaljsan_protismerni_primer)
        if not podaljsan_protismerni_primer == Seq(""):
            
            
            
            podaljsan_smerni_primer_hetero_test = podaljsan_smerni_primer
            podaljsan_protismerni_primer_hetero_test = podaljsan_protismerni_primer
            
            preklop = True           
            iterator_po_krajšanju_oligov = 0
            
            if (test_heterodimer (str(podaljsan_smerni_primer_hetero_test), 
                                    str(podaljsan_protismerni_primer_hetero_test), 
                                    tm)) is False:
    
                while (test_heterodimer (str(podaljsan_smerni_primer_hetero_test), 
                                        str(podaljsan_protismerni_primer_hetero_test), 
                                        tm) is False) :

                    test_za_premajhno_dozino_oligov =(len(podaljsan_smerni_primer_hetero_test) < len(prvotni_smerni_primer) 
                                                      and
                                                      (len(podaljsan_protismerni_primer_hetero_test) < len(prvotni_protismerni_primer)))
                    
                    
                    print("test_primer2")
                    print(podaljsan_smerni_primer_hetero_test)
                    print(podaljsan_protismerni_primer_hetero_test)
                    if test_za_premajhno_dozino_oligov:
                        print("test_fail")
                        return None
                    
                    
                    
                    if preklop:
                        if len(podaljsan_smerni_primer_hetero_test) > len(prvotni_smerni_primer): 
                            podaljsan_smerni_primer_hetero_test = podaljsan_smerni_primer[iterator_po_krajšanju_oligov:]
                            
                        preklop = False
                    
                    else:
                        if len(podaljsan_protismerni_primer_hetero_test) > len(prvotni_protismerni_primer):
                            podaljsan_protismerni_primer_hetero_test = podaljsan_protismerni_primer[iterator_po_krajšanju_oligov:]
    
                        preklop = True
    
                    iterator_po_krajšanju_oligov = iterator_po_krajšanju_oligov + 1
                
                
                if not preklop:
                    while not ((test_heterodimer(str(podaljsan_smerni_primer_hetero_test), 
                                             str(podaljsan_protismerni_primer_hetero_test), 
                                             tm))) or not (len(podaljsan_smerni_primer_hetero_test) == len(podaljsan_smerni_primer)):
                        
                        
                        podaljsan_smerni_primer_hetero_test = podaljsan_smerni_primer[iterator_po_krajšanju_oligov:]
                        iterator_po_krajšanju_oligov = iterator_po_krajšanju_oligov - 1

                elif preklop:
                    while not(test_heterodimer(str(podaljsan_smerni_primer_hetero_test), 
                                             str(podaljsan_protismerni_primer_hetero_test), 
                                             tm)) or (len(podaljsan_protismerni_primer_hetero_test) == len(podaljsan_protismerni_primer)):

                        podaljsan_protismerni_primer_hetero_test = podaljsan_protismerni_primer[iterator_po_krajšanju_oligov:]
                        iterator_po_krajšanju_oligov = iterator_po_krajšanju_oligov - 1

            podaljsani_oligi = [podaljsan_smerni_primer_hetero_test,podaljsan_protismerni_primer_hetero_test]
            
            return(podaljsani_oligi)


def konstrukcija_primerjev(deli, konstrukt, gradnik, izvor, inicializacija,potreben):
    if potreben:
        input_sekvenca = deli[1]
        dolžina_dela = len(input_sekvenca)
        seq_args = {
            'SEQUENCE_TEMPLATE': str(input_sekvenca),
            'SEQUENCE_FORCE_RIGHT_START': [dolžina_dela-1],
            'SEQUENCE_FORCE_LEFT_START': [0],
        }
        primer_params = {        
        'PRIMER_OPT_SIZE': 18,
        'PRIMER_MIN_SIZE': 15,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_MAX_DIFF_TM': 5.0,
        'PRIMER_MAX_HAIRPIN_TH':100,
        'PRIMER_MIN_TM': 54.5,
        'PRIMER_OPT_TM': 60,
        'PRIMER_MAX_TM': 75,
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_MIN_GC': 0,
        'PRIMER_MAX_GC': 100,                        
        'PRIMER_PRODUCT_SIZE_RANGE': [40,dolžina_dela],
        'PRIMER_EXPLAIN_FLAG': 1,
        'PRIMER_DNA_CONC': 500.0,
        'PRIMER_DNTP_CONC': 0.8
        }                   
        result = bindings.design_primers(primer_params,seq_args)
        
        for i in range(10):
            
            try:
                lprimer_str = "PRIMER_LEFT_"+str(i)+"_SEQUENCE"
                rprimer_str = "PRIMER_RIGHT_"+str(i)+"_SEQUENCE"
                
                smerni_primer_tm = float(result["PRIMER_LEFT_"+str(i)+"_TM"]) + 273
                protismerni_primer_tm = float(result["PRIMER_RIGHT_"+str(i)+"_TM"]) + 273
                
                print("A")
                
                smerni_primer_bind = Seq(result[lprimer_str])
                protismerni_primer_bind = Seq(result[rprimer_str])
                
                tm = min(smerni_primer_tm,protismerni_primer_tm)+3
        
                podaljsan_neobdelan_sprednji_primer = deli[0] + smerni_primer_bind
                podaljsan_neobdelan_protismerni_primer = deli[2].reverse_complement() + protismerni_primer_bind
                
                optimizirana_oliga_lista = optimizacija_konstrukcije_oligov(podaljsan_neobdelan_sprednji_primer,
                                                 smerni_primer_bind,
                                                 podaljsan_neobdelan_protismerni_primer,
                                                 protismerni_primer_bind,
                                                 tm)
                
                
                print("B")
                print(optimizirana_oliga_lista)
                if not optimizirana_oliga_lista is None:
                    
                    podaljsek_smerni = optimizirana_oliga_lista[0][:-len(smerni_primer_bind)]
                    podaljsek_protismerni = optimizirana_oliga_lista[1][:-len(protismerni_primer_bind)]
                    
                    if inicializacija:
                        izvor_nov = globals()["identifikator_števnik"]
                    else:
                        izvor_nov = str(izvor) + "__" + str(globals()["identifikator_števnik"])
                    
                    
                    
                    
                    
                    
                    primer_list_dataframe = pd.DataFrame({
                                                        
                                                        "matrica":[deli[1]],
                                                           
                                                        "smerni_oligo":[optimizirana_oliga_lista[0]],
                                                        "smerni_vezava1":[smerni_primer_bind], 
                                                        "smerni_tm1":[smerni_primer_tm], 
                                                        "smerni_5podaljsek":[podaljsek_smerni], 
                                                        
                                                        
                                                        "protismerni_oligo":[optimizirana_oliga_lista[1]], 
                                                        "protsmerni_vezava1":[protismerni_primer_bind], 
                                                        "protismerni_tm1":[protismerni_primer_tm], 
                                                        "protismerni_5podaljsek":[podaljsek_protismerni],
                                                        
                                                        
                                                        "gradnik":[gradnik],
                                                        "izvor": [izvor_nov],

                                                        "id":[globals()["identifikator_števnik"]],
                                                        "PodaljsanaMatrica": [podaljsek_smerni + deli[1] + podaljsek_protismerni.reverse_complement()],
                                                        
                                                        "PreostalPodaljsekF":[deli[0][:-len(podaljsek_smerni)]],
                                                        "PreostalPodaljsekR":[deli[2][len(podaljsek_protismerni):]],
                                                        "Enak_produkt_id":[""]
                                                        })
                    
                    
                    print(primer_list_dataframe)
                    
                    
                    if f"seznam_primerjev_za_konstrukt:{konstrukt}" in globals():
                        
                        row = primer_list_dataframe.iloc[0]
                        
                        df = globals()[f"seznam_primerjev_za_konstrukt:{konstrukt}"].copy()
                        
                        exists = (
                            (df["smerni_oligo"] == row["smerni_oligo"]) &
                            (df["smerni_vezava1"] == row["smerni_vezava1"]) &
                            (df["protismerni_oligo"] == row["protismerni_oligo"]) &
                            (df["protsmerni_vezava1"] == row["protsmerni_vezava1"]) &
                            (df["PodaljsanaMatrica"] == row["PodaljsanaMatrica"])
                        ).any()
                        
                        print("CV")
                        print(exists)
                        
                        if not exists:
                            podaljšek_že_obstaja = ((df["PodaljsanaMatrica"] == row["PodaljsanaMatrica"]) & (df["matrica"] == row["matrica"])).any()
                            
                            if podaljšek_že_obstaja:
                                
                                primer_list_dataframe["Enak_produkt_id"] = df[(df["PodaljsanaMatrica"] == row["PodaljsanaMatrica"])
                                                                              & (df["matrica"] == row["matrica"])]["id"].min()
                                
                            globals()[f"seznam_primerjev_za_konstrukt:{konstrukt}"] = pd.concat([globals()[f"seznam_primerjev_za_konstrukt:{konstrukt}"],primer_list_dataframe])
                            
                            globals()["lista_id_za_obdelavo"].append(globals()["identifikator_števnik"])
                            globals()["identifikator_števnik"] = globals()["identifikator_števnik"] + 1
                         
                        print(globals()["lista_id_za_obdelavo"])
                        print("..__..")
   
    
                    else:
                        globals()[f"seznam_primerjev_za_konstrukt:{konstrukt}"] = primer_list_dataframe
                        
                        globals()["lista_id_za_obdelavo"].append(globals()["identifikator_števnik"])
                        globals()["identifikator_števnik"] = globals()["identifikator_števnik"] + 1
                        
            except:
               print("ni_primerjev")
               continue
                  
#%%____________________________________________________________________________

globals()["identifikator_števnik"] = 0

#print(globals()["lista_sestavljenih_sekvenc"])



for key, value in globals()["lista_sestavljenih_sekvenc"].items():
    
    safe_key = (
                key.replace(":", "_")
                   .replace("/", "_")
                   .replace("\\", "_")
                   .replace("*", "_")
                   .replace("?", "_")
                   .replace("[", "_")
                   .replace("]", "_")
                   .replace("|","_")
            )
    
    
    
    lokacija_excela = os.path.join(test_program_lok,safe_key + ".xlsx")
    
    

    with pd.ExcelWriter(lokacija_excela, engine="openpyxl") as pisi_excel:
    
#%%___________pridobi_podatke_o_vseh_rezalnih_mestih___________________________

     ime_prvega_res_encima = value[-1][0]
     ime_zadnjega_res_encima = value[-1][1]
     
     Prvi_res_encim = Seq(ime_prvega_res_encima.site)
     Drugi_res_encim = Seq(ime_zadnjega_res_encima.site)
     
     relevantni_lista = value[1].copy()
     relevantni_lista.insert(0,Prvi_res_encim)
     
     relevantni_lista.append(Drugi_res_encim)
     
     relevantno_zaporedje_za_skupino = imena_elementov.copy()
     
     relevantno_zaporedje_za_skupino.insert(0,str(ime_prvega_res_encima))
     relevantno_zaporedje_za_skupino.append(str(ime_zadnjega_res_encima))
     
     #_____________________sestavljanje_delov_______________________________________
     
     i=0
     
     primer_serija = {}
     
     indeksi_delov_s_primerji_v_listi_vseh_delov = []
     lista_vseh_delov = []
     
     sestavljanje_dela = Seq("")
     iterator_za_določevanje_indexov_delov_za_primanje_v_listi_vseh_delov = 0
     
     count_primer_deli = 0
     lista_imen_primer_delov = []
     
     for part in relevantni_lista:
         ime_relevantnega_dela = relevantno_zaporedje_za_skupino[i]
         i = i + 1
         
         if len(part)>= 100:
             
             primer_serija[ime_relevantnega_dela] = [sestavljanje_dela,part]
             lista_imen_primer_delov.append(ime_relevantnega_dela)
             
             if count_primer_deli > 0:
                 prejšnji_primer_del = lista_imen_primer_delov[count_primer_deli-1]
                 
                 primer_serija[prejšnji_primer_del].append(sestavljanje_dela)
                 
     
             count_primer_deli = count_primer_deli+1
             
             if not sestavljanje_dela == Seq(""):
                 iterator_za_določevanje_indexov_delov_za_primanje_v_listi_vseh_delov = iterator_za_določevanje_indexov_delov_za_primanje_v_listi_vseh_delov+1
                 lista_vseh_delov.append(sestavljanje_dela)
                 sestavljanje_dela = Seq("")
     
             indeksi_delov_s_primerji_v_listi_vseh_delov.append(iterator_za_določevanje_indexov_delov_za_primanje_v_listi_vseh_delov)
             
             lista_vseh_delov.append(part)
             iterator_za_določevanje_indexov_delov_za_primanje_v_listi_vseh_delov = iterator_za_določevanje_indexov_delov_za_primanje_v_listi_vseh_delov+1
             
         
         else:
             sestavljanje_dela = sestavljanje_dela + part
             
             if i == len(relevantni_lista):
                 lista_vseh_delov.append(sestavljanje_dela)
     
                 if count_primer_deli > 0:
                     
                     prejšnji_primer_del = lista_imen_primer_delov[count_primer_deli-1]
                     primer_serija[prejšnji_primer_del].append(sestavljanje_dela)
     
     
     
     globals()["identifikator_števnik"] = 0
     for key_gradnik,value_gradnik in primer_serija.items():
         
         
         globals()["lista_id_za_obdelavo"] = []
         globals()["lista_id_za_obdelavo_aktivna"] = []
         globals()["lista_id_obdelani"] = []
         
         inicializacija = True
         potreben = True
         
         konstrukcija_primerjev(value_gradnik,key,key_gradnik,0,inicializacija,potreben)
         
         globals()["lista_id_za_obdelavo_aktivna"] = globals()["lista_id_za_obdelavo"]
    
         globals()["lista_id_za_obdelavo"] = []
         
         
         
         
         inicializacija = False
         for i in range(3):
             print("::--::")
             print(i)
             print(globals()["lista_id_za_obdelavo_aktivna"])
             for id_za_obdelavo in globals()["lista_id_za_obdelavo_aktivna"]:
    
                 maska_za_obdelavo = globals()[f"seznam_primerjev_za_konstrukt:{key}"].loc[ globals()[f"seznam_primerjev_za_konstrukt:{key}"]["id"] == id_za_obdelavo]
    
                 matrica_nova = maska_za_obdelavo["PodaljsanaMatrica"].iloc[0]
                 podaljsek_f = maska_za_obdelavo["PreostalPodaljsekF"].iloc[0]
                 podaljsek_r = maska_za_obdelavo["PreostalPodaljsekR"].iloc[0]
                 deli_nova_matrica = [podaljsek_f,matrica_nova,podaljsek_r]
                 
                 izvor = maska_za_obdelavo["izvor"].iloc[0]
                 
                 #print("__")
                 #print(matrica_nova)
                 #print("--,,")
    
                 if podaljsek_f == Seq("") and podaljsek_r == Seq(""):
                     potreben = False
                 else:
                     potreben = True
                 
                 konstrukcija_primerjev(deli_nova_matrica, key, key_gradnik, izvor,inicializacija,potreben)
                 
                 
                 
             #print(globals()[f"seznam_primerjev_za_konstrukt:{key}"])
             globals()["lista_id_za_obdelavo_aktivna"] = globals()["lista_id_za_obdelavo"]
             globals()["lista_id_za_obdelavo"] = []
           
         #print(globals()[f"seznam_primerjev_za_konstrukt:{key}"])
         
        #print()
         dataframe_gradnik = globals()[f"seznam_primerjev_za_konstrukt:{key}"][
                                                                             globals()[f"seznam_primerjev_za_konstrukt:{key}"]["gradnik"] == key_gradnik
                                                                             ]
    
         
         #print(key_gradnik)
         
         dataframe_gradnik.to_excel(pisi_excel, sheet_name=key_gradnik, index=False)
    
         
             
    
    # print(primer_serija.keys())
    # print(globals()[f"seznam_primerjev_za_konstrukt:{key}"])
     
     kopija_liste_OE_PCR = globals()[f"seznam_primerjev_za_konstrukt:{key}"].copy()
     kopija_liste_OE_PCR = kopija_liste_OE_PCR[(kopija_liste_OE_PCR["PreostalPodaljsekF"] == Seq("")) &
                                            (kopija_liste_OE_PCR["PreostalPodaljsekR"] == Seq(""))
                                            ]
     
     
     print("::")
     print(kopija_liste_OE_PCR)
     print("..")    
     
     
     
     
     
     
         
         
     

 
    