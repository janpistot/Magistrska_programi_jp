
#%% importi
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

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

import random
#%% Inicializacija vsega
root = tk.Tk()
root.title("A2_assembly")
root.geometry("800x800")

#začetne spremenljivke_________________________________________________________
barve_vrstni_red = []

velikost_okna_sestavljanje_konstrukta = 600

#generator barve_______________________________________________________________
for i in range(100):
    # generate a random color in hex format (#RRGGBB)
    color = "#{:02x}{:02x}{:02x}".format(
        random.randint(0, 255),
        random.randint(0, 255),
        random.randint(0, 255)
    )
    barve_vrstni_red.append(color)

#____________________________Konec inicializacije______________________________


#center frame: prostor za prikaz podatkov______________________________________
 
center = tk.Frame(root,height=700)
center.pack(expand=True,fill="both", side="bottom")


#orodna vrstica: widgeti za uvajanje podatkov__________________________________

orodna_vrstica = tk.Frame(root,height=50)
orodna_vrstica.pack(expand=True,fill="both", side="top")


#-->Podatki: meni z možnostmi za uvajanje podatkov v program___________________

podatki = tk.Menubutton(orodna_vrstica, text= "podatki", relief=tk.RAISED)
podatki_meni = tk.Menu(podatki, tearoff=0)
podatki.config(menu=podatki_meni)
podatki.pack(side= "left", fill="both",expand=True, padx= 5, pady= 5)


#------>Komanda_sestavljanje_konstrukta: novo toplevel okno znotraj istega roota za sestavljanje 2A konstrukta
#%% Sestavljanje konstrukta in pripadajoče globalne spremenljivke

globals()["skupine_popolni_record"] = {"vektor":[0,0,""]} #dicitionary z listami, ki vsebujejo vse podatke o uporabljenih sekvencah

#skupine_popolni_record  = {zaporedje:[zaporedje], vektor:[vektor_insert_od,vektor_insert_do,SeqRecord(vektor)],skupina1[seq_record11,seq_record12...],skupina2[seq_record21,seq_record22...]}

pridobljeni_podatki_flag = tk.BooleanVar(value=False)

def uvozi_podatke_gumb_root():
    
    global pridobljeni_podatki_flag

    file_path = filedialog.askopenfilename(title="izberi_datoteko", 
                                           filetypes=[("Pickle files", "*.pkl")])
    #fasta file, prvi heading je zaporedje in nima sekvence, skupine so ločene z "__" (A__B__C za zaporedje ["A","B","C"])
    
    if file_path:
        
        with open(file_path,"rb") as file_s_sekvencami:                   
 
            sk_ime = pickle.load(file_s_sekvencami) 
            
            if sk_ime == "brez_restrikcije":
                globals()["skupine_popolni_record"] = pickle.load(file_s_sekvencami)
                
                pridobljeni_podatki_flag.set(True)

            elif sk_ime == "z_restrikcijo":
                globals()["skupine_popolni_record"] = pickle.load(file_s_sekvencami)
                globals()["lista_sestavljenih_sekvenc"] = pickle.load(file_s_sekvencami)
                
                
                
                pridobljeni_podatki_flag.set(True)
                globals()["obdelava_z_restrikcijskimi_zakljucena"].set(True)
    
def sestavljanje_konstrukta_toplevel(reset,skupine_popolni_record):
    #__________________________________________________________________________
#%% Inicializacija toplevel
    Sestavljanje_konstrukta = tk.Toplevel(root)
    Sestavljanje_konstrukta.title("Sestavljanje konstrukta")
    Sestavljanje_konstrukta.geometry(str(velikost_okna_sestavljanje_konstrukta) +"x" +str(velikost_okna_sestavljanje_konstrukta))
    
    #Frame z vrstico za vnos podatkov
    vrstica_za_sestavljanje = tk.Frame(Sestavljanje_konstrukta,height= velikost_okna_sestavljanje_konstrukta/10)
    vrstica_za_sestavljanje.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    #Frame z vrstico z navodili za vnos podatkov
    vrstica_navodila_vnos_skupin = tk.Frame(Sestavljanje_konstrukta,height= velikost_okna_sestavljanje_konstrukta/10)
    vrstica_navodila_vnos_skupin.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    #Frame z gumbi za vnos sekvenc
    vrstica_vnos_sekvenc= tk.Frame(Sestavljanje_konstrukta,height= velikost_okna_sestavljanje_konstrukta/10)
    vrstica_vnos_sekvenc.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    #Frame z gumboma za naprej in prekliči
    vrstica_gumb_naprej_cancel = tk.Frame(Sestavljanje_konstrukta,height= velikost_okna_sestavljanje_konstrukta/10)
    vrstica_gumb_naprej_cancel.pack(side="bottom",fill="both",expand=True, padx= 5, pady= 5)
    
    #Frame s sliko
    vrstica_slika = tk.Frame(Sestavljanje_konstrukta,height= (velikost_okna_sestavljanje_konstrukta/10) *4)
    vrstica_slika.pack(side="bottom",fill="both",expand=True, padx= 5, pady= 5)
    
    #Frame z gumbom za vektor
    Vektor_import_gumb_frame = tk.Frame(Sestavljanje_konstrukta,height= (velikost_okna_sestavljanje_konstrukta/10))
    Vektor_import_gumb_frame.pack(side="bottom",fill="both",expand=True, padx= 5, pady= 5)
    
    #napis z navodili
    label_skupine = ttk.Label(vrstica_navodila_vnos_skupin, text="Vnesi imena skupin genskih elementov")
    label_skupine.grid(row= 0,column= 0, columnspan=1, sticky="nsew")
    
    zaporedje = []

#%% widgeti za definicijo vektorja    
    
    def vektor_textbox_od_do_komanda(event):
        widget = event.widget
        if widget.edit_modified:
            widget.edit_modified(False)
            
            stolpec = int(widget.grid_info()["column"])
            
            globals()["skupine_popolni_record"]["vektor"][stolpec] = widget.get("1.0","end-1c")

    def vektor_import_gumb_komanda():
        
        file_path = filedialog.askopenfilename(parent = Sestavljanje_konstrukta,title="izberi_datoteko", 
                                               filetypes=[("Text files", "*.fasta")])

        if file_path:
            with open(file_path,"r") as file_s_sekvencami:
                for record in SeqIO.parse(file_s_sekvencami,"fasta"):
                    
                    globals()["skupine_popolni_record"]["vektor"][2] = record
                    
                    vektor_tekst.set(file_path.split("/")[-1])


    vektor_textbox_od = tk.Text(Vektor_import_gumb_frame, width=10, height=5)
    vektor_textbox_do = tk.Text(Vektor_import_gumb_frame, width=10, height=5)
    
    vektor_textbox_od.bind("<<Modified>>",vektor_textbox_od_do_komanda)
    vektor_textbox_do.bind("<<Modified>>",vektor_textbox_od_do_komanda)
    
    vektor_od_label = ttk.Label(Vektor_import_gumb_frame, text="začetno mesto inserta v vektor")
    vektor_do_label = ttk.Label(Vektor_import_gumb_frame, text="končno mesto inserta v vektor")
    
    if isinstance(globals()["skupine_popolni_record"]["vektor"][2],SeqRecord):
        vektor_tekst = tk.StringVar(value= globals()["skupine_popolni_record"]["vektor"][2].id.split("__")[0])
   
    else:
        vektor_tekst = tk.StringVar(value="išči vektor")
    
    vektor_import_gumb = ttk.Button(Vektor_import_gumb_frame,textvariable= vektor_tekst, command= vektor_import_gumb_komanda)
    
    
    vektor_od_label.grid(row= 0,column= 0, columnspan=1, sticky="nsew")
    vektor_do_label.grid(row= 0,column= 1, columnspan=1, sticky="nsew")
    
    vektor_textbox_od.grid(row= 1,column= 0, columnspan=1, sticky="nsew")
    vektor_textbox_od.insert("1.0", str(globals()["skupine_popolni_record"]["vektor"][0]))
    
    vektor_textbox_do.grid(row= 1,column= 1, columnspan=1, sticky="nsew")
    vektor_textbox_do.insert("1.0", str(globals()["skupine_popolni_record"]["vektor"][1]))
    
    vektor_import_gumb.grid(row= 1,column= 2, columnspan=1, sticky="nsew")
    
    

#%% widgeti za shranjevanje in preklic v toplevel    
    def shrani_kot_command():
        nonlocal zaporedje
        globals()["skupine_popolni_record"]["zaporedje"] = zaporedje
                
        shrani_kot = filedialog.asksaveasfilename(
            title="shrani kot ...",
            defaultextension=".pkl",
            filetypes=[("Pickle files", "*.pkl")],
            initialdir=".",          # starting directory
            initialfile="datoteka"   # suggested filename
        )
        
        if not shrani_kot:
            return
        
        with open(shrani_kot,"ab") as sk:
            pickle.dump("brez_restrikcije",sk)
            pickle.dump(globals()["skupine_popolni_record"],sk)
           
           
        pridobljeni_podatki_flag.set(True)
        Sestavljanje_konstrukta.destroy()
   
    def gumb_naprej_command():
        nonlocal zaporedje
        globals()["skupine_popolni_record"]["zaporedje"] = zaporedje
        pridobljeni_podatki_flag.set(True)

        Sestavljanje_konstrukta.destroy()
  
    #ttk_spremenljivki_gumbov_naprej/prekliči
    gumb_naprej = ttk.Button(vrstica_gumb_naprej_cancel,text ="Potrdi izbiro", command= gumb_naprej_command)
    gumb_shrani_kot = ttk.Button(vrstica_gumb_naprej_cancel,text ="Shrani_kot", command= shrani_kot_command)
    gumb_prekliči = ttk.Button(vrstica_gumb_naprej_cancel,text ="prekliči",command= lambda: Sestavljanje_konstrukta.destroy())
    
    #pozicija_gumbov_naprej_prekliči   	     
    gumb_naprej.pack(side="right",fill="both",expand=True, padx= 5, pady= 5)    
    gumb_prekliči.pack(side="left",fill="both",expand=True, padx= 5, pady= 5)
    
    gumb_shrani_kot.pack(side="left",fill="both",expand=True, padx= 5, pady= 5)
    
    global pridobljeni_podatki_flag
    

# %% gumbi v toplevel
    
    zadnji_stolpec = 0
    textbox_list = []
    gumbi_sekvence_list =[]
    
    povprečje_list = [0]
    
    tekstbox = tk.Text(vrstica_za_sestavljanje, width=10, height=5)
    tekstbox.grid(row= 0,column= zadnji_stolpec, sticky="nsew", padx= 2, pady= 2)
    textbox_list.append(tekstbox)  
    
#%% tekstboxi
    if reset:
        pridobljeni_podatki_flag.set(False)
        globals()["skupine_popolni_record"] = {"vektor":[0,0,""]}
        
    def tekstbox_komanda(event):
        widget = event.widget
        zadnji_stolpec = len(textbox_list) -1

        def delanje_slik(povprečje_list, lista_textbox_podatkov):
                
                count = 0
                suma_dolžin_slik = 0 # inicaializacija spremenljivke, ki spremlja dolžino med izgradnjo slike za določanje pozicije novega elementa
                skupna_dolžina_elementov = sum(povprečje_list) # skupna dolžina za določanje relatvnih dolžin fragmentov glede na velikost okna
                
                
                for element in povprečje_list:
                    
                    if not element == 0: #zadnji dodan je vedno prvotno 0, da lahko spreminjam elemente v kateremkoli vrstnem redu
                        try:
                            element_s_prilagojeno_dolžino = (element/skupna_dolžina_elementov) * velikost_okna_sestavljanje_konstrukta
                        
                        except:
                            element_s_prilagojeno_dolžino = 0

                        napis = lista_textbox_podatkov[count]
                        
                        locals()[napis +"__slika"] = tk.Label(vrstica_slika, text=napis, bg= barve_vrstni_red[count])                
                        locals()[napis +"__slika"].place(x=suma_dolžin_slik, y= 0, width=element_s_prilagojeno_dolžino, relheight=0.3)
                                        
                        suma_dolžin_slik = suma_dolžin_slik + element_s_prilagojeno_dolžino
                        
                    count = count + 1
        
        if widget.edit_modified():
            widget.edit_modified(False)
            
            stolpec = int(widget.grid_info()["column"])
            
            
            if zadnji_stolpec == stolpec:

                povprečje_list.append(0)
                zaporedje.append("")
                
                
                nov_tekstbox = tk.Text(vrstica_za_sestavljanje, width=10, height=5)
                nov_tekstbox.grid(row= 0,column= zadnji_stolpec + 1, sticky="nsew", padx= 2, pady= 2)
                textbox_list.append(nov_tekstbox)

                nov_tekstbox.bind("<<Modified>>",tekstbox_komanda)
#%%gumbi za iskanje sekvenc
                def iskanje_sekvenc():
                    nonlocal povprečje_list
                    nonlocal stolpec
                    nonlocal zaporedje
                    
                    
                    gumb_sekvence_full_path= tk.StringVar()
                    
                    
                    file_path = filedialog.askopenfilename(parent = Sestavljanje_konstrukta,
                                                           title="izberi_datoteko",
                                                           filetypes=[("Text files", "*.fasta")])
                    
                    if file_path:
                        
                        gumb_išči_sekvence_tekst.set(file_path.split("/")[-1])
                        gumb_sekvence_full_path.set(file_path)
                        
                        with open(file_path,"r") as sekvenca_file:
                            vse_seq_dolžine =[]
                            
                            for record in SeqIO.parse(sekvenca_file, "fasta"):
                                vse_seq_dolžine.append(len(record))

                                ime_skupine = widget.get("1.0","end-1c")
                                
                                if ime_skupine in globals()["skupine_popolni_record"].keys():
                                    globals()["skupine_popolni_record"][ime_skupine].append(record)
                                else:
                                    globals()["skupine_popolni_record"][ime_skupine] = [record]
                                    
                            povprečna_dolžina_elementa = sum(vse_seq_dolžine)/len(vse_seq_dolžine)
                            
                       
                        povprečje_list[stolpec] = povprečna_dolžina_elementa
                        
                        delanje_slik(povprečje_list, zaporedje)                      

   
                gumb_išči_sekvence_tekst = tk.StringVar(value="išči sekvenco") 
                
                gumb_išči_sekvence = tk.Button(vrstica_vnos_sekvenc,textvariable= gumb_išči_sekvence_tekst, command=iskanje_sekvenc)                
                gumb_išči_sekvence.grid(row= 3,column=stolpec, sticky="nsew", padx= 2, pady= 2)
                
                gumbi_sekvence_list.append(gumb_išči_sekvence)              

#%% brisanje stolpcev, če se zadnji izprazne 
           
            if widget.get("1.0","end-1c") == "" and zadnji_stolpec-1 == stolpec:
                
                textbox_list[-1].destroy()
                gumbi_sekvence_list[-1].destroy()
                
                del textbox_list[-1]
                del gumbi_sekvence_list[-1]
                del povprečje_list[-1]
                del zaporedje[-1]

            zaporedje[stolpec] = widget.get("1.0","end-1c")
            delanje_slik(povprečje_list, zaporedje)
           
    tekstbox.bind("<<Modified>>",tekstbox_komanda)
    

#%% uvajanje uvoženih zaporedij v sestavljanje konstrukta UI
    
    if not reset and "zaporedje" in globals()["skupine_popolni_record"] and pridobljeni_podatki_flag:
        
        pridobljeni_podatki_flag.set(False)
        
        def uvajanje_uvoženih_zaporedij(i=0):
            nonlocal povprečje_list            
            zaporedje = globals()["skupine_popolni_record"]["zaporedje"]

            if i < len(zaporedje):

                element = zaporedje[i]
                textbox_list[i].insert("1.0", element)
                Sestavljanje_konstrukta.after(0, uvajanje_uvoženih_zaporedij, i + 1)

                povprečje_temp_lista =[]
                record_temp_lista = []
 
                for record in globals()["skupine_popolni_record"][element]:
                    
                    povprečje_temp_lista.append(len(record.seq))
                    record_temp_lista.append(record)

                povprečje_skupina = sum(povprečje_temp_lista)/len(povprečje_temp_lista)

                povprečje_list[i] = povprečje_skupina

                  
        uvajanje_uvoženih_zaporedij()
        
        if not pridobljeni_podatki_flag.get():
            pridobljeni_podatki_flag.set(True)
    

#%% dodajanje gumbov in spremljanje sprememb

podatki_meni.add_command(label="uvozi podatke", command= uvozi_podatke_gumb_root)
podatki_meni.add_command(label="sestavi konstrukt", command= lambda reset = True: sestavljanje_konstrukta_toplevel(reset,globals()["skupine_popolni_record"]))

globals()["obdelava_z_restrikcijskimi_zakljucena"] = tk.BooleanVar(value= False)

def vektor_restrikcija_komanda():
    #%% sestavljanje polnih konstruktov brez podaljškov: pomembno za restrikcijo
    
    Vektor_restrikcija_toplevel = tk.Toplevel(root)
    Vektor_restrikcija_toplevel.title("Sestavljanje konstrukta")
    Vektor_restrikcija_toplevel.geometry(str(velikost_okna_sestavljanje_konstrukta*2) +"x" +str(velikost_okna_sestavljanje_konstrukta))
    
    restrikcijski_encimi_določeni_check_frame = tk.Frame(Vektor_restrikcija_toplevel,height= velikost_okna_sestavljanje_konstrukta, width= 20)
    restrikcijski_encimi_določeni_check_frame.pack(side="right",fill="both", expand=True, padx= 5, pady= 5)
    
    restrikcijski_encimi_določeni_check_frame_podframe_oznake = tk.Frame(restrikcijski_encimi_določeni_check_frame,height= 0.5 * velikost_okna_sestavljanje_konstrukta, width= 20)
    restrikcijski_encimi_določeni_check_frame_podframe_oznake.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    restrikcijski_encimi_določeni_check_frame_podframe_gumb_naprej = tk.Frame(restrikcijski_encimi_določeni_check_frame,height= 0.5 * velikost_okna_sestavljanje_konstrukta, width= 20)
    restrikcijski_encimi_določeni_check_frame_podframe_gumb_naprej.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    #navodila
    label_z_navodili_frame = tk.Frame(Vektor_restrikcija_toplevel,height= velikost_okna_sestavljanje_konstrukta/10)
    label_z_navodili_frame.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    navodila_label = ttk.Label(label_z_navodili_frame, text="Vnesi restrikcijske encime, ki so na voljo ali prenesi seznam s temi encimi")
    navodila_label.grid(row= 0,column= 0, columnspan=1, sticky="nsew")
    
    #tekstboxi za vnos podatkov
    combobox_frame = tk.Frame(Vektor_restrikcija_toplevel,height= velikost_okna_sestavljanje_konstrukta/10)
    combobox_frame.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    shrani_restrikcijske_encime_frame = tk.Frame(Vektor_restrikcija_toplevel,height= velikost_okna_sestavljanje_konstrukta/10)
    shrani_restrikcijske_encime_frame.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    slika_vektor_frame = tk.Frame(Vektor_restrikcija_toplevel,height= velikost_okna_sestavljanje_konstrukta/10 * 4)
    slika_vektor_frame.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    combobox_izbira_ustreznih_encimov_frame = tk.Frame(Vektor_restrikcija_toplevel,height= velikost_okna_sestavljanje_konstrukta/10)
    combobox_izbira_ustreznih_encimov_frame.pack(side="bottom",fill="both", expand=True, padx= 5, pady= 5)
    
    
    Restrikcijski_encimi = list(CommOnly)
    Imena_restrikcijskih_encimov = []
    
    flag_podatki_res_encimi_pridobljeni = tk.BooleanVar(value= False)
    
    
    
    for Restrikcijski_encim in Restrikcijski_encimi:
        
        Imena_restrikcijskih_encimov.append(Restrikcijski_encim.__name__)
    
    #%%comboboxi z restrikcijskimi encimi
    Imena_restrikcijskih_encimov = sorted(Imena_restrikcijskih_encimov, key=str.lower)

    začetni_combobox_za_izbiro_res_encimov = ttk.Combobox(combobox_frame, values= Imena_restrikcijskih_encimov)
    začetni_combobox_za_izbiro_res_encimov.grid(row= 0,column= 0, columnspan=1, sticky="nsew")
    
    combobox_lista = [začetni_combobox_za_izbiro_res_encimov]
    
    Lista_encimov_na_voljo =[""]
    
    določeni_restrikcijski_encimi_dic = {}
    
    def restrikcija_z_uporabnimi_encimi():
        
        nonlocal Lista_encimov_na_voljo
        
        Lista_encimov_na_voljo_kopija =  [x for x in Lista_encimov_na_voljo if x != ""]
        uporabni_restrikcijski_encimi = RestrictionBatch(Lista_encimov_na_voljo_kopija)
        
        print(globals()["skupine_popolni_record"]["vektor"])
        sekvenca_vektorja = globals()["skupine_popolni_record"]["vektor"][2].seq

        območje_vstavitve_v_vektor = [globals()["skupine_popolni_record"]["vektor"][0],globals()["skupine_popolni_record"]["vektor"][1]]
        
        vektor_analiza = Analysis(uporabni_restrikcijski_encimi, sekvenca_vektorja)
        rezultati_vektor = vektor_analiza.full()

        uporabni_encimi_vektor = []
        uporabni_encimi_vektor_lokacije ={}
        
        for k, v in rezultati_vektor.items():
            break_flag=True
            for mesto in v:                
                if not (int(območje_vstavitve_v_vektor[0]) <= mesto <= int( območje_vstavitve_v_vektor[1])) and break_flag:
                    break_flag=False
                    break
            
            if break_flag and not v ==[]:
                
                uporabni_encimi_vektor.append(k)
                uporabni_encimi_vektor_lokacije[k] = v

        uporabni_encimi_fragment = {}

        
        for key, value in  globals()["lista_sestavljenih_sekvenc"].items():
            sekvenca_združenega_fragmenta = value[0]

            
            insert_analiza = Analysis(uporabni_restrikcijski_encimi, sekvenca_združenega_fragmenta)
            rezultati_insert = insert_analiza.full()
            
            lista_uporabnih_encimov_za_fragment = []
            
            for k, v in rezultati_insert.items():
                if v == []:
                    lista_uporabnih_encimov_za_fragment.append(k)
                    
            uporabni_encimi_fragment[key] = lista_uporabnih_encimov_za_fragment.copy()
        
        for key,value in uporabni_encimi_fragment.items():
            dejanski_uporabni_encimi_za_fragment =[]
            encimi_lokacije = []
            
            for encim_insert in value:
                
                if encim_insert in uporabni_encimi_vektor:
                    dejanski_uporabni_encimi_za_fragment.append(encim_insert)
                    
                    encimi_lokacije.append(uporabni_encimi_vektor_lokacije[encim_insert])
                

            if len(globals()["lista_sestavljenih_sekvenc"][key]) < 4:
                globals()["lista_sestavljenih_sekvenc"][key].append("")
                globals()["lista_sestavljenih_sekvenc"][key].append("")
                
            globals()["lista_sestavljenih_sekvenc"][key][2] = dejanski_uporabni_encimi_za_fragment
            globals()["lista_sestavljenih_sekvenc"][key][3] = encimi_lokacije
            
    slika_checkbox_label_list =[]
    slika_checkbox_res_mesta_label_list = []
    lista_checkov_za_restrikcijska_mesta = []
    
    def Gumb_vse_doloceno_komanda():
        for key, value in  globals()["lista_sestavljenih_sekvenc"].items():
            for k, v in določeni_restrikcijski_encimi_dic.items():
                if key == k and not v == []:
                    try:
                        value[4] = v
                    except:
                        value.append(v)

        globals()["obdelava_z_restrikcijskimi_zakljucena"].set(True)
        
        
        Vektor_restrikcija_toplevel.destroy()
    
    def Gumb_shrani_vse_doloceno_kot_komanda():
        nonlocal Lista_encimov_na_voljo,checkbox_check
        
        Lista_encimov_na_voljo_kopija = Lista_encimov_na_voljo.copy()
        
        for key, value in  globals()["lista_sestavljenih_sekvenc"].items():
            for k, v in določeni_restrikcijski_encimi_dic.items():
                if key == k and not v == []:
                    try:
                        value[4] = v
                    except:
                        value.append(v)
        
        shrani_kot = filedialog.asksaveasfilename(
            title="shrani kot ...",
            defaultextension=".pkl",
            filetypes=[("Pickle files", "*.pkl")],
            initialdir=".",          
            initialfile="datoteka",   
            parent= Vektor_restrikcija_toplevel,
        )
        
        
        if not shrani_kot:
            return
        
        #pikl
        with open(shrani_kot,"ab") as sk:
            pickle.dump("z_restrikcijo",sk)
            pickle.dump(globals()["skupine_popolni_record"],sk)
            pickle.dump(globals()["lista_sestavljenih_sekvenc"],sk)
    
            
        globals()["obdelava_z_restrikcijskimi_zakljucena"].set(True)
        Vektor_restrikcija_toplevel.destroy()
        
    
    Gumb_vse_doloceno = tk.Button(restrikcijski_encimi_določeni_check_frame_podframe_gumb_naprej,
                                  text= "Vsem kombinacijam še niso bili določeni ustrezni restrikcijski encimi", 
                                  command= Gumb_vse_doloceno_komanda,
                                  state="disabled")
    
    Gumb_vse_doloceno.pack(side="left",fill="both",expand=True, padx= 5, pady= 5)
    
    Gumb_shrani_vse_doloceno_kot = tk.Button(restrikcijski_encimi_določeni_check_frame_podframe_gumb_naprej,
                                             text= "Shrani določene encime kot", 
                                             command= Gumb_shrani_vse_doloceno_kot_komanda,
                                             state="disabled")
    
    Gumb_shrani_vse_doloceno_kot.pack(side="left",fill="both",expand=True, padx= 5, pady= 5)

    
    def slika_checkbox_komanda(event):
        nonlocal slika_checkbox_label_list, slika_checkbox_res_mesta_label_list            
        
        for label in slika_checkbox_res_mesta_label_list:
            label.destroy()
            
        for label in slika_checkbox_label_list:
            label.destroy()
        
        slika_checkbox_label_list =[]
        slika_checkbox_res_mesta_label_list = []
        
        
        izbran_sestavljen_konstrukt = vektor_slika_checkbox.get()
        
        if not izbran_sestavljen_konstrukt == "izberi si konstrukt, ki ga želiš analizirati":
            restrikcijski_encimi = globals()["lista_sestavljenih_sekvenc"][izbran_sestavljen_konstrukt][2]
            lokacije_rezov = globals()["lista_sestavljenih_sekvenc"][izbran_sestavljen_konstrukt][3]
            ime_vektorja = globals()["skupine_popolni_record"]["vektor"][2].id
            dolžina_sekvence_vektorja = len(globals()["skupine_popolni_record"]["vektor"][2].seq)
            izbran_MCS = [int(globals()["skupine_popolni_record"]["vektor"][0]), int(globals()["skupine_popolni_record"]["vektor"][1])]
    
            en_nukleotid_dolžina_na_zaslonu = velikost_okna_sestavljanje_konstrukta/dolžina_sekvence_vektorja
    
            slika_vektor = tk.Label(slika_vektor_frame, text= ime_vektorja,bg= barve_vrstni_red[0])
            slika_vektor.place(x= 0, y=0, width= velikost_okna_sestavljanje_konstrukta, height= 30)
            
            slika_MCS = tk.Label(slika_vektor_frame, text= "izbran_MCS", bg= barve_vrstni_red[1])
            slika_MCS.place(x=round(izbran_MCS[0]*en_nukleotid_dolžina_na_zaslonu), y=0, width= round((izbran_MCS[1]-izbran_MCS[0])*en_nukleotid_dolžina_na_zaslonu), height= 30)
            
            slika_checkbox_label_list.append(slika_vektor)
            slika_checkbox_label_list.append(slika_MCS)
            
            mesto_ne_liste = []
            count = 0
            
            for mesto_restrikcijski_encim in lokacije_rezov:
                count1 = 0
                for lokacija_rezov in mesto_restrikcijski_encim:
                    mesto_ne_liste.append(lokacija_rezov)
                    locals()["slika_restrikcijsko_mesto_"+str(count1)+ "_" + restrikcijski_encimi[count].__name__] = tk.Label(slika_vektor_frame, text= "", bg= barve_vrstni_red[count+3])
                    locals()["slika_restrikcijsko_mesto_"+str(count1)+ "_" + restrikcijski_encimi[count].__name__].place(x= round(int(lokacija_rezov) *en_nukleotid_dolžina_na_zaslonu) , y=0, width=1, height= 30)
                    
                    slika_checkbox_label_list.append(locals()["slika_restrikcijsko_mesto_"+str(count1)+ "_" + restrikcijski_encimi[count].__name__])
                    
                    count1 = count1 +1
                    
                    counter_y = 1
                    for i in mesto_ne_liste:
                        if i-200 < lokacija_rezov < i+200:
                            counter_y =counter_y + 1
                    
                    locals()["restrikcijsko_mesto_tag"+restrikcijski_encimi[count].__name__+ "_" + str(lokacija_rezov)] = tk.Label(slika_vektor_frame, text= restrikcijski_encimi[count].__name__ + "__" + str(lokacija_rezov), font=("Arial", 5))
                    locals()["restrikcijsko_mesto_tag"+restrikcijski_encimi[count].__name__+ "_" + str(lokacija_rezov)].place(x= round(int(lokacija_rezov) *en_nukleotid_dolžina_na_zaslonu), y=counter_y*8+10, height= 5)
                    
                    slika_checkbox_res_mesta_label_list.append(locals()["restrikcijsko_mesto_tag"+restrikcijski_encimi[count].__name__+ "_" + str(lokacija_rezov)])
                count = count +1  
            
             
                      
            widget = event.widget
            comobox_tekst = widget.get()
            
            
            izbrana_encima = ["",""]
            lokacija_ena_lista = []
            lista_encimov_s_številom_lokacij = []
            
            listi_možnih_vnosov=[globals()["lista_sestavljenih_sekvenc"][comobox_tekst][2],""]
            
            def combobox_encim_določanje(event,combobox_index):
                
                nonlocal comobox_tekst, slika_checkbox_res_mesta_label_list,lokacija_ena_lista,lista_encimov_s_številom_lokacij,combobox_res_encim1,combobox_res_encim2,listi_možnih_vnosov
                
                lista_encimov_s_številom_lokacij = []
                
                
                combobox_encimi_izbira = event.widget
                
                cb_text = combobox_encimi_izbira.get()
    
                izbrana_encima[combobox_index] = cb_text
                
                lista_uporabnih_encimov = [str(encim) for encim in globals()["lista_sestavljenih_sekvenc"][comobox_tekst][2]]
                lokacije_uporabnih_encimov = globals()["lista_sestavljenih_sekvenc"][comobox_tekst][3]
                
                
                index = -1            
                for lokacije in lokacije_uporabnih_encimov:
                    index = index+ 1
                    for lokacija in lokacije:
                        lista_encimov_s_številom_lokacij.append(lista_uporabnih_encimov[index])
                        lokacija_ena_lista.append(lokacija)
                
                index = -1
                for slika in slika_checkbox_res_mesta_label_list:
                    index = index + 1
                    if not lista_encimov_s_številom_lokacij[index] in izbrana_encima:
                        slika.configure(bg= "#F0F0F0")
                    
                    else:
                        slika_checkbox_res_mesta_label_list[index].configure(bg= "light green")
    
                
    
                widget = event.widget
                
                if widget.get() not in [encim.__name__ for encim in globals()["lista_sestavljenih_sekvenc"][comobox_tekst][2]]:
                    izbrana_encima[combobox_index] = ""
                
                lista_najdenih_encimov = []
                lista_encimov_ne_ok =  []
                
                if not lokacija_ena_lista == []:
                    index_encimi_1 = -1
                    lokacije_izbran_encim_1 = []
        
                        
                    for enz in lista_encimov_s_številom_lokacij:
                        
                        index_encimi_1 = index_encimi_1 + 1
                        
                        
                        if enz == widget.get():
                            
                            lokacija_reza_izbran_encim_1 = lokacija_ena_lista[index_encimi_1]
                            lokacije_izbran_encim_1.append(lokacija_reza_izbran_encim_1)                    
                            
                    
                    for encim in globals()["lista_sestavljenih_sekvenc"][comobox_tekst][2]:
    
                        index_encima_v_zaporedju_lokacij = lista_encimov_s_številom_lokacij.index(str(encim.__name__))
                        
                        for lokacija_1 in lokacije_izbran_encim_1:
                            if combobox_index == 1:
                                if lokacija_ena_lista[index_encima_v_zaporedju_lokacij] >= lokacija_1:
                                    
                                    lista_encimov_ne_ok.append(encim.__name__)
                                    
                            elif combobox_index == 0:
                                if lokacija_ena_lista[index_encima_v_zaporedju_lokacij] <= lokacija_1:
                                    
                                    lista_encimov_ne_ok.append(encim.__name__)
      
                else:
                    lista_encimov_ne_ok = []
                
                
                for encim in globals()["lista_sestavljenih_sekvenc"][comobox_tekst][2]:   
                    if str(encim.__name__) not in lista_encimov_ne_ok:
                        lista_najdenih_encimov.append(encim)
                    
                
                if widget == combobox_res_encim1:
                    combobox_res_encim2.config(value=lista_najdenih_encimov)
                    listi_možnih_vnosov[1]=lista_najdenih_encimov
                    
                elif widget == combobox_res_encim2:
                    combobox_res_encim1.config(value=lista_najdenih_encimov)
                    listi_možnih_vnosov[0]=lista_najdenih_encimov
                
            
            
            def combobox_lista_config(event,combobox_index):
                
                nonlocal listi_možnih_vnosov,slika_checkbox_res_mesta_label_list,lista_encimov_s_številom_lokacij,izbrana_encima
                
                widget_config= event.widget
                cbb_text = widget_config.get()
                
                izbrana_encima[combobox_index]= cbb_text

                lista_najdenih_encimov = []

                for encim in listi_možnih_vnosov[combobox_index]:
                    if cbb_text == "":
                        lista_najdenih_encimov = listi_možnih_vnosov[combobox_index]
                    
                    elif (cbb_text in str(encim)):
                        lista_najdenih_encimov.append(encim)

                        if cbb_text == encim.__name__:
                            
                            widget_config.event_generate("<<ComboboxSelected>>")
                        
                        else:
                            index = -1
                            
                            for encim_lok in lista_encimov_s_številom_lokacij:
                                index = index + 1
                                
                                if not encim_lok in izbrana_encima:
                                    slika_checkbox_res_mesta_label_list[index].config(bg= "#F0F0F0")
                                    
                 
                                
                                                        
                widget_config.config(value=lista_najdenih_encimov)               
                    

            def gumb_izbira_encimov_command():
                nonlocal vektor_slika_checkbox,combobox_konstrukt_tekst, določeni_restrikcijski_encimi_dic
                
                neurejeni_konstrukti = []
                for key, value in globals()["lista_sestavljenih_sekvenc"].items():
                    
                    encimi_list =[]
                    for vredu_encim in value[2]:
                        for izbran_encim in izbrana_encima:
                            
                            if vredu_encim.__name__ == izbran_encim:
                                encimi_list.append(vredu_encim)
                                break
                    
                    if len(encimi_list) == 2 and not (key in list(določeni_restrikcijski_encimi_dic.keys())):
                        
                        določeni_restrikcijski_encimi_dic[key] = encimi_list
                        
                        vsi_kljuci = list(globals()["lista_sestavljenih_sekvenc"].keys())
                        lista_checkov_za_restrikcijska_mesta[vsi_kljuci.index(key)].configure(bg= "green") 
                    

                for key in globals()["lista_sestavljenih_sekvenc"].keys():
                    if not (key  in določeni_restrikcijski_encimi_dic.keys()):
                        neurejeni_konstrukti.append(key)
                        break
                        
                        
                if neurejeni_konstrukti == []:
                    Gumb_vse_doloceno.configure(state= "normal", text= "Vsem kombinacijam so bili določeni ustrezni encimi: Nadaljuj")
                    Gumb_shrani_vse_doloceno_kot.configure(state= "normal")
                    
                combobox_konstrukt_tekst.set("konstrukt urejen")
                vektor_slika_checkbox.configure(values= neurejeni_konstrukti)
                
            
            if encimi_pridobljeni.get():
    
                combobox_res_encim1 = ttk.Combobox(combobox_izbira_ustreznih_encimov_frame, values= globals()["lista_sestavljenih_sekvenc"][comobox_tekst][2])
                combobox_res_encim1.grid(row= 0,column=0, columnspan=1, sticky="nsew")
                
                combobox_res_encim1.bind("<<ComboboxSelected>>", lambda event : combobox_encim_določanje(event, 0))
                combobox_res_encim1.bind("<KeyRelease>", lambda event:combobox_lista_config(event,0),  add="+")
                
                combobox_res_encim2 = ttk.Combobox(combobox_izbira_ustreznih_encimov_frame, values= globals()["lista_sestavljenih_sekvenc"][comobox_tekst][2])
                combobox_res_encim2.grid(row= 0,column=1, columnspan=1, sticky="nsew")
                
                combobox_res_encim2.bind("<<ComboboxSelected>>", lambda event : combobox_encim_določanje(event, 1))
                combobox_res_encim2.bind("<KeyRelease>", lambda event:combobox_lista_config(event,1),  add="+")

                gumb_izbira_encimov = tk.Button(combobox_izbira_ustreznih_encimov_frame,text= "Potrdi izbiro restrikcijskih encimov", command= gumb_izbira_encimov_command)
                gumb_izbira_encimov.grid(row= 0,column=2, columnspan=1, sticky="nsew")
                
    encimi_pridobljeni = tk.BooleanVar(value= False)
    
    def combobox_komanda(event):
        
        nonlocal Imena_restrikcijskih_encimov
        nonlocal Lista_encimov_na_voljo
        nonlocal combobox_lista
        nonlocal checkbox_check
        widget = event.widget
        
        if Lista_encimov_na_voljo ==[x.__name__ for x in CommOnly]:
            Lista_encimov_na_voljo = [""]
        
        encimi_pridobljeni.set(True)
        
        stolpec = int(widget.grid_info()["column"])
        vnos_v_combobox = widget.get()
        
        nov_drop_list = []
        for encim in Imena_restrikcijskih_encimov:
            if vnos_v_combobox in encim:
                nov_drop_list.append(encim)
                
        widget.config(values= nov_drop_list)
        
        widget.after(0, lambda: widget.focus_set())
        
        if widget.get() in Imena_restrikcijskih_encimov:
            Lista_encimov_na_voljo[stolpec] = vnos_v_combobox
        
        restrikcija_z_uporabnimi_encimi()
        
        if (len(combobox_lista) >= 2):           
            if (widget == combobox_lista[-2]) and (widget.get() == ""):
                combobox_lista[-1].destroy()
                del combobox_lista[-1]
            
        vektor_slika_checkbox.event_generate("<<ComboboxSelected>>")          
        
        if widget == combobox_lista[-1] and not widget.get() == "":
                            
            if not Lista_encimov_na_voljo[-1] == "":
    
                Lista_encimov_na_voljo.append("")
            nov_combobox = ttk.Combobox(combobox_frame, values= Imena_restrikcijskih_encimov)
            nov_combobox.grid(row= 0,column= stolpec + 1, columnspan=1, sticky="nsew")
            
            nov_combobox.bind("<KeyRelease>",combobox_komanda)
            nov_combobox.bind("<<ComboboxSelected>>", combobox_komanda)

            combobox_lista.append(nov_combobox)
        
    začetni_combobox_za_izbiro_res_encimov.bind("<KeyRelease>",combobox_komanda)
    začetni_combobox_za_izbiro_res_encimov.bind("<<ComboboxSelected>>", combobox_komanda)
    
    
    #%% shrani podatke v txt file, ena vrstica, v kateri so encimi ločeni z "__"
    def shrani_restrikcijske_encime_kot():
        nonlocal Lista_encimov_na_voljo
        Lista_encimov_na_voljo_kopija = [x for x in Lista_encimov_na_voljo if x != ""]
        
        shrani_kot = filedialog.asksaveasfilename(
            title="shrani kot ...",
            parent= Vektor_restrikcija_toplevel,
            defaultextension=".txt",
            filetypes=[("Text files")],
            initialdir=".",          
            initialfile="datoteka"
        )
        
        if not shrani_kot:
            return
        
        with open(shrani_kot,"a") as sk:
            tekst_shrani_kot = "__".join(Lista_encimov_na_voljo_kopija[:-1])
            sk.write(tekst_shrani_kot)
    
    gumb_shrani_restrikcijske_encime_kot = tk.Button(shrani_restrikcijske_encime_frame,text= "shrani listo restrikcijskih encimov kot", command= shrani_restrikcijske_encime_kot)
    gumb_shrani_restrikcijske_encime_kot.grid(row= 0,column= 0, columnspan=1, sticky="nsew")
        

    #%%upload že vnešenih restrikcijskih encimov
    lista_encimov_upload_gumb_frame = tk.Frame(Vektor_restrikcija_toplevel,height= velikost_okna_sestavljanje_konstrukta/10)
    lista_encimov_upload_gumb_frame.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    gumb_išči_sekvence_tekst = tk.StringVar(value="vnesi shranjen set restrikcijskih encimov na voljo") 
    
    #encimi ločeni z "__" v eni vrstici txt fil-a.
    def komanda_lista_encimov_upload_gumb():
        nonlocal checkbox_check
            
        encimi_pridobljeni.set(True)
        
        nonlocal Lista_encimov_na_voljo
        nonlocal gumb_išči_sekvence_tekst
        
        file_path = filedialog.askopenfilename(parent = Vektor_restrikcija_toplevel,
                                               title="izberi_datoteko",
                                               filetypes=[("Text files","*.txt")])
        if file_path:
                            
            gumb_išči_sekvence_tekst.set(file_path.split("/")[-1])
        
            with open(file_path,"r") as sekvence_upload_file:
                for line in sekvence_upload_file:
                    Lista_encimov_na_voljo = line.split("__")
                    
            
            Lista_encimov_na_voljo = [x for x in Lista_encimov_na_voljo if x != ""]
               
            def rekurzivna_funkcija_za_vnašanje_encimov_v_comboboxe(i = 0):
                nonlocal combobox_lista
                nonlocal Lista_encimov_na_voljo
    
                tekst_insert = Lista_encimov_na_voljo[i]
                trenutni_combobox = combobox_lista[i]
                
                trenutni_combobox.set("")
                trenutni_combobox.insert(0, tekst_insert)
    
                event = tk.Event()
                event.widget = trenutni_combobox
                combobox_komanda(event)            
                
                if i < len(Lista_encimov_na_voljo)-1:
                    Vektor_restrikcija_toplevel.after(0,rekurzivna_funkcija_za_vnašanje_encimov_v_comboboxe,i+1)
                    
            
            for combobox in combobox_lista[1:]:
                combobox.destroy()
            
            del(combobox_lista[1:])
                
            rekurzivna_funkcija_za_vnašanje_encimov_v_comboboxe()
            restrikcija_z_uporabnimi_encimi()
         
                    
    lista_encimov_upload_gumb = tk.Button(lista_encimov_upload_gumb_frame,textvariable= gumb_išči_sekvence_tekst, command=komanda_lista_encimov_upload_gumb)                
    lista_encimov_upload_gumb.grid(row= 0,column=0, sticky="nsew", padx= 2, pady= 2)
        
    checkbox_check = tk.IntVar()
    vsi_encimi_checkbox = tk.Checkbutton(lista_encimov_upload_gumb_frame, text="namesto liste encimov uporabi vse encime v repozitoriju", variable=checkbox_check)
    vsi_encimi_checkbox.grid(row= 0,column=1, sticky="nsew", padx= 2, pady= 2)
    
    def vsi_encimi_uporaba_checkbox(*args):
        
        nonlocal Lista_encimov_na_voljo, combobox_lista, vektor_slika_checkbox
                
        if checkbox_check.get() == 1:
            
            encimi_pridobljeni.set(True)
            
            Lista_encimov_na_voljo = []
            
            for encim in CommOnly:
                Lista_encimov_na_voljo.append(encim.__name__)

            restrikcija_z_uporabnimi_encimi()
            
            
            for combobox in combobox_lista:
                for combobox in combobox_lista[1:]:
                    combobox.destroy()
                
                del(combobox_lista[1:])
            
            combobox_lista[0].state(["disabled"])
            combobox_lista[0].set("")
            gumb_shrani_restrikcijske_encime_kot.config(state= "disabled")
            lista_encimov_upload_gumb.config(state= "disabled")
            
            izbran_sestavljen_konstrukt = vektor_slika_checkbox.get()

 
        elif  checkbox_check.get() == 0:
            combobox_lista[0].state(["!disabled"])
            gumb_shrani_restrikcijske_encime_kot.config(state= "normal")
            lista_encimov_upload_gumb.config(state= "normal")
            
            Lista_encimov_na_voljo = [""]
            restrikcija_z_uporabnimi_encimi()
            
            izbran_sestavljen_konstrukt = vektor_slika_checkbox.get()
            
        if not izbran_sestavljen_konstrukt == "izberi si konstrukt, ki ga želiš analizirati":
            vektor_slika_checkbox.event_generate("<<ComboboxSelected>>")
            

    checkbox_check.trace_add("write", vsi_encimi_uporaba_checkbox)
    
    #%% sestavljanje sekvence
    STOP_kodoni = [Seq("TGA"),Seq("TAA"),Seq("TAG")]
    START_kodon= Seq("ATG")

    zaporedje = globals()["skupine_popolni_record"]["zaporedje"]
    
    začasna_lista1 = {}
    
    for record in globals()["skupine_popolni_record"][zaporedje[0]]:
        začasna_lista1[record.id] = [record.seq[:-3],[record.seq[:-3]]]
    
    def sekvenca_assembly(dic_s_podatki):
        dic_za_assembly = {}
        #print(skupina)
        #print(dic_s_podatki)
        #print(globals()["skupine_popolni_record"])
        for key, value in dic_s_podatki.items():
            
            for record in globals()["skupine_popolni_record"][skupina]:
                
                lista_sekvenc = value[1].copy()
                dodana_sekvenca = record.seq

                while dodana_sekvenca[:3] in START_kodon:
                    
                    dodana_sekvenca = dodana_sekvenca[3:]

                if not skupina in zaporedje[-1]:

                    while dodana_sekvenca[-3:] in STOP_kodoni:
                        dodana_sekvenca = dodana_sekvenca[:-3]
                
                lista_sekvenc.append(dodana_sekvenca)
                
                dic_za_assembly[key +"____" + record.id] = [value[0] + dodana_sekvenca,lista_sekvenc]
                
        return dic_za_assembly
    
    
    for skupina in zaporedje[1:]:
        začasna_lista1 = sekvenca_assembly(začasna_lista1)

    globals()["lista_sestavljenih_sekvenc"] = začasna_lista1

    imena_sestavljenih_konstruktov = [x for x in začasna_lista1.keys()]
    
    vektor_slika_checkbox_frame = tk.Frame(Vektor_restrikcija_toplevel,height= velikost_okna_sestavljanje_konstrukta/10)
    vektor_slika_checkbox_frame.pack(side="top",fill="both", expand=True, padx= 5, pady= 5)
    
    vektor_slika_checkbox = ""
    combobox_konstrukt_tekst = tk.StringVar(value= "izberi si konstrukt, ki ga želiš analizirati")
    def trace_podatki_restrikcijski_encimi_pridobljeni(*args):
        if encimi_pridobljeni.get():
            nonlocal vektor_slika_checkbox, combobox_konstrukt_tekst
            
            vektor_slika_checkbox=  ttk.Combobox(vektor_slika_checkbox_frame, values= imena_sestavljenih_konstruktov, textvariable= combobox_konstrukt_tekst)
            vektor_slika_checkbox.place(x=0, y= 0, relwidth=1, relheight=0.3)
            
            vektor_slika_checkbox.bind("<<ComboboxSelected>>",slika_checkbox_komanda)

            encimi_pridobljeni.trace_remove('write', vektor_slika_checkbox_tracing)
    
    vektor_slika_checkbox_tracing = encimi_pridobljeni.trace_add('write', trace_podatki_restrikcijski_encimi_pridobljeni)
        
    
    for key in globals()["lista_sestavljenih_sekvenc"].keys():
        
        locals()["label_check__" + key] = tk.Label(restrikcijski_encimi_določeni_check_frame_podframe_oznake, text= key, font=("Arial", 5))
        locals()["label_check__" + key].pack(side = "top", fill="both",expand=True, padx= 5, pady= 5)
        
        lista_checkov_za_restrikcijska_mesta.append(locals()["label_check__" + key])
        
    
    
                  
#%% preostanek komand
def podatki_pridobljeni(*args):
    
    global pridobljeni_podatki_flag

    if pridobljeni_podatki_flag.get():

            obdelava_podatkov = tk.Menubutton(orodna_vrstica, text= "obdelava_podatkov", relief=tk.RAISED)
            obdelava_podatkov_meni = tk.Menu(obdelava_podatkov, tearoff= 0)
            obdelava_podatkov.config(menu=obdelava_podatkov_meni)
    
            obdelava_podatkov.pack(side= "left", fill="both",expand=True, padx= 5, pady= 5)
            podatki_meni.add_command(label="uredi podatke", command=lambda reset = False: sestavljanje_konstrukta_toplevel(reset,globals()["skupine_popolni_record"]))
            

            obdelava_podatkov_meni.add_command(label= "Overlap extension", command= vektor_restrikcija_komanda)
            
            pridobljeni_podatki_flag.trace_remove('write', pridobljeni_podatki_tracing)

pridobljeni_podatki_tracing = pridobljeni_podatki_flag.trace_add('write', podatki_pridobljeni)



koncentracije_primerjev_v_M = 0.5*(1e-6)
zeljen_delez = 0.05
    
zelen_delez_homodimerov = koncentracije_primerjev_v_M * zeljen_delez
delez_nehomodimeriziranih = koncentracije_primerjev_v_M - koncentracije_primerjev_v_M*(1-zeljen_delez)   
razmerje_za_homodimere = zelen_delez_homodimerov /(pow(delez_nehomodimeriziranih,2))

hp_zelen_lnKxR = -8.314*math.log(zeljen_delez/(1-zeljen_delez))
homodimer_zelen_lnKxR = -8.314 * math.log(razmerje_za_homodimere)


def iterator_za_konstrukcijo_optimalnih_podaljškov(podaljšan_primer, prvotni_primer, tm):
    
    hp_zelen_dG = hp_zelen_lnKxR * tm
    homodimer_zelen_dG = homodimer_zelen_lnKxR*tm
    
    if len(podaljšan_primer)>=60:
        
        primer_iterator = str(podaljšan_primer[-59:])
        #print(primer_iterator)
    else:
        primer_iterator = str(podaljšan_primer)
    
    podaljšan_primer_hp_dG = -10e9
    podaljšan_primer_homo_dG = -10e9
    
    # print("PRIMER")
    while True:
        
        primer_iterator = primer_iterator[1:]
        
        podaljšan_primer_hp_TD = calc_hairpin(primer_iterator)
        podaljšan_primer_hp_dG = (podaljšan_primer_hp_TD.dh-podaljšan_primer_hp_TD.ds*tm)* 4.184
        
        podaljšan_primer_homo_TD = calc_homodimer(primer_iterator)
        podaljšan_primer_homo_dG = (podaljšan_primer_homo_TD.dh - podaljšan_primer_homo_TD.ds*tm)* 4.184

        if ((podaljšan_primer_hp_dG >= hp_zelen_dG) and (podaljšan_primer_homo_dG >= homodimer_zelen_dG)) or (len(primer_iterator) <= len(prvotni_primer)):
            return(primer_iterator)
            break

def konstrukcija_primerjev (ime_elementa,input_sekvenca):
    
    primerji = []
    
    f_pr = []
    r_pr = []
    f_pr_tm = []
    r_pr_tm = []
    PCR_tm = []
    
    dolžina_dela = len(input_sekvenca)
    
    seq_args = {
        'SEQUENCE_ID': str(ime_elementa),
        'SEQUENCE_TEMPLATE': str(input_sekvenca),
        
        'SEQUENCE_FORCE_RIGHT_START': [dolžina_dela-1],
        'SEQUENCE_FORCE_LEFT_START': [0],
    }
    
    primer_params = {        
    'PRIMER_OPT_SIZE': 18,
    'PRIMER_MIN_SIZE': 15,
    'PRIMER_MAX_SIZE': 20,
    'PRIMER_MAX_DIFF_TM': 5.0,
    
    
    'PRIMER_MAX_HAIRPIN_TH':100,
    
    'PRIMER_MIN_TM': 54.5,
    'PRIMER_OPT_TM': 60,
    'PRIMER_MAX_TM': 75,
    'PRIMER_NUM_RETURN': 10,
    'PRIMER_MIN_GC': 0,
    'PRIMER_MAX_GC': 100,                        
    'PRIMER_PRODUCT_SIZE_RANGE': [40,dolžina_dela],
    'PRIMER_EXPLAIN_FLAG': 1,
    
    'PRIMER_DNA_CONC': 500.0,
    'PRIMER_DNTP_CONC': 0.8
    }
                        
    result = bindings.design_primers(primer_params,seq_args)

    for i in range(9):
        lprimer_str = "PRIMER_LEFT_"+str(i)+"_SEQUENCE"
        rprimer_str = "PRIMER_RIGHT_"+str(i)+"_SEQUENCE"
        
        f_tm = float(result["PRIMER_LEFT_"+str(i)+"_TM"]) + 273
        r_tm = float(result["PRIMER_RIGHT_"+str(i)+"_TM"]) + 273
        
        primer_left_seq = Seq(result[lprimer_str])
        primer_right_seq = Seq(result[rprimer_str])
        
        f_pr.append(primer_left_seq)
        r_pr.append(primer_right_seq)
        f_pr_tm.append(f_tm)
        r_pr_tm.append(r_tm)
        PCR_tm.append(min(f_tm,r_tm)+3)

    primerji = {}
    primerji["f_pr"] = f_pr
    primerji["r_pr"] = r_pr
    primerji["f_tm"] = f_pr_tm
    primerji["r_tm"] = r_pr_tm
    primerji["PCR_tm"] = PCR_tm
        
    return(primerji)
    
def Sestavljanje_primerjev_komanda():
    
    Sestavljanje_oligov_lista = {}

    for key, value in globals()["lista_sestavljenih_sekvenc"].items():
        
        print(key)
        
        imena_elementov = globals()["skupine_popolni_record"]["zaporedje"].copy()
        
        relevantni_lista = value[1]
        ime_prvega_res_encima = value[-1][0]
        ime_zadnjega_res_encima = value[-1][1]
        
        Prvi_res_encim = Seq(ime_prvega_res_encima.site)
        Drugi_res_encim = Seq(ime_zadnjega_res_encima.site)
        
        relevantni_lista.insert(0,Prvi_res_encim)
        relevantni_lista.insert(len(relevantni_lista),Drugi_res_encim)
        
        relevantno_zaporedje_za_skupino = imena_elementov.copy()
                
        relevantno_zaporedje_za_skupino.insert(0,ime_prvega_res_encima)
        relevantno_zaporedje_za_skupino.insert(len(relevantno_zaporedje_za_skupino), ime_zadnjega_res_encima)
        
        sestavljen_del = Seq("")
        deli_z_primerji =[]
        
        primer_serija = {}
        imena_delov_s_primerji = []   
        
        i=0
        
        indeks_prejšnjega_dela_s_primerji = 0
        
        deli_s_primerji_in_sestavljeni = []
        
        #prva zanka za determinacijo delov, za katere bodo zasnovani oligonukleotidi
        for part in relevantni_lista:
            
            ime_relevantnega_dela = relevantno_zaporedje_za_skupino[i]
            i = i+1
            if len(part)>= 100:
                imena_delov_s_primerji.append(ime_relevantnega_dela)
                deli_z_primerji.append(part)

        i = 0
        deli_s_primerji_iterator = 0
        
        #druga zanka za determinacijo f_pr in r_pr
        for part in relevantni_lista:
            
            ime_relevantnega_dela = relevantno_zaporedje_za_skupino[i]
            
            if len(part) >= 100:
                #primer_bind_deli_design
                
                primer_bind_rezultati = konstrukcija_primerjev(ime_relevantnega_dela,part)
                primer_serija[ime_relevantnega_dela] = primer_bind_rezultati
                
                #f_pr design
               
                
                primer_serija[ime_relevantnega_dela]["f_pr_obdelan"] = []
                primer_serija[ime_relevantnega_dela]["r_pr_obdelan"] = []
                
                
                n = 0
                for f_pr in primer_serija[ime_relevantnega_dela]["f_pr"]:
                    tm = primer_bind_rezultati["PCR_tm"][n]
                    
                    if ime_relevantnega_dela == imena_delov_s_primerji[0]:
                         neobdelan_podaljšan_f_pr = sestavljen_del + f_pr 
                    else:
                        neobdelan_podaljšan_f_pr = relevantni_lista[indeks_prejšnjega_dela_s_primerji] + sestavljen_del + f_pr
                    
                    obdelan_primer_f = iterator_za_konstrukcijo_optimalnih_podaljškov(neobdelan_podaljšan_f_pr,f_pr,tm)
                    
                    primer_serija[ime_relevantnega_dela]["f_pr_obdelan"].append(obdelan_primer_f)
                    
                    n = n+1

                #r_pr razen_zadnji_design
                
                n = 0
                if not part in deli_z_primerji[0]:
                    
                    ime_prejsnjega_relavantnega_dela = imena_delov_s_primerji[deli_s_primerji_iterator-1]
                    
                    for r_pr in primer_serija[ime_prejsnjega_relavantnega_dela]["r_pr"]:

                        tm = primer_serija[ime_prejsnjega_relavantnega_dela]["PCR_tm"][n]
                        
                        neobdelan_podaljšan_r_pr = part.reverse_complement() + sestavljen_del.reverse_complement() + r_pr
                        obdelan_primer_r = iterator_za_konstrukcijo_optimalnih_podaljškov(neobdelan_podaljšan_r_pr,r_pr,tm)
                        primer_serija[ime_prejsnjega_relavantnega_dela]["r_pr_obdelan"].append(obdelan_primer_r)
                        
                        n = n+1
                        
                sestavljen_del = Seq("")
                
                deli_s_primerji_iterator = deli_s_primerji_iterator + 1
                indeks_prejšnjega_dela_s_primerji = i
   
            else:
                sestavljen_del = sestavljen_del + part            

            i= i+1
        
        #zadnji_rpr_design
        lista_elementov_vključenih_v_zadnji_podaljšek =relevantni_lista[indeks_prejšnjega_dela_s_primerji+1:]
        
        print("lista_zadnji_pod")
        print(lista_elementov_vključenih_v_zadnji_podaljšek)
        
        zadnji_podaljsek = Seq("").join(lista_elementov_vključenih_v_zadnji_podaljšek)
        
        print(zadnji_podaljsek)
        
        ime_zadnjega_dela_za_primanje = imena_delov_s_primerji[-1]
        podaljsan_obdelan_zadnji_primer_list = []
        
        n = 0 
        for zadnji_primer_bind_rpr in primer_serija[ime_zadnjega_dela_za_primanje]["r_pr"]:

            zadnji_primer_bind_tm = primer_serija[ime_zadnjega_dela_za_primanje]["PCR_tm"][n]            
            zadnji_podaljsan_primer = zadnji_podaljsek.reverse_complement() + zadnji_primer_bind_rpr
            
            obdelan_podaljsan_zadnji_rpr = iterator_za_konstrukcijo_optimalnih_podaljškov(zadnji_podaljsan_primer,zadnji_primer_bind_rpr,zadnji_primer_bind_tm)
            podaljsan_obdelan_zadnji_primer_list.append(obdelan_podaljsan_zadnji_rpr)
        
        primer_serija[ime_zadnjega_dela_za_primanje]["r_pr_obdelan"] = podaljsan_obdelan_zadnji_primer_list         
        Sestavljanje_oligov_lista[key] = primer_serija
        
        imena_elementov.append(key)     

        print()
    print(Sestavljanje_oligov_lista)
        
   
    
   
    #print("___")
    #print(globals()["lista_sestavljenih_sekvenc"])
    #print()
    #print(globals()["skupine_popolni_record"])
    
    
            
    
    
        
    
    #for key, value in relavantini_dic.items():
        
    #    for v in value[1]:
            


#     def dolocanje_minimalnih_dolzin_podaljskov(sekvenca,minimalna_tm):
        
#         kopija_sekvence = str(sekvenca)
        
#         while calc_tm(kopija_sekvence) > minimalna_tm:
#             print(calc_tm(kopija_sekvence))
#             kopija_sekvence = kopija_sekvence[1:]
            
#         return(Seq(kopija_sekvence))
            
        
    
#     def homo_hairpin_iterator(primer_prileganje_seq,tm,minimalni_primer,zacetni_podaljsek):
#         nonlocal modificirani_primerji
        
#         hp_dg_tm = 0     
        
#         count = 0
#         n = 0
        
#         obdelan_primer = str(zacetni_podaljsek+primer_prileganje_seq)
        
#         homodimer_termo = calc_homodimer(obdelan_primer)
#         homo_ds = homodimer_termo.ds
#         homo_dh = homodimer_termo.dh
        
#         homo_dg_tm = -(homo_dh - (homo_ds*tm))
        
#         primerjalna_homo_tm_dg = homodimer_K_R*tm
#         primerjalna_hp = hp_konstanta*tm 

#         print(obdelan_primer)
        
#         while (hp_dg_tm < primerjalna_hp or homo_dg_tm < primerjalna_homo_tm_dg) or not count == 0:

#             count = count + 1
            
#             if len(obdelan_primer) < len(minimalni_primer) and n < 6:
#                 count = 0
#                 n = n+1
#                 neurejen_primer = zacetni_podaljsek[n:] + Seq(n*"A") + primer_left_seq
            
#             obdelan_primer = str(neurejen_primer[count:])

#             hairpin_termo = calc_hairpin(obdelan_primer)
#             hp_ds = hairpin_termo.ds * 4.184 
#             hp_dh = hairpin_termo.dh * 4.184 
            
#             hp_dg_tm = -(hp_dh - hp_ds*tm)
            
#             homodimer_termo = calc_homodimer(obdelan_primer)
            
#             homo_ds = homodimer_termo.ds* 4.184 
#             homo_dh = homodimer_termo.dh* 4.184 
#             homo_dg_tm = (homo_dh - (homo_ds*tm))
            
            
#         return(obdelan_primer)


    
#     Sekvenca_vektor = globals()["skupine_popolni_record"]["vektor"][2].seq
#     koncentracije_primerjev_v_M = 0.5*(1e-6)
#     zeljen_delez = 0.01
#     hp_konstanta = 8.314*math.log(zeljen_delez/(1-zeljen_delez))
    
#     zelen_delez_homodimerov = koncentracije_primerjev_v_M * zeljen_delez
#     delez_nehomodimeriziranih = koncentracije_primerjev_v_M - koncentracije_primerjev_v_M*(1-zeljen_delez)
    
#     razmerje_za_homodimere = zelen_delez_homodimerov /(pow(delez_nehomodimeriziranih,2))
#     homodimer_K_R = -8.314 * math.log(razmerje_za_homodimere)


#     for key, value in globals()["lista_sestavljenih_sekvenc"].items():
        
#         counter = 0
#         dejanske_lokacije_rezov = ["",""]
        
#         for dejanska_encima in value[4]:
            
#             index_lokacije = value[2].index(dejanska_encima)
#             lokacije_enega_encima = value[3][index_lokacije]
            
#             dejanske_lokacije_rezov[counter] = lokacije_enega_encima[counter]
#             counter = counter - 1
                

#         mesto_5_5_podaljsek = value[4][0].fst3
#         mesto_5_3_podaljsek = value[4][0].fst5
        
        
#         if not mesto_5_5_podaljsek+ mesto_5_3_podaljsek == 0:
#             konec_mesto_5 = dejanske_lokacije_rezov[0] + abs(mesto_5_5_podaljsek+ mesto_5_3_podaljsek)-1
#         else:
#             konec_mesto_5 = dejanske_lokacije_rezov[0] + max(abs(mesto_5_5_podaljsek), abs(mesto_5_3_podaljsek))-1
        
#         začetek_mesto_5 = konec_mesto_5 - 55

#         mesto_3_5_podaljsek = value[4][1].fst3
#         mesto_3_3_podaljsek = value[4][1].fst5
        
#         if not mesto_3_3_podaljsek + mesto_3_5_podaljsek == 0:
#             začetek_mesto_3 = dejanske_lokacije_rezov[1] - abs(mesto_3_3_podaljsek + mesto_3_5_podaljsek)-1 
#         else:
#             začetek_mesto_3 = dejanske_lokacije_rezov[1] -1 - max(abs(mesto_3_3_podaljsek),abs(mesto_3_5_podaljsek))
            
            
#         konec_mesto_3 = začetek_mesto_3 + 55
        
#         mesto_5_podaljsek_neurejen = Sekvenca_vektor[začetek_mesto_5:konec_mesto_5]
#         mesto_3_podaljsek_neurejen = Sekvenca_vektor[začetek_mesto_3:konec_mesto_3]
        
#         minimalni_zahtevani_tm = 60
        
#         minimalna_mesto_5_podaljsek_neurejen_seq = mesto_5_podaljsek_neurejen
#         minimalna_mesto_3_podaljsek_neurejen_seq = mesto_3_podaljsek_neurejen

#         minimalna_mesto_5_podaljsek_neurejen_seq = dolocanje_minimalnih_dolzin_podaljskov(minimalna_mesto_5_podaljsek_neurejen_seq,minimalni_zahtevani_tm)
#         minimalna_mesto_3_podaljsek_neurejen_seq = dolocanje_minimalnih_dolzin_podaljskov(minimalna_mesto_3_podaljsek_neurejen_seq.reverse_complement(),minimalni_zahtevani_tm)
        
#         #n= 0
#         #while calc_tm(str(minimalna_mesto_5_podaljsek_neurejen_seq)) > minimalni_zahtevani_tm:
#         #    n= n+1
#         #    minimalna_mesto_5_podaljsek_neurejen_seq = minimalna_mesto_5_podaljsek_neurejen_seq[n:]
        
        
        
#         #n=-1
#         #while calc_tm(str(minimalna_mesto_3_podaljsek_neurejen_seq)) > minimalni_zahtevani_tm:
#         #    n= n-1
#         #    minimalna_mesto_3_podaljsek_neurejen_seq = minimalna_mesto_3_podaljsek_neurejen_seq[:n]
        
#         lista_primerjev = []
        
        
        
        
#         count = 0
#         for v in value[1]:
#             dolžina_dela = len(v)
            
#             if dolžina_dela >= 100:
                
#                 primerji = []
                
#                 seq_args = {
#                     'SEQUENCE_ID': "primer__" + key,
#                     'SEQUENCE_TEMPLATE': v,
                    
#                     'SEQUENCE_FORCE_RIGHT_START': [dolžina_dela-1],
#                     'SEQUENCE_FORCE_LEFT_START': [0],
#                 }
                
#                 primer_params = {        
#                 'PRIMER_OPT_SIZE': 18,
#                 'PRIMER_MIN_SIZE': 15,
#                 'PRIMER_MAX_SIZE': 20,
#                 'PRIMER_MAX_DIFF_TM': 5.0,
                
                
#                 'PRIMER_MAX_HAIRPIN_TH':100,
                
#                 'PRIMER_MIN_TM': 54.5,
#                 'PRIMER_OPT_TM': 60,
#                 'PRIMER_MAX_TM': 75,
#                 'PRIMER_NUM_RETURN': 10,
#                 'PRIMER_MIN_GC': 0,
#                 'PRIMER_MAX_GC': 100,                        
#                 'PRIMER_PRODUCT_SIZE_RANGE': [40,dolžina_dela],
#                 'PRIMER_EXPLAIN_FLAG': 1,
                
#                 'PRIMER_DNA_CONC': 500.0,
#                 'PRIMER_DNTP_CONC': 0.8
#                 }
                                    
#                 result = bindings.design_primers(primer_params,seq_args)
                
#                 #try:
#                 for i in range(9):
#                     lprimer_str = "PRIMER_LEFT_"+str(i)+"_SEQUENCE"
#                     rprimer_str = "PRIMER_RIGHT_"+str(i)+"_SEQUENCE"
                    
#                     l_tm = float(result["PRIMER_LEFT_"+str(i)+"_TM"]) + 273
#                     r_tm = float(result["PRIMER_RIGHT_"+str(i)+"_TM"]) + 273
                    
#                     primer_left_seq = Seq(result[lprimer_str])
#                     primer_right_seq = Seq(result[rprimer_str])
#                     primerji.append([[primer_left_seq,primer_right_seq], [l_tm,r_tm]])
                    
#                     modificirani_primerji = ["",""]
                    
                    
                    
#                     if count == 0:
                        
#                         zacetni_f_podaljsek = mesto_5_podaljsek_neurejen[len(mesto_5_podaljsek_neurejen)+ len(primer_left_seq)-59:]
#                         minimalni_f_primer = minimalna_mesto_5_podaljsek_neurejen_seq + primer_left_seq

#                         modificirani_primerji[0] = homo_hairpin_iterator(primer_left_seq,l_tm,minimalni_f_primer,zacetni_f_podaljsek)
                        
#                     elif count == len(globals()["skupine_popolni_record"]["zaporedje"])-1:
                        
#                         zacetni_r_podaljsek = mesto_3_podaljsek_neurejen.reverse_complement()[len(mesto_3_podaljsek_neurejen)+len(primer_right_seq)-59:]
#                         minimalni_r_primer = minimalna_mesto_3_podaljsek_neurejen_seq.reverse_complement() + primer_right_seq
                        
#                         modificirani_primerji[1] = homo_hairpin_iterator(primer_right_seq,r_tm,minimalni_r_primer,zacetni_r_podaljsek)

                    
#                     else:
#                         #element = globals()["skupine_popolni_record"]["zaporedje"][count]
                        
#                         print("OSTALO")
#                         #kolikšen je zadosten overlap med tem delom in naslednjim
#                         n = 0
#                         prejšnji_del_za_f_primer = Seq("")
#                         while len(prejšnji_del_za_f_primer) <= 60:
#                             n = n+1
#                             prejšnji_del_za_f_primer = prejšnji_del_za_f_primer + globals()["lista_sestavljenih_sekvenc"][key][1][count-n]
                        
#                         nasledjni_del_za_r_primer = Seq("")
#                         n = 0
#                         while len(nasledjni_del_za_r_primer) <= 60:
#                             n = n+1
#                             nasledjni_del_za_r_primer = globals()["lista_sestavljenih_sekvenc"][key][1][count+n]
                        
#                         dejanski_del = globals()["lista_sestavljenih_sekvenc"][key][1][count]
#                         #print(prejšnji_del_za_f_primer)
#                         #print(nasledjni_del_za_r_primer)
#                         #print(dejanski_del)
                        
                        
#                         minimalni_f_podaljsek = prejšnji_del_za_f_primer
#                         minimalni_r_podaljsek = nasledjni_del_za_r_primer
                        
#                         print(primer_left_seq)
#                         print(primer_right_seq)
                        
                        
#                         print(":::")
                        
#                         dolocanje_minimalnih_dolzin_podaljskov
#                         while calc_tm(str(minimalni_f_podaljsek)) > minimalni_zahtevani_tm:
#                             minimalni_f_podaljsek = minimalni_f_podaljsek[1:]
#                             print(minimalni_f_podaljsek)
 
                        
#                         n=1   
#                         while calc_tm(str(minimalni_r_podaljsek)) > minimalni_zahtevani_tm:
#                             minimalni_r_podaljsek = minimalni_r_podaljsek[1:]
#                             print(minimalni_r_podaljsek)

                            
                        
#                         #def homo_hairpin_iterator(primer_prileganje_seq,tm,minimalni_primer,zacetni_podaljsek):
                        
#                         print("#####")
                        
                        
                        
                        
#                     primerji.append(modificirani_primerji)
#                     print(modificirani_primerji)
     
#                 # except:
#                 #     print(".................................")
#                 #     print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
#                 #     print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
#                 #     print(".................................")
                
                
                        
#                 lista_primerjev.append(primerji)
#                 lista_primerjev.append(modificirani_primerji)
                
#                 print(lista_primerjev)
 
#             count = count + 1
            
            
#         try:
#             globals()["lista_sestavljenih_sekvenc"][key][5]= lista_primerjev
#             #globals()["lista_sestavljenih_sekvenc"][key][6]= lista_modificiranih_primerjev
#         except:
#             globals()["lista_sestavljenih_sekvenc"][key].append(lista_primerjev)
#             #globals()["lista_sestavljenih_sekvenc"][key].append(lista_modificiranih_primerjev)
#     print("#####")    
#     print(globals()["lista_sestavljenih_sekvenc"])
#     print()
    

globals()["Primer_design_gumb"] = ""

def restrikcija_končana(*args):
    
    if globals()["obdelava_z_restrikcijskimi_zakljucena"].get():
        
        if not globals()["Primer_design_gumb"] == "":
            globals()["Primer_design_gumb"].destroy()
        
        globals()["Primer_design_gumb"] = ttk.Button(orodna_vrstica,text= "Sestavi začetne oligoukleotide", command= Sestavljanje_primerjev_komanda)
        globals()["Primer_design_gumb"].pack(side= "left", fill="both",expand=True, padx= 5, pady= 5)

        

restrikcija_končana_tracing = globals()["obdelava_z_restrikcijskimi_zakljucena"].trace_add('write', restrikcija_končana)
        

root.mainloop()  
    
