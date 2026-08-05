import os
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from primer3 import bindings

#wdir
wdir = r"C:\Users\HP\Desktop\MAGISTRSKA\sekevence\wdir"

#Primerji
Primerji = os.path.join(wdir,"Primerji")
primer_input = os.path.join(Primerji,"primer3input")
pr_out = os.path.join(Primerji, "primer3output")

#Biodeli
Biodeli = os.path.join(wdir,"Biodeli")
zapi = os.path.join(Biodeli,"zaporedje","zaporedje.txt")
terminator = os.path.join(Biodeli, "terminator")

programi_dir = r"C:\Users\HP\Desktop\working"

#------------------------------------------------------------------------------
def purge_folder(folder_path):
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isfile(item_path) or os.path.islink(item_path):
            os.remove(item_path)


for file in os.listdir(Primerji):
    if file in ("primer3output"):
        purge_folder(os.path.join(Primerji, file))
          
#--------------------------primer design---------------------------------------
zaporedje = []
with open(os.path.join(zapi)) as rd:
    for line in rd:
        a = line.strip()
        zaporedje.append(a)
print(zaporedje)

def pri (inpt, konc):
    num = 1
    print(inpt)
    loos = [] 
    with open(inpt,"r") as rd:
        for line in rd:
            for record in SeqIO.parse(rd,"fasta"):
                
                if "|" in record.id:
                    abb = record.id
                    bbe = abb.replace("|","-")
                    vak =""
                else:
                    abb = record.id
                    bbe = konc
                    vak = record.id
                print(bbe)
                with open(os.path.join(pr_out,bbe +"__" + vak + ".fasta"),"a") as primer:
                    
                    sequence = record.seq                    
                    target_length = len(sequence)
                    primer_params = {
                    'PRIMER_OPT_SIZE': 18,
                    'PRIMER_MIN_SIZE': 17,
                    'PRIMER_MAX_SIZE': 23,
                    'PRIMER_MIN_TM': 52,
                    'PRIMER_OPT_TM': 62,
                    'PRIMER_MAX_TM': 75,
                    'PRIMER_NUM_RETURN': num,
                    'PRIMER_MIN_GC': 0,
                    'PRIMER_MAX_GC': 100,                        
                    'PRIMER_PRODUCT_SIZE_RANGE': [target_length, target_length],
                    'PRIMER_PAIR_MAX_DIFF_SIZE' : 10,
                    'PRIMER_PAIR_MAX_DIFF_TM' : 3.0,
                    'PRIMER_EXPLAIN_FLAG': 1,
                    
                    'PRIMER_DNA_CONC': 500.0,
                    'PRIMER_DNTP_CONC': 0.8,
                    
                    'PRIMER_MAX_HAIRPIN_TH': 66,
                    }
                   
                    seq_args = {
                        'SEQUENCE_ID': abb,
                        'SEQUENCE_TEMPLATE': sequence,
                    }
                    
                    result = bindings.design_primers(primer_params,seq_args)
                    try:
                        print(result["PRIMER_LEFT_0_SEQUENCE"])
                        
                    except:
                        print(".................................")
                        print(record.id)
                        print("PRIMER_LEFT_EXPLAIN:", result.get("PRIMER_LEFT_EXPLAIN", "No explanation available"))
                        print("PRIMER_RIGHT_EXPLAIN:", result.get("PRIMER_RIGHT_EXPLAIN", "No explanation available"))
                        loos.append(record.id)
                        print(".................................")
                        continue
                    
                    primer_results = bindings.design_primers(primer_params,seq_args)
                    
                    for i in range(num):
                        
                        a = Seq(primer_results[f'PRIMER_LEFT_{i}_SEQUENCE']) 
                        b = Seq(primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE'])
                        fTm = str(primer_results[f'PRIMER_LEFT_{i}_TM'])
                        rTm = str(primer_results[f'PRIMER_RIGHT_{i}_TM'])
                        
                        smak = str(record.id)
                        kak = smak.replace(">", "")
                        bak = kak.strip()
                        print(bak)
                        
                        primer.write(f">forward_primer{i}" +"__"+ fTm +"__" + bak +"__"+ "CDS" +"\n")
                        primer.write(str(a)+ "\n" )
                        primer.write(f">reverse_primer{i}" +"__"+ rTm + "__" + bak +"__"+ "CDS" + "\n")
                        primer.write(str(b)+ "\n" )


encm = os.path.join(primer_input,"primer3inpt.fasta")
pri(encm,"")

    
for file in os.listdir(terminator):
    if "terminator" in zaporedje:
        pri(os.path.join(terminator,file),"ter")
    else:
        break


subprocess.run(["python", os.path.join(programi_dir, "5_rez_pretvorba_vektor_fasta.py")])
      
                

                                