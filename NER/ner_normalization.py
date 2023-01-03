import pandas as pd
from tqdm import tqdm

# Read annotations
gene_metadata_df = pd.read_csv("data/gene_names.csv")
gene_disease_assn_df = pd.read_csv("data/gene_disease_assn.csv")

# Read BioBERT results
diseases_df = pd.read_csv("results/ner_diseases.csv")
genes_df = pd.read_csv("results/ner_genes.csv")

copy_genes_df = genes_df.copy()
copy_diseases_df = diseases_df.copy()

# NER Disambiguation: Genes
print("NER Disambiguation: Genes ...")
for idx, rows in tqdm(copy_genes_df.iterrows(), total=copy_genes_df.shape[0]):
    word = rows.word
    if "mutyh" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MUTYH"
    elif ("apc" in str(word).lower()) or ("adenomatous polyposis coli" in str(word).lower()):
        genes_df.loc[idx, "gene"] = "APC"   
    elif "pten" in str(word).lower():
        genes_df.loc[idx, "gene"] = "PTEN"   
    elif "atm" in str(word).lower():
        genes_df.loc[idx, "gene"] = "ATM"   
    elif "chek2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "CHEK2"   
    elif "brca1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "BRCA1"   
    elif "nf1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "NF1"   
    elif "rad50" in str(word).lower():
        genes_df.loc[idx, "gene"] = "RAD50"   
    elif "atr" in str(word).lower():
        genes_df.loc[idx, "gene"] = "ATR"
    elif "mlh1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MLH1"
    elif "msh2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MSH2"
    elif "msh6" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MSH6"
    elif "tp53" in str(word).lower():
        genes_df.loc[idx, "gene"] = "TP53"
    elif "abraxas1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "ABRAXAS1"
    elif "brca2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "BRCA2"
    elif "brip1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "BRIP1"
    elif "cdh1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "CDH1"
    elif "mre11" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MRE11"
    elif "palb2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "PALB2"
    elif "pms2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "PMS2"
    elif "rad51c" in str(word).lower():
        genes_df.loc[idx, "gene"] = "RAD51C"
    elif "rad51d" in str(word).lower():
        genes_df.loc[idx, "gene"] = "RAD51D"
    elif "nbn" in str(word).lower():
        genes_df.loc[idx, "gene"] = "NBN"
    elif ("cadherin" in str(word).lower()) or ("cdh1" in str(word).lower()):
        genes_df.loc[idx, "gene"] = "CDH1"
    elif ("mmr" in str(word).lower()) or ("mismatch" in str(word).lower()):
        genes_df.loc[idx, "gene"] = "Lynch Gene"
    elif "wt1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "WT1"
    elif "cdc20" in str(word).lower():
        genes_df.loc[idx, "gene"] = "CDC20"
    elif "pbcdc20" in str(word).lower():
        genes_df.loc[idx, "gene"] = "PBCDC20"
    elif "smarce1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "SMARCE1"
    elif "pik3ca" in str(word).lower():
        genes_df.loc[idx, "gene"] = "PIK3CA"
    elif "kras" in str(word).lower():
        genes_df.loc[idx, "gene"] = "KRAS"
    elif "p53" in str(word).lower():
        genes_df.loc[idx, "gene"] = "TP53"
    elif "smad4" in str(word).lower():
        genes_df.loc[idx, "gene"] = "SMAD4"
    elif "rb1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "RB1"
    elif ("p16ink4a" in str(word).lower()) or ("p14" in str(word).lower()):
        genes_df.loc[idx, "gene"] = "CDKN2A"
    elif "ar" == str(word).lower():
        genes_df.loc[idx, "gene"] = "AR"
    elif "vhl" in str(word).lower():
        genes_df.loc[idx, "gene"] = "VHL"
    elif "arf" == str(word).lower():
        genes_df.loc[idx, "gene"] = "ARF"
    elif "ataxin" in str(word).lower():
        genes_df.loc[idx, "gene"] = "ATXN1"
    elif "ctbp2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "CTBP2"
    elif "pot1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "POT1"
    elif "txnip" in str(word).lower():
        genes_df.loc[idx, "gene"] = "TXNIP"
    elif ("ccnd1" in str(word).lower()) or ("cyclin d1" in str(word).lower()):
        genes_df.loc[idx, "gene"] = "CCND1"
    elif "muc6" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MUC6"
    elif "muc2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MUC2"
    elif "at1r" in str(word).lower():
        genes_df.loc[idx, "gene"] = "AGTR1"
    elif "at2r" in str(word).lower():
        genes_df.loc[idx, "gene"] = "AGTR2"
    elif "irf1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "IRF1"
    elif "igf2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "IGF2"
    elif "hogg1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "OGG1"
    elif "msr1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MSR1"
    elif "crlf2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "CRLF2"
    elif "pher2tyr" in str(word).lower():
        genes_df.loc[idx, "gene"] = "PTEN"
    elif "egfr" in str(word).lower():
        genes_df.loc[idx, "gene"] = "EGFR"
    elif "her2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "HER2"
    elif "her3" in str(word).lower():
        genes_df.loc[idx, "gene"] = "HER3"       
    elif "braf" in str(word).lower():
        genes_df.loc[idx, "gene"] = "BRAF"     
    elif "rnasel" in str(word).lower():
        genes_df.loc[idx, "gene"] = "RNASEL"     
    elif "gstp1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "GSTP1"     
    elif "nthl1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "NTHL1"     
    elif "neil1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "NEL1"     
    elif "neil2" in str(word).lower():
        genes_df.loc[idx, "gene"] = "NEL2"     
    elif "tdg" in str(word).lower():
        genes_df.loc[idx, "gene"] = "TDG"     
    elif "ung" in str(word).lower():
        genes_df.loc[idx, "gene"] = "UNG"     
    elif "mpg" in str(word).lower():
        genes_df.loc[idx, "gene"] = "MPG"     
    elif "wnt" in str(word).lower():
        genes_df.loc[idx, "gene"] = "WNT"     
    elif "galnt12" in str(word).lower():
        genes_df.loc[idx, "gene"] = "GALNT12"     
    elif "edar" in str(word).lower():
        genes_df.loc[idx, "gene"] = "EDAR"    
    elif "rnf43" in str(word).lower():
        genes_df.loc[idx, "gene"] = "RNF43"  
    elif ("patatin" in str(word).lower()) or ("pnpla3" in str(word).lower()):
        genes_df.loc[idx, "gene"] = "PNPLA3"   
    elif ("tm6sf2" in str(word).lower()) or ("transmembrane 6" in  str(word).lower()):
        genes_df.loc[idx, "gene"] = "TM6SF2"  
    elif "ndrg1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "NDRG1"  
    elif "kir" in str(word).lower():
        genes_df.loc[idx, "gene"] = "KIR"  
    elif "cd3" in str(word).lower():
        genes_df.loc[idx, "gene"] = "CD3"  
    elif "erg" in str(word).lower():
        genes_df.loc[idx, "gene"] = "ERG"  
    elif "bap1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "BAP1"  
    elif "smarca4" in str(word).lower():
        genes_df.loc[idx, "gene"] = "SMARCA4"  
    elif "atg16l1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "ATG16L1"  
    elif "pnpla3" in str(word).lower():
        genes_df.loc[idx, "gene"] = "PNPLA3"  
    elif "atg16l1" in str(word).lower():
        genes_df.loc[idx, "gene"] = "ATG16L1"  
    elif "plk4" in str(word).lower():
        genes_df.loc[idx, "gene"] = "PLK4"  
    elif "fa" in str(word).lower():
        genes_df.loc[idx, "gene"] = "FA"  
    elif "brca" in str(word).lower():
        genes_df.loc[idx, "gene"] = "BRCA1"  
    elif "akt" in str(word).lower():
        genes_df.loc[idx, "gene"] = "AKT"
    elif "aata" in str(word).lower():
        genes_df.loc[idx, "gene"] = "astA"  
    elif "asta" in str(word).lower():
        genes_df.loc[idx, "gene"] = "aatA"  
    elif ("emg" in str(word).lower()) or ("crc" in str(word).lower()) or ("t300a" in str(word).lower()) or \
        ("hpgds" in str(word).lower()) or ("glutathione" in str(word).lower()) or ("nan" in str(word).lower()) or \
        ("fap" in str(word).lower()) or ("ttr" in str(word).lower()) or ("dispersin" in str(word).lower()) or \
        ("hba1c" in str(word).lower()) or ("lynch" in str(word).lower()) or ("plasmid" in str(word).lower()) or \
        ("catenin" in str(word).lower()) or ("growth" in str(word).lower()) or \
        ("heterotrimeric" in str(word).lower()) or ("gnas1" in str(word).lower()) or \
        ("chd1" in str(word).lower()) or ("rarbeta" in str(word).lower()) or ("abl" in str(word).lower()) or \
        ("scf" in str(word).lower()) or ("human" in str(word).lower()) or ("factor" in str(word).lower()) or \
        ("g20210a" in str(word).lower()) or ("prostate" in str(word).lower()) or ("hfe" in str(word).lower()):
        pass
    elif str(word).lower() not in ["pca", "hgpin", "ho - 1", "se", "t", "at", "lp", "p", "associated", "ls", "alt",
                                   "bard1", "protein c", "c", "fibrinogen", "hemoglobin", "deletion", "l", "h",
                                   "gs - g", "gs - e", "x - linked g6pd g202a variant", "aj", "end", "nhl", "hbeag",
                                   "aj founder mutations", "g", "mipa - d159", "cb", "rarb", "zeb family members",
                                   "glucocorticoid - induced leucine zipper", "txnip", "i", "r", "ms", "t", "el", 
                                   "in", "gene", "as", "rh"]:
        gene_names = list(set([row.GeneMasterName for _, row in gene_metadata_df.iterrows() if (str(word).lower() in str(row.GeneMasterName).lower()) | (str(row.GeneMasterName).lower() in str(word).lower())]))
        if len(gene_names) != 0:
            genes_df.loc[idx, "gene"] = gene_names[0]
        elif len(gene_names) == 0:
            gene_names = list(set([row.GeneMasterName for _, row in gene_metadata_df.iterrows() if (str(word).lower() in str(row.AltNameGene).lower()) | (str(row.AltNameGene).lower() in str(word).lower())]))
            if len(gene_names) != 0:
                genes_df.loc[idx, "gene"] = gene_names[0]


genes_df.to_csv("results/ner_disambiguated_genes.csv", index=None)

# NER Disambiguation: Diseases

diseases_dict = {"hcc": "Hepatobiliary Cancer",
                "hboc": "Breast and Ovarian Cancer",
                "tsc": "Tuberous Sclerosis complex (TSC)",
                "fap": "Colorectal Neoplasia",
                "tc": "Testicular Cancer",
                "cm": "Cutaneous Melanoma",
                "md": "Menetrier's disease",
                "mtc": "Thyroid Cancer",
                "sls": "Lynch Syndrome",
                "ls": "Lynch Syndrome",
                 "aml": "Leukemia",
                 "bc": "Breast Cancer", 
                 "oc": "Ovarian Cancer",
                 "ibc": "Breast Cancer", 
                 "bcs": "Breast Cancer", 
                 "bcc": "Skin Cancer (Non-Melanoma)",
                 "bccs": "Skin Cancer (Non-Melanoma)",
                 "fpc":  "Pancreatic Cancer",
                 "pdac": "Pancreatic Cancer",
                 "gc": "Gastric Cancer",
                 "gc tumor": "Gastric Cancer",
                 "gcs": "Gastric Cancer",
                 "melanoma": "Melanoma",
                 " melanoma": "Melanoma",
                 "malignant melanoma": "Melanoma",
                 "hereditary melanoma": "Melanoma",
                 "in situ melanoma": "Melanoma",
                 "familial melanoma": "Melanoma",
                 "childhood melanoma":" Melanoma",
                 ", melanoma": "Melanoma",
                 "mnac": "Mesonephric Adenocarcinoma",
                 "hsa": "Blood (Benign)",
                 "as": "Blood (Benign)"
                 }

print("NER Disambiguation: Diseases ...")
for idx, rows in tqdm(copy_diseases_df.iterrows(), total=copy_diseases_df.shape[0]):
    word = rows.word
    if str(word).lower() in list(diseases_dict.keys()):
        diseases_df.loc[idx, "diseases"] = diseases_dict[str(word).lower()]

    elif "mpnsts" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Peripheral Nervous System Benign"  

    elif "down syndrome" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Down Dyndrome"  

    elif "mesonephric" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Mesonephric Adenocarcinoma"  

    elif "myoepithelioma" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Myoepithelioma"  

    elif "ataxia" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Ataxia Telangiectasia"  

    elif ("lynch" in str(word).lower()) or ("ls -" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Lynch Syndrome"  

    elif "hepatitis" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Hepatitis"  

    elif "hdgc" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Gastric Cancer"  

    elif "embryon" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Embryonal Cancer"  

    elif ("crc" in str(word).lower()) or ("colorectal" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Colorectal Cancer"

    elif ("teeth" in str(word).lower()) or ("dental" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Teeth (Benign)"

    elif ("testicular" in str(word).lower()) or ("testes" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Testicular Cancer"

    elif ("polyposis" in str(word).lower()) or ("polyps" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Colorectal Neoplasia"

    elif ("breast cancer") in str(word).lower() or ("mammary" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Breast Cancer"

    elif "ovarian cancer" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Ovarian Cancer"

    elif ("lung cancer") in str(word).lower() or ("laryn" in str(word).lower()) or ("pulmonary" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Lung Cancer"

    elif "cervical cancer" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Cervical Cancer"

    elif "kidney cancer" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Kidney (Benign)"

    elif "lymphoma" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Lymphoma"

    elif "liver cancer" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Liver (Benign)"

    elif "gastric cancer" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Gastric Cancer"

    elif ("colon cancer") in str(word).lower() or  ("rectum" in str(word).lower()) or \
         ("rectal" in str(word).lower()) or ("clorectal" in str(word).lower()) or ("colorectal" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Colorectal Cancer"

    elif ("gastrointestinal" in str(word).lower()) or ("digestive" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Gastrointestinal (Benign)"

    elif ("leukemia" in str(word).lower()) or ("leukaemia" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Leukemia"

    elif ("sarcoma" in str(word).lower()) and ("epithelioid" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Sarcoma"

    elif ("gallstone" in str(word).lower()) or ("hepatocellular" in str(word).lower()) or \
        ("biliary" in str(word).lower()) or ("bile duct" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Hepatobiliary Cancer"

    elif "stomach" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Gastrointestinal (Benign)"

    elif ("gi cancer" in str(word).lower()) or ("gi - " in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "GI Neoplasm"

    elif "renal" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Kidney Cancer"

    elif "breast" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Breast Cancer"

    elif "ovarian" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Ovarian Cancer"

    elif "pancreatic" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Pancreatic Cancer"
    
    elif "bladder" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Bladder Cancer"

    elif "pancreas" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Pancreas (Benign)"

    elif ("uterus" in str(word).lower()) or ("##ova" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Uterus Cancer"

    elif "sclerosis" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Multiple Sclerosis"

    elif ("cervix" in str(word).lower()) or ("cervical" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Cervical Cancer"

    elif "esophageal" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Esophageal Cancer"

    elif "prostate" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Prostate Cancer"

    elif ("muscle" in str(word).lower()) or ("myopathy" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Myopathy"

    elif "pituitary" in str(word).lower() or "brain" in str(word).lower() or "glioma" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Brain Cancer"

    elif "thyroid" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Thyroid Cancer"

    elif ("myelodysplasia" in str(word).lower()) or ("blood" in str(word).lower()) or \
         ("hemangiosarcoma" in str(word).lower()) or ("angiosarcoma" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Blood (Benign)"

    elif ("cortical" in str(word).lower()) or ("cortical dysplasia" in str(word).lower()) or ("cns" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "CNS (Benign)"

    elif ("osteofibrous" in str(word).lower()) or ("osteofibrous dysplasia" in str(word).lower()) or ("bone" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Bone (Benign)"

    elif "dysplasia" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = None

    elif "adrenal" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Adrenal Cortical Carcinoma"

    elif "vaginal" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Vaginal Cancer"

    elif "bone" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Bone (Benign)"

    elif ("brain" in str(word).lower()) or ("astrocytoma" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Brain Tumor"

    elif ("heart" in str(word).lower()) or ("cardiovascular" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Cardiovascular Neoplasm"

    elif ("eye" in str(word).lower()) or ("retinal" in str(word).lower()) or ("chrpe" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "CHRPE"

    elif "cirrhosis" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Cirrhosis"

    elif "cutaneous" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Cutaneous Melanoma"

    elif ("skin" in str(word).lower()) or ("basal" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Skin (Benign)"

    elif ("head and neck" in str(word).lower()) or ("squamous" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Head And Neck Cancer"

    elif "coffin" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Coffin-Siris syndrome"

    elif ("developmental" in str(word).lower()) or ("intellectual" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Developmental"

    elif ("autism" in str(word).lower()) or ("intellectual" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Psychiatric"

    elif "neuroendocrine" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Neuroendocrine"

    elif ("child" in str(word).lower()) or ("pediatric" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Childhood Cancer"

    elif "asept" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Aseptic Loosening"

    elif "duodenum" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Duodenum Cabcer"

    elif "dysmorphic" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Dysmorphic Feature"

    elif "ear" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Ear (Benign)"

    elif "endocrine" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Endocrine"

    elif "endometrial" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Endometrial Cancer"

    elif ("esophageal" in str(word).lower()) or ("esophagus" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Esophageal Cancer"

    elif "orofacial" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Orofacial Clefting"

    elif "facial palsy" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Facial Palsy"

    elif "facial papules" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Facial Papules"

    elif "fibroid" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Fibroid tumors"

    elif "fumarase" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Fumarase Deficiency"

    elif "genital" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Genitourinary (Benign)"

    elif ("small bowel" in str(word).lower()) or ("small - bowel" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Small Intestine Cancer"

    elif "hair" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Hair (Benign)"

    elif "hamartoma" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Hamartoma"

    elif "lipoma" in str(word).lower() or "lipid" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Lipoma"

    elif ("lymphatic" in str(word).lower()) or ("lymphangiomyomatosis" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Lymphangiomyomatosis"

    elif ("nodes" in str(word).lower()) or ("lymphoma" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Lymphoma"

    elif "metabolic" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Metabolic"

    elif "nail" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Nails"

    elif "urinary tract" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Urinary Tract Cancer"

    elif "ureter" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Ureter Cancer"

    elif "Uterus" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Uterus Cancer"

    elif ("nerve" in str(word).lower()) or ("paraganglioma" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Pheochromocytoma/Paraganglioma"

    elif ("wilms" in str(word).lower()) or ("nephroblastoma" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Wilms Tumor"

    elif ("mandibular" in str(word).lower()) or ("mdpl" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "MDPL syndrome"

    elif "oropharynx" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Oropharynx"

    elif "parathyroid" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Parathyroid"

    elif "prolactinoma" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Prolactinoma"

    elif ("pns" in str(word).lower()) or ("neuroblastoma" in str(word).lower()):
        diseases_df.loc[idx, "diseases"] = "Peripheral Nervous System Cancer"

    elif "neural tissue" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Peripheral Nervous System Neoplasm"

    elif "plasma" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Plasma Cell Cancer"

    elif "sebaceous" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Sebaceous Cancer"

    elif "soft tissue" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Soft Tissue Neoplasm"

    elif "spleen" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Spleen (Benign)"

    elif "tuberous sclerosis complex" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Tuberous sclerosis complex (TSC)"

    elif "rhabdoid" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Rhabdoid Tumour"

    elif "endometrioid" in str(word).lower():
        diseases_df.loc[idx, "diseases"] = "Endometrioid Carcinoma"

    elif str(word).lower() not in ["cancer", "tumor", "h", "cancers", "nan cancer", "malignancy", "carcinoma",
                                  "ps", "e.", "malignant tumours", "s", "sa", "and", "all", "p", "c", "b", "neo",
                                  "t cell malignancy", "fhit", "ap", "em", "adenosquamous carcinoma", "adenocarcinoma",
                                  "malignant tumors", "non - malignant disease", "non - cns solid tumors", "neop",
                                  "dysplasia", "malignancies", "second malignancies", "bilateral", "tumour", "asm",
                                  "adenoma", "moderate adenoma", "hpc", "non", "poor", "pm", "hh", "e", "afad", "os",
                                  "scc", "meta", "malignant neoplasms", "advanced cancers", "gallstone", "o", "ii",
                                  "ma", "cfs", "bi", "ngs tumor", "vhl", "bhd", "tsc", "amls", "ant", "fa", "ch", "dp",
                                  "chole", "-", "- induced malignancies", "g", "pbi", "the", "mps", "mp", "duodenal",
                                  "tumors", "hyperplasia", "mi", "ben", "benign", "hands and", "gene", "complex", "general",
                                  "dys", "im", "cag", "sg", "laryngeal",  "main cancers", "cll", "f", "ds", "rh d",
                                  "duo", "pro", "primary tumors", "u", "cvid", "or", "cy", "npc", "n", "fdrs", "ec",
                                  "es", "acc", "accs", "dcis", "hhc", "lc", "escc", "uc", "rp", "npcs",
                                  "eoc", "mm", "vus", "hl", "ad", "als", "ald", "at", "a", "di"]:
        diseases_names = list(set([row.DiseaseName for _, row in gene_disease_assn_df.iterrows() if (str(word).lower() in str(row.DiseaseName).lower()) | (str(row.DiseaseName).lower() in str(word).lower()) |
                                                                                                    (str(word).lower() in str(row.OldDiseaseName).lower()) | (str(row.OldDiseaseName).lower() in str(word).lower()) |
                                                                                                    (str(word).lower() in str(row.Z_Old_DiseaseName).lower()) | (str(row.Z_Old_DiseaseName).lower() in str(word).lower()) |
                                                                                                    (str(word).lower() in str(row.Organ).lower()) | (str(row.Organ).lower() in str(word).lower()) ]))
        if (len(diseases_names) != 0):
            diseases_df.loc[idx, "diseases"] = diseases_names[0]

diseases_df.to_csv("results/ner_disambiguated_diseases.csv", index=None)
