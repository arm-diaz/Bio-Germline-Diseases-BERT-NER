from transformers import pipeline
from transformers import AutoTokenizer, AutoModelForTokenClassification
import pandas as pd
from tqdm import tqdm


def clean_output(outputs, tokenizer, pubmed_id, sentence_num):
    results = []
    current = []
    last_idx = 0
    # make to sub group by position
    for output in outputs:
        if output["index"]-1==last_idx:
            current.append(output)
        else:
            results.append(current)
            current = [output, ]
        last_idx = output["index"]
    if len(current)>0:
        results.append(current)
    
    # from tokens to string
    strings = []
    for c in results:
        tokens = []
        starts = []
        ends = []
        for o in c:
            tokens.append(o['word'])
            starts.append(o['start'])
            ends.append(o['end'])
        new_str = tokenizer.convert_tokens_to_string(tokens)
        if new_str!='':
            strings.append(dict(
                pubmed_id=pubmed_id,
                sentence_num=sentence_num,
                word=new_str,
                start = min(starts),
                end = max(ends),
                entity = c[0]['entity']
            ))
    return strings

def entity_table(pipeline, tokenizer, pubmed_id, sentence_num, **pipeline_kw):
    def create_table(text):
        return clean_output(
                pipeline(text, **pipeline_kw), tokenizer=tokenizer, pubmed_id=pubmed_id, sentence_num=sentence_num
            )
    return create_table

def save_entity_table(pipeline, tokenizer, data, file_name):
    entities = []
    for _, row in tqdm(data.iterrows(), total=data.shape[0]):
        sentence = row.sentence
        pubmed_id = row.pubmed_id
        sentence_num = row["sentence#"]
        entities.extend(entity_table(pipeline=pipeline, tokenizer=tokenizer, pubmed_id=pubmed_id, sentence_num=sentence_num)(sentence))
    entity_df = pd.DataFrame(entities)
    entity_df.to_csv(f"results/{file_name}.csv", index=None)


# Read file
data = pd.read_csv("results/pubmed_sentences.csv")

# Diseases: Named entity recognition pipeline, passing in a specific model and tokenizer
tokenizer_disease = AutoTokenizer.from_pretrained("drAbreu/bioBERT-NER-NCBI_disease", model_max_length=512) # dmis-lab/biobert-v1.1, drAbreu/bioBERT-NER-NCBI_disease", drAbreu/bioBERT-NER-BC2GM_corpus
model_disease = AutoModelForTokenClassification.from_pretrained("drAbreu/bioBERT-NER-NCBI_disease") 
ner_disease = pipeline(task="token-classification", model=model_disease, tokenizer=tokenizer_disease) # pass device=0 if using gpu

# Genes: Named entity recognition pipeline, passing in a specific model and tokenizer
model_gene = AutoModelForTokenClassification.from_pretrained("drAbreu/bioBERT-NER-BC2GM_corpus")
tokenizer_gene = AutoTokenizer.from_pretrained("drAbreu/bioBERT-NER-BC2GM_corpus", model_max_length=512)
ner_gene = pipeline('token-classification', model=model_gene, tokenizer=tokenizer_gene)

print("NER: Genes ...")
save_entity_table(ner_gene, tokenizer_gene, data, "genes")

print("NER: Diseases ...")
save_entity_table(ner_disease, tokenizer_disease, data, "diseases")

print("DONE!")