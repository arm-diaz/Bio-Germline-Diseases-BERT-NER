import pandas as pd
import nltk
nltk.download('punkt')
nltk.download('stopwords')
# nltk.download('corpus')

# Read abstracts
abstracts_df = pd.read_csv("data/pubmed_abstracts.csv")

# Breakdown abstract into sentences
sentences = []
count = 1
for idx, row in abstracts_df.iterrows():
    tokens = nltk.sent_tokenize(row.abstract_text)
    for sentence in tokens:
        sentence_num = f"sentence: {count}"
        pubmed_id = row.pubmed_id
        sentences.append([pubmed_id, sentence_num, sentence])
        count += 1

print(sentences[0:5])

# Convert into DataFrame
sentence_df = pd.DataFrame(sentences, columns=["pubmed_id", "sentence#", "sentence"])
sentence_df.to_csv("results/pubmed_sentences2.csv", index=None)
