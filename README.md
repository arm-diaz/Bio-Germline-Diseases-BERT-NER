# Bio-Germline-Diseases-BERT-NER

## Workflow

![Figure: Workflow](images/workflow.png "Workflow")

## Dependencies

`Python 3.7`

## Setup

Install packages: `pip install  -r requirements.txt`

## Create Jupyter Notebook Kernel

`python -m ipykernel install --user --name=bio-germline-diseases-bert`

## Raw Data

Data can be found in the folder `data`.

Files:

+ `pubmed_abstracts`: CSV file contains Pubmed ID, and Abstract.
+ `gene_names`: CSV file contains list of genes and alternative gene names.
+ `gene_disease_assn`: CSV file contains list of genes, associated diseases, and alternative disease names.

## Run Pipeline

### Data Pre-processing

`nohup python pre-processing/tokenization.py &`

### Name Entity Recognition

`nohup python NER/ner_model.py &`

### Name Entity Disambiguation

`nohup python NER/ner_disambiguation.py &`

### Knowledge Graph

Execute jupyter notebook `knowledge-graph/create-graph`

You can also build the graph by running the cql query `knowledge-graph/crqCreateGraph.txt` in your Neo4j Browser.

### Charts

Execute jupyter notebook `charts/chart-wordcloud`

#### WordCloud Genes: Name Entity Recognition

![wordCloud1](images/wordCloudGene.png "wordCloud1")

#### WordCloud Genes: Name Entity Disambiguation

![wordCloud2](images/wordCloudDisambiguatedGene.png "wordCloud2")

#### WordCloud Diseases: Name Entity Recognition

![wordCloud3](images/wordCloudDisease.png "wordCloud3")

#### WordCloud Diseases: Name Entity Disambiguation

![wordCloud2](images/wordCloudDisambiguatedDisease.png "wordCloud2")


### General Statistics

Execute jupyter notebook `charts/stats`

## Results

Results can be found in the folder `results`.

Files:

+ `pubmed_sentences`: CSV file contains Pubmed ID, Sentence#, and Sentences from abstracts
+ `ner_genes`: CSV file contains Pubmed ID, Sentence#, Word (NER Gene), Start (Position), End (Position), Entity (Type)
+ `ner_diseases`: CSV file contains Pubmed ID, Sentence#, Word (NER Disease), Start (Position), End (Position), Entity (Type)
+ `ner_disambiguated_genes`: CSV file contains disambiguated entities for the `ner_genes` file
+ `ner_disambiguated_diseases`: CSV file contains disambiguated entities for the `ner_diseases` file
  
## Graphs

### Graph 1: Gene <- PubMedID -> Disease {name: 'Teeth (Benign)'}

![Graph1](images/graph1.png "Graph1")

### Graph 2: PubMedID -> Disease {name: 'Teeth (Benign)'}

![Graph2](images/graph2.png "Graph2")

### Graph 3: PubMedID -> Gene {name: 'AXIN2'}

![Graph3](images/graph3.png "Graph3")

### Graph 4: Gene <- PubMedID {name: '9024708'} -> Disease

![Graph4](images/graph4.png "Graph4")
