# UniProt Scraper
Get UniProt summaries for genes appended to dataframe    

Example:   
```
UniProtID   ...
P#####      ...
P#####      ...
```
Will output:   
```
UniProtID   ...   Summary
P#####      ...   Info1
P#####      ...   Info2  
```

### Installation
Install dependencies:
```
pip install pandas numpy html2text
```

### Run
Navigate to `uniprot_scraper` script and execute the following:
```
python uniprot_scraper.py
```
Then follow the prompts...

### Notes
- If the column where the UniProt IDs are found has other characters besides the UniProt ID, the prompter will ask you what the characters before the ID are. 
Example:   
```
sp|P19283|gene_name
```
- For the example above you would provide `sp\|` as input. Note that special characters, such as `|` must be pre-pended with a `\`
