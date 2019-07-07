"""
Name: UniProt Scraper
Author: Jordan A. Berg
Affiliation: Rutter Lab, University of Utah
Version: 0.0.1
Last Update: 7 Jul 2019

MIT License

Copyright (c) 2019 Jordan Berg

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import sys
import pandas as pd
pd.options.mode.chained_assignment = None
import csv
from math import ceil

from urllib.request import Request, urlopen
import html2text

import multiprocessing # For debugging
from multiprocessing import cpu_count, Pool
import gc

def get_file(file):

    if file.endswith('.csv'):
        print('Error: comma-separated tables not allowed. Please provide a tab-delimited table with the suffix .txt or .tsv')
        print('Exiting...')
        sys.exit(1)
    else:
        sep = '\t'

    df = pd.read_csv(file, sep=sep)

    return df, sep

def scrape(uniprot_id):

    url = 'https://www.uniprot.org/uniprot/' + str(uniprot_id)

    try:
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        webpage = urlopen(req).read()
        page = webpage.decode('utf-8').splitlines()

        for line in range(len(page) + 1):
            if '</html>' in page[line]:
                break

            if 'This section provides any useful information about the protein, mostly biological knowledge' in page[line] or 'This subsection of the ‘Pathology and Biotech’ section provides information' in page[line]:
                return html2text.html2text(page[line].split("Function<sup>i</sup>")[1])
    except:
        print('Could not find annotations for ' + str(uniprot_id))

def run_scraper(chunk):
    for index, row in chunk.iterrows():

        scapped = scrape(row[id_loc])
        chunk.at[index,'summary'] = scapped

    return chunk

def make_chunks(data):

    cpu = cpu_count()
    chunks = []
    chunk_size = int(ceil(len(data.index.tolist()) / cpu))
    chunk_counter = 0

    for x in range(cpu):

        if (chunk_counter + chunk_size - 1) > data.index[-1] or x == max(range(cpu)):
            new_chunk = data.iloc[chunk_counter:data.index[-1] + 1]
            chunks.append(new_chunk)
            return chunks

        else:
            new_chunk = data.iloc[chunk_counter:chunk_counter + chunk_size - 1]
            chunks.append(new_chunk)
            chunk_counter += chunk_size - 1

    return chunks

if __name__ == 'uniprot_scraper':

    # Get file
    file = input('Please provide the full path and filename of the file you would like to work with: ')
    column = input('What is the name of the column containing the UNIPROT IDs (case-sensitive)?: ')
    solo = input('Does this column contain any other annotations beside the UNIPROT ID? (yes/no): ').lower()

    data, sep = get_file(str(file))

    if solo == 'yes':
        remove = input('Please provide the text ahead of the UNIPROT ID (note, any special characters must be proceded by a \"\\\"): ')

        data[column] = data[column].str.replace(remove,'')
        data[column] = data[column].str[:6]

    # Get UNIPROT gene summary
    id_loc = data.columns.get_loc(column)
    data['summary'] = ''

    # Split dataframe up
    chunks = make_chunks(data)
    chunks = [x for x in chunks if x is not None]

    # Parallel process
    cores = len(chunks)
    pool = Pool(cores)
    chunks = pool.map(run_scraper, chunks)
    pool.close()
    pool.join()
    gc.collect()

    # Merge
    data = pd.concat(chunks)
    data.tail()

    data['summary'] = data['summary'].str.replace('\n', ' ')

    # Output
    data.to_csv(
        str(file)[:-4] + '_annotated' + str(file)[-4:],
        sep = sep,
        index = False)
