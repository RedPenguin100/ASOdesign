import pandas as pd
from Bio import SeqIO
import re


def get_expression_of_cell_line(cell_line, expression_file):
    # read header to get gene names
    with open(expression_file, 'r') as f:
        header = f.readline().strip().split(',')

    # find the cell line's row
    target_row = None
    with open(expression_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            if line.startswith(cell_line + ','):
                target_row = line.strip().split(',')
                break

    if target_row is None:
        raise ValueError(f"Cell line '{cell_line}' not found in file.")

    # build df from the row
    gene_names = header[1:]  # first column is the cell line name
    expression_values = target_row[1:]
    expression_series = pd.Series(expression_values, index=gene_names, name=cell_line).astype(float)
    expression_df = expression_series.reset_index()
    expression_df.columns = ['Gene', f'{cell_line}_expression_norm']

    # calculate TPM
    expression_df['expression_TPM'] = (2 ** expression_df[f'{cell_line}_expression_norm']) - 1

    # filter and sort
    expression_df = expression_df[expression_df['expression_TPM'] > 0]
    expression_df = expression_df.sort_values(by='expression_TPM', ascending=False).reset_index(drop=True)

    # save and return
    expression_df.to_csv(f'{cell_line}_expression.csv', index=False)

    return expression_df


def get_mutations_of_cell_line(cell_line, mutation_file):
    necessary_cols = ['ModelID', 'Ref', 'Alt', 'VariantType', 'VariantInfo', 'DNAChange', 'HugoSymbol']

    # read data
    mutation_data = pd.read_csv(mutation_file, usecols=necessary_cols)

    # filter by cell line
    cell_line_mutation = mutation_data[mutation_data['ModelID'] == cell_line]
    cell_line_mutation.to_csv(f'{cell_line}_mutations.csv', index=False)
    return cell_line_mutation


def get_record_by_hugo(indexed_fasta, hugo_symbol):
    for record in indexed_fasta:
        header = record.description
        if f"gene_symbol:{hugo_symbol}" in header:
            return record
    return None


def extract_cds_from_gtf(gtf_path):
    cds_rows = []

    with open(gtf_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # skip header
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            if feature != 'CDS':
                continue

            # extract transcript_id from attributes field
            attr_dict = {}
            for attr in attributes.strip().split(';'):
                if attr.strip():
                    key_value = attr.strip().split(' ', 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        attr_dict[key] = value.strip('"')

            transcript_id = attr_dict.get('transcript_id')
            if not transcript_id:
                continue

            cds_rows.append({
                'transcript_id': transcript_id,
                'start': int(start),
                'end': int(end),
                'strand': strand
            })

    return pd.DataFrame(cds_rows)


def prepare_gtf_dataframe(gtf_file):
    """
    Load and preprocess GTF file, keeping only 'exon' and 'CDS' entries with parsed transcript IDs.
    """
    col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=col_names)

    # Filter only 'exon' and 'CDS'
    df = df[df['feature'].isin(['exon', 'CDS'])]

    # Extract transcript_id from attributes
    df['transcript_id'] = df['attribute'].str.extract(r'transcript_id "([^"]+)"')

    # Drop entries without transcript_id (should be rare or malformed)
    df = df.dropna(subset=['transcript_id'])

    return df


def find_shift(tid_dict, mut_idx):
    mut_idx = int(mut_idx)
    if mut_idx <= 0:
        print(f"Warning: mut_idx must be 1-based (>= 1), got {mut_idx}")
        return None

    accu_len = 0
    s = 1
    n = -1
    cds_list = tid_dict['cds']

    while accu_len < mut_idx:
        n += 1
        if n >= len(cds_list):
            print(f"Warning: mut_idx {mut_idx} is out of bounds for total CDS length.")
            return None
        s, e = cds_list[n]
        accu_len += e - s + 1

    shift = s - tid_dict['start']

    return shift


def get_distance_to_start_codon_from_df(df, transcript_id):
    """
    Compute CDS start position (in transcript/mRNA coordinates) by flattening exons and locating CDS start.
    """
    # get exons and CDS for this transcript
    exons = df[(df['transcript_id'] == transcript_id) & (df['feature'] == 'exon')]
    cds = df[(df['transcript_id'] == transcript_id) & (df['feature'] == 'CDS')]

    if exons.empty or cds.empty:
        raise ValueError(f"Transcript {transcript_id} not found or missing exons/CDS in GTF.")

    strand = exons.iloc[0]['strand']

    # sort exons and CDS by genomic position (according to strand)
    if strand == '+':
        exons = exons.sort_values(by='start')
        cds_start = cds['start'].min()
    else:
        exons = exons.sort_values(by='start', ascending=False)
        cds_start = cds['end'].max()

    # compute CDS start position in transcript (mRNA) coordinates
    transcript_pos = 0
    for _, exon in exons.iterrows():
        exon_start = exon['start']
        exon_end = exon['end']
        exon_len = exon_end - exon_start + 1

        if strand == '+':
            if exon_start <= cds_start <= exon_end:
                return transcript_pos + (cds_start - exon_start)
        else:  # reverse strand
            if exon_start <= cds_start <= exon_end:
                return transcript_pos + (exon_end - cds_start)

        transcript_pos += exon_len

    raise ValueError(f"Could not find CDS start position in exons for transcript {transcript_id}")


# returns a dict that contains the orders for mutation

def mutation_dict(mut_row):
    # Skip if DNAChange is missing or empty
    if pd.isna(mut_row.get('DNAChange')) or mut_row['DNAChange'].strip() == '':
        return None

    try:
        # Split into ID and mutation string
        transcript_id, mutation = mut_row['DNAChange'].split(':')
    except ValueError:
        return None  # if it can't unpack into two parts

    if mut_row['VariantType'] == 'SNV':
        match = re.search(r'\d+', mutation)
        if match:
            start, end = match.span()
            location = mutation[start:end]
            from_nt, to_nt = mutation[end:].split('>')
            result = {'id': transcript_id, 'start': location, 'end': location, 'ref': from_nt, 'alt': to_nt, 'type': mut_row['VariantType']}
        else:
            result = {'id': transcript_id, 'start': None, 'end': None, 'ref': None, 'alt': None, 'type': None}

    elif mut_row['VariantType'] == 'deletion':
        numbers = re.findall(r'\d+', mutation)
        if len(numbers) == 1:
            result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[0], 'type': mut_row['VariantType']}
        elif len(numbers) >= 2:
            result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[1], 'type': mut_row['VariantType']}
        else:
            result = {'id': transcript_id, 'start': None, 'end': None, 'type': None}

    elif mut_row['VariantType'] == 'substitution':
        numbers = re.findall(r'\d+', mutation)
        if len(numbers) == 1:
            result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[0]}
        elif len(numbers) >= 2:
            result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[1]}
        else:
            result = {'id': transcript_id, 'start': None, 'end': None, 'type': None}

        if 'delins' in mutation:
            result['type'] = 'delins'
            result['subs'] = mutation.split('delins')[-1]
        elif 'inv' in mutation:
            result['type'] = 'inversion'

    elif mut_row['VariantType'] == 'insertion':
        numbers = re.findall(r'\d+', mutation)
        if len(numbers) == 1:
            result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[0]}
        elif len(numbers) >= 2:
            result = {'id': transcript_id, 'start': numbers[0], 'end': numbers[1]}
        else:
            result = {'id': transcript_id, 'start': None, 'end': None, 'type': None}

        if 'ins' in mutation:
            result['type'] = 'insertion'
            result['ins'] = mutation.split('ins')[-1]
        elif 'dup' in mutation:
            result['type'] = 'duplication'
    else:
        return None

    result['info'] = mut_row.get('VariantInfo')
    return result


def comp_strand(seq):
    seq = seq[::-1]
    new_seq = ''
    for n in range(0,len(seq)):
        if seq[n] == 'A' or seq[n] =='a':
            new_seq += 'T'
        elif seq[n] == 'T' or seq[n] =='t' :
            new_seq += 'A'
        elif seq[n] == 'C' or seq[n] == 'c':
            new_seq += 'G'
        elif seq[n] == 'G' or seq[n] == 'g':
            new_seq += 'C'
        else:
            new_seq += seq[n]
    return new_seq


def mutate(mut_dict, shift, sequence):
    if 'splice' in mut_dict['info'][:5]:
        return "Splice variant cannot be predicted"

    else:
        mut_seq = list(sequence)

        if mut_dict['type'] == 'SNV':

            loc = int(mut_dict['start']) + shift - 1
            ref = mut_dict['ref']
            alt = mut_dict['alt']

            if sequence[loc] == ref:
                mut_seq[loc] = alt

            else:
                mut_seq = 'mismatch'

        elif mut_dict['type'] == 'deletion':

            start = int(mut_dict['start']) + shift - 1
            end = int(mut_dict['end']) + shift
            del mut_seq[start:end]

        elif mut_dict['type'] == 'insertion':

            start = int(mut_dict['start']) + shift - 1
            mut_seq[start] += mut_dict['ins']

        elif mut_dict['type'] == 'duplication':

            start = int(mut_dict['start']) + shift - 1
            end = int(mut_dict['end']) + shift
            dup = mut_seq[start:end]
            mut_seq[end-1] += ''.join(dup)

        elif mut_dict['type'] == 'delins':

            start = int(mut_dict['start']) + shift - 1
            end = int(mut_dict['end']) + shift
            del mut_seq[start:end]
            mut_seq[start-1] += mut_dict['subs']

        elif mut_dict['type'] == 'inversion':

            start = int(mut_dict['start']) + shift - 1
            end = int(mut_dict['end']) + shift
            sub = ''.join(mut_seq[start:end])
            new = comp_strand(sub)
            del mut_seq[start:end]
            mut_seq[start - 1] += new

    return ''.join(mut_seq)


def final_func(seq_fasta, expression_file, mutation_file, cds_gtf_file, cell_line, top_n):
    # read the mutation and expression data
    mut_data = get_mutations_of_cell_line(cell_line, mutation_file)
    print('read mutation data')
    exp_data = get_expression_of_cell_line(cell_line, expression_file)
    exp_data = exp_data.head(top_n)
    print('read expression data')

    # add new columns
    exp_data['Transcript_ID'] = None
    exp_data['Mutated Transcript sequence'] = None
    exp_data['Original Transcript sequence'] = None

    records = list(SeqIO.parse(seq_fasta, "fasta"))
    record_by_hugo = {}
    for record in records:
        match = re.search(r'gene_symbol:([^\s]+)', record.description)
        if match:
            hugo = match.group(1)
            if hugo not in record_by_hugo:  # only take the first occurrence
                record_by_hugo[hugo] = record
    # add to exp_data the id
    for idx, row in exp_data.iterrows():
        hugo_symbol = row['Gene'].split()[0]

        record = record_by_hugo.get(hugo_symbol)
        if record:
            exp_data.at[idx, 'Transcript_ID'] = record.id
            exp_data.at[idx, 'Original Transcript sequence'] = str(record.seq)
        else:
            exp_data.at[idx, 'Transcript_ID'] = None
            exp_data.at[idx, 'Original Transcript sequence'] = None
        print('got sequence')
    print('finished getting sequences')

    ############################################################################################################

    # Parse GTF file
    cds_data = prepare_gtf_dataframe(cds_gtf_file)
    print('parsed GTF')

    # Convert exp_data to dict indexed by transcript ID for O(1) lookup
    exp_data_indexed = exp_data.set_index('Transcript_ID')

    for idx, row in mut_data.iterrows():
        mut_dict = mutation_dict(row)

        if mut_dict is not None:
            tid = mut_dict['id']

            if tid in exp_data_indexed.index:
                try:
                    shift = get_distance_to_start_codon_from_df(cds_data, tid)
                    print(f"{tid} shift: {shift}")
                    seq = exp_data_indexed.at[tid, 'Original Transcript sequence']
                    mutated_seq = mutate(mut_dict, shift, seq)
                    exp_data_indexed.at[tid, 'Mutated Transcript sequence'] = mutated_seq
                    print(f'mutated {tid}')
                except Exception as e:
                    print(f"Skipping {tid} due to error: {e}")

    ##############################################################################################################

    exp_data = exp_data_indexed.reset_index()
    exp_data.to_csv(f'{cell_line}_transcriptome_top{top_n}.csv', index=False)
    print('finished!')
    return exp_data


celline_list = ['ACH-001328',
                'ACH-000463',
                'ACH-001188',
                'ACH-001086',
                'ACH-000739',
                'ACH-000232',
                ]

# for line in celline_list:
#     final_func(
#         'Homo_sapiens.GRCh38.cdna.all.fa',
#         'OmicsExpressionProteinCodingGenesTPMLogp1.csv',
#         'OmicsSomaticMutations.csv',
#         'gencode.v48.chr_patch_hapl_scaff.annotation.gtf',
#         line, 500
#     )
