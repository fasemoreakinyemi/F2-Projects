# F2-Notebook
## Week 1
### Papers
- dRNA review http://www.sciencedirect.com/science/article/pii/S1369527414000800
- dRNA Seq of VC http://www.pnas.org/content/112/7/E766.abstract
- Pre print http://mbio.asm.org/content/8/3/e00438-17.full
### Biopython exercises
- Translation of DNA sequences
```
def amino_maker (file,file_type):
    try:
        from Bio import SeqIO
        my_seq = SeqIO.parse(file,file_type)
    except:
        print("If you have properly installed Biopython, ensure the file matches with the type ")
    amn_dic = dict()
    for seq_records in my_seq:
        if len(seq_records.seq) % 3 == 0:
            sid = seq_records.id
            amn = seq_records.seq.translate()
            amn_dic[sid] = amn
           
        else: 
            print(seq_records.id ,": This sequence is not a multiple of 3")
    return amn_dic   
   ```
 - Fasta Parser
   ```
   def fasta_parser2(file):
    file_open = open(file)
    file_read = file_open.read()
    sequence = []
    header = []
    line_holder = []
    gene_dic = dict()
    file_split = file_read.strip().split('>')
    for lines in file_split:
        if len(lines) == 0: continue
        lines = lines.split('\n')
        line_holder.append(lines)
    for items in line_holder:
        if len(items) == 0: continue
        header_ = items[0]
        header.append(items[0])
        sequence_ = ''.join(items[1:])
        sequence.append(sequence_)
        gene_dic[header_] = sequence_
        for values in gene_dic:
            print(values +'\n'+ gene_dic[values])
  ```
## Week 2
### Transcription Start Site annotation
- data : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2263196
- Genome : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/965/GCF_000026965.1_ASM2696v1/
- annogesic

