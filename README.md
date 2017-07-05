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
### TSS Prediction
- TSS Predator
- Tunning :  java -Xmx1G -jar
- MasterTable and TSS statistics

### TSS Prediction WITH Tss Predator and IGB
- installed IGB
- dataset = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1519465
- Genome file = ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Vibrio_cholerae/latest_assembly_versions/GCF_000006745.1_ASM674v1/
## Week 3
### ANNOgesic
- vibrio dataset
```
TEX_LIBS="/home/mandela/Downloads/vibrio/GSM1519456_D47E_S1_0.1_minus_TEX_div_by_10536538.0_multi_by_5191739.0_forward.wig:notex:1:a:+ \
/home/mandela/Downloads/vibrio/GSM1519456_D47E_S1_0.1_minus_TEX_div_by_10536538.0_multi_by_5191739.0_reverse.wig:notex:1:a:- \
/home/mandela/Downloads/vibrio/GSM1519457_D47E_S1_0.1_plus_TEX_div_by_8186187.0_multi_by_5191739.0_forward.wig:tex:1:a:+ \
/home/mandela/Downloads/vibrio/GSM1519457_D47E_S1_0.1_plus_TEX_div_by_8186187.0_multi_by_5191739.0_reverse.wig:tex:1:a:+ \
/home/mandela/Downloads/vibrio/GSM1519464_WT_S1_0.1_minus_TEX_div_by_11591377.0_multi_by_5191739.0_forward.wig:notex:2:a:+ \
/home/mandela/Downloads/vibrio/GSM1519464_WT_S1_0.1_minus_TEX_div_by_11591377.0_multi_by_5191739.0_reverse.wig:notex:2:a:+ \
/home/mandela/Downloads/vibrio/GSM1519465_WT_S1_0.1_plus_TEX_div_by_5191739.0_multi_by_5191739.0_forward.wig:tex:2:a:+ \
/home/mandela/Downloads/vibrio/GSM1519465_WT_S1_0.1_plus_TEX_div_by_5191739.0_multi_by_5191739.0_reverse.wig:tex:2:a:-"

annogesic tss_ps \
--tsspredator_path /home/mandela/Downloads/apps/TSSpredator-1.06.jar \
--fasta_files /home/mandela/Downloads/vibrio/GCF_000006745.1_ASM674v1_genomic.fna \
--annotation_files /home/mandela/Downloads/vibrio/GCF_000006745.1_ASM674v1_genomic.gff \
--tex_notex_libs $TEX_LIBS \
--condition_names nes_pred \
--validate_gene \
--replicate_tex all_1 \
--project_path ANNOgesic
```
## Week 4 TSS Comparison
### Datasets
- KAI (Vibrio) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1519465)
- Lars (Neisseria) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85252)
- E. coli (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78041)
- Campylobacter (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38883)
- Helicobacter (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67564)
- Upstream sequece extraction code
```
import pandas as pd
anno_file = pd.read_csv('/home/mandela/ANNOgesic/output/TSSs/MasterTables/MasterTable_NC_000915.1/MasterTable.tsv', sep = '\t', comment = '#')
k = open('heli.txt', 'w')
d_anno_file = anno_file[anno_file['detected'] == 1]
pos = d_anno_file['SuperPos']
seqe = d_anno_file['Sequence -50 nt upstream + TSS (51nt)']
gen = d_anno_file['Genome']
loc = d_anno_file['Locus_tag']
prod = d_anno_file['Product']
seq_lst = []
id_lst = []
gen_lst = []
loc_lst = []
prod_lst = []
for items in seqe:
    seq_lst.append(items)
for ele in pos:
    id_lst.append(ele)
for genes in gen:
    gen_lst.append(genes)
for products in prod:
    prod_lst.append(products)
for locus in loc:
    loc_lst.append(locus)
i = 0
while i <= (len(seq_lst)-1):
    ide = ('>' + str(id_lst[i]) + '|' + str(gen_lst[i]) + '|' + str(loc_lst[i]) + '|' + str(prod_lst[i]) + '|' + str(i))
    seq= str(seq_lst[i])
    k.write(ide + '\n')
    k.write(seq + '\n')
    i = i + 1
k.close()
```
- Wget Script for downloading genome files
```
source=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/905/GCF_000017905.1_ASM1790v1/
wget_dir=wget_folder
if [! -d wget_folder]; then
        mkdir -p $wget_dir
fi
rm -rf $wget_dir
wget -P $wget_dir ${source}/*fna.gz
wget -P $wget_dir ${source}/*gff.gz
gunzip $wget_dir/*gz

```
