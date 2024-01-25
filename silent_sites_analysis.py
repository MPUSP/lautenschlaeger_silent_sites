import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import pysam
import subprocess
import datetime
import gzip
from multiprocessing import Process, Queue
import shutil

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-c", "--config", dest="config", help="Config file [None]", default=None)
parser.add_argument("-n", "--cores", dest="cores", type=int, help="Number of cores to use [1]", default=1)
parser.add_argument("-i", '--input', dest="input", help="Input table [None]", default=None)
parser.add_argument('-t', '--temp', dest="temp", help="Temporary folder [./temp]", default='./temp')
parser.add_argument('-o', '--output', dest="output", help="Output folder [./output]", default='./output')
parser.add_argument('-r', '--region_size', dest="region_size", type=int, help="Minimum region size [100]", default=100)
parser.add_argument('-p', '--promoter_region', dest="promoter_region", type=int, help="Promoter region [250]", default=250)
parser.add_argument('--anno', dest="anno", help="Annotation files and info [None]", default='')
parser.add_argument('--coverage_threshold', dest="coverage_threshold", type=float, help="Normalized coverage threshold [0.02]", default=0.02)
parser.add_argument('--verbose', dest='verbose', action='store_true', help="Verbose mode [False]", default=False)
parser.add_argument('-k', '--keep', dest='keep', action='store_true', help="Keep temporary files [False]", default=False)

args = parser.parse_args()

#classes
class defaults:
    cores = None
    temp_folder = ''
    output_folder = ''
    strains = set()
    min_region_size = None
    promoter_region = None
    coverage_threshold = None
    
    class mapping:
        bwa = ''
        samtools = ''

#functions
def subc(command):
    print(command, flush=True)
    ctr = True
    try:
        # Execute the command
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Wait for the command to complete
        stdout, stderr = process.communicate()

        # Decode the output and error (if any)
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')

        if process.returncode == 0:
            print("Output:\n", stdout)
        else:
            print("Error:\n", stderr)

    except Exception as e:
        print("An error occurred: ", str(e))
        ctr = False
    return ctr

def parse_input(input_file):
    dataset = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            if not line.startswith('#'):
                line = line.strip().split()#ref and fwd and rev if available
                for file in line:
                    if not os.path.exists(file):
                        print(F'File {file} not found. Aborting.')
                        sys.exit(1)
                if not line[0] in dataset: dataset[line[0]] = []
                dataset[line[0]].append(line[1:])
    return dataset

def read_gff(ref_gff):
    annotation = {}
    if ref_gff.endswith('.gz'):
        infile = gzip.open(ref_gff, 'rt')
    else:
        infile = open(ref_gff, 'r')

    for line in infile:
        if line[0]=='#':
            continue
        else:
            line=line.strip().split('\t')
            if line[2]=='gene':
                ID=(int(line[3])-1,int(line[4])-1,line[6])
                annotation[ID]=dict([v.strip().split('=',1) for v in line[8].split(';')])
    infile.close()
    return annotation

def parse_anno(anno_file):
    dataset = {}
    print(F"Parsing annotation file {anno_file} ...", flush=True)
    with open(anno_file, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            
            line[0] = defaults.temp_folder+'reference/'+line[0].split('/')[-1]
            if not line[0].endswith('gz'): line[0]+='.gz'
            
            if not line[0] in dataset: dataset[line[0]]={}
            if len(line)==2:
                if not os.path.exists(line[1]):
                    print(F"Annotation file {line[1]} not found. Aborting.")
                    sys.exit(1)
                dataset[line[0]]['anno']=read_gff(line[1])
            elif len(line)==3:
                if not line[1] in dataset[line[0]]: dataset[line[0]][line[1]] = []
                dataset[line[0]][line[1]].append([int(p) for p in line[2].split(':')])
    return dataset

def parse_config(config):
    
    def find_executable(name):
        try:
            subprocess.run(['whereis', name], stdout=subprocess.PIPE, text=True)
            paths = result.stdout.strip().split()
            if len(paths) > 1:
                return paths[1]  # The first element is the command itself, the second is the path
            else:
                return None
        except subprocess.SubprocessError as e:
            print(f"An error occurred: {e}")
            return None
    
    print(F"Parsing config file and validating arguments ...", flush=True)
    config_required=0
    if config=='':
        for prog in ['bwa', 'samtools']:
            res = find_executable(prog)
            if res:
                config_required+=1
                if prog=='bwa': defaults.mapping.bwa = res
                elif prog=='samtools': defaults.mapping.samtools = res
                print(F"Required program {prog} found at {res}.", flush=True)
        if config_required<2:
            print('No config file specified and programs missing. Aborting. Please read the README.')
            sys.exit(1)
    elif not os.path.exists(config):
        print(F'Config file {config} not found. Aborting.')
        sys.exit(1)
    else:
        with open(config, 'r') as infile:
            for line in infile:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    if line[0].lower()=='bwa':
                        if not os.path.exists(line[1]):
                            print(F'BWA executable {line[1]} not found. Aborting.')
                            sys.exit(1)
                        else:
                            defaults.mapping.bwa = line[1]
                            if args.verbose: print(F"Using BWA executable {line[1]}", flush=True)
                    elif line[0].lower()=='samtools':
                        if not os.path.exists(line[1]):
                            print(F'Samtools executable {line[1]} not found. Aborting.')
                            sys.exit(1)
                        else:
                            defaults.mapping.samtools = line[1]
                            if args.verbose: print(F"Using Samtools executable {line[1]}", flush=True)
                    else:
                        print(F"Unknown parameter {line[0]} in config file. Aborting.")
                        sys.exit(1)
    #check program parameters
    if args.cores<1: defaults.cores = 1
    else: defaults.cores = args.cores
    
    if args.region_size<1:
        print(F"Minimum region size must be >=1. Aborting.")
        sys.exit(1)
    else:
        defaults.min_region_size = args.region_size
    
    if args.promoter_region<0:
        print(F"Promoter region must be >=0. Aborting.")
        sys.exit(1)
    else:
        defaults.promoter_region = args.promoter_region
    
    if args.coverage_threshold<0:
        print(F"Coverage threshold must be >=0. Aborting.")
        sys.exit(1)
    else:
        defaults.coverage_threshold = args.coverage_threshold
    
    if args.input=='':
        print('No input table specified. Aborting.')
        sys.exit(1)
    elif not os.path.exists(args.input):
        print(F'Input table {args.input} not found. Aborting.')
        sys.exit(1)
    else:
        input_dataset = parse_input(args.input)
    
    if args.temp=='': defaults.temp_folder = './temp/'
    else: defaults.temp_folder = args.temp
    if not os.path.exists(defaults.temp_folder): os.makedirs(defaults.temp_folder)
    if not defaults.temp_folder[-1]=='/': defaults.temp_folder+='/'
    for subfolder in ['mappings', 'reference']:
        if not os.path.exists(defaults.temp_folder+subfolder): os.makedirs(defaults.temp_folder+subfolder)
    
    if args.output=='': defaults.output_folder = './output/'
    else: defaults.output_folder = args.output
    if not os.path.exists(defaults.output_folder): os.makedirs(defaults.output_folder)
    if not defaults.output_folder[-1]=='/': defaults.output_folder+='/'
    
    if args.anno=='':
        print('No annotation files specified. Aborting.')
        sys.exit(1)
    elif not os.path.exists(args.anno):
        print(F'Annotation file {args.anno} not found. Aborting.')
        sys.exit(1)
    else:
        anno = parse_anno(args.anno)
    
    return input_dataset, anno
    
def copying_ref_files(input_dataset):
    for ref in list(input_dataset.keys()):
        if not os.path.exists(defaults.temp_folder+'reference/'+ref.split('/')[-1]+'.fas.gz'):
            command = F'cp {ref} {defaults.temp_folder}reference/'
            res = subc(command)
            if not res: sys.exit(1)
        
        if not ref.endswith('.gz'):
            command = F"gzip {defaults.temp_folder}reference/{ref.split('/')[-1]}"
            res = subc(command)
            if not res: sys.exit(1)
            input_dataset[defaults.temp_folder+'reference/'+ref.split('/')[-1]+'.gz'] = input_dataset.pop(ref)
            ref = defaults.temp_folder+'reference/'+ref.split('/')[-1]+'.gz'
        else:
            input_dataset[defaults.temp_folder+'reference/'+ref.split('/')[-1]] = input_dataset.pop(ref)
            ref = defaults.temp_folder+'reference/'+ref.split('/')[-1]
        
        strain = ref.split('/')[-1].rsplit('.',2)[0]
        
        defaults.strains.add((strain, ref))
    
    return input_dataset

def prepare_references(input_dataset):
    for key in input_dataset:
        if not os.path.exists(F"{key}.amb"):
            command = f'{defaults.mapping.bwa} index {key}'
            if args.verbose: print(command, flush=True)
            res = subc(command)
            if not res: sys.exit(1)
    return True

def run_mapping(input_dataset):
    print('Starting mappings ...', flush=True)
    for strain, ref in defaults.strains:
        mapping_folder = defaults.temp_folder+'mappings/'+strain+'/'
        if not os.path.exists(mapping_folder):
            os.makedirs(mapping_folder)
        for read_file_list in input_dataset[ref]:
            if read_file_list[0].endswith('.gz'):
                ID = read_file_list[0].split('/')[-1].rsplit('.',2)[0]
            else:
                ID = read_file_list[0].split('/')[-1].rsplit('.',1)[0]
            print(F"Processing {strain} - {ID}", flush=True)
            
            if not os.path.exists(F"{mapping_folder}{ID}.dedup.bam"):
                fq_files = ' '.join(read_file_list)
                command = F'{defaults.mapping.bwa} mem -t {defaults.cores} -R "@RG\\tID:ID\\tLB:lib1\\tPL:illumina\\tSM:None" {ref} {fq_files}'
                command += F" | {defaults.mapping.samtools} view -b --threads {defaults.cores} -"
                command += F" | {defaults.mapping.samtools} sort -@ {defaults.cores} -T 'tmp' -"
                command += F" | {defaults.mapping.samtools} rmdup --reference {ref} - {mapping_folder}/{ID}.dedup.bam"
                command += F" && {defaults.mapping.samtools} index {mapping_folder}/{ID}.dedup.bam"
        
                res = subc(command)
                if not res: sys.exit(1)
    return True

def extract_coverage_to_file():
    cov_files = {}
    for strain, ref in defaults.strains:
        cov_files[strain]=[]
        mapping_folder = defaults.temp_folder+'mappings/'+strain+'/'
        for bamfile in os.listdir(mapping_folder):
            if bamfile.endswith('.dedup.bam'):
                ID = bamfile.rsplit('.',2)[0]
                if not os.path.exists(F"{mapping_folder}{ID}.cov.gz"):
                    print(F"Extracting coverage for {ID} ...", flush=True)
                    xtr = extract_coverage(F"{mapping_folder}{ID}.dedup.bam")
                    with gzip.open(F"{mapping_folder}{ID}.cov.gz", 'wt') as outfile:
                        for i in xtr:
                            outfile.write(F"{i}\n")
                cov_files[strain].append(F"{mapping_folder}{ID}.cov.gz")
    return cov_files

def extract_coverage(datei):#apply_UQ_normalization(datei):
    bamfile = pysam.AlignmentFile(datei, 'rb')
    ref_length = bamfile.lengths[0]
    dataset = [0 for i in range(ref_length)]
    for pileupcolumn in bamfile.pileup():
        dataset[pileupcolumn.reference_pos-1] = pileupcolumn.n
    bamfile.close()
    return dataset

def norm_coverage(coverage, type='max'):
    """
    Normalizes coverage
    """
    if type=='max':
        max_cov = max(coverage)
        return [c/max_cov for c in coverage]
    elif type=='median':
        median_cov = np.median(coverage)
        return [c/median_cov for c in coverage]
    elif type=='mean':
        mean_cov = np.mean(coverage)
        return [c/mean_cov for c in coverage]
    elif type=='trimUQ':
        trim = np.median(sorted(coverage)[-len(coverage)//4:])
        if trim==0: trim = 1.0
        return [c/trim if c<=trim else 1.0 for c in coverage]
    else:
        print("Unknown normalization type!")
        sys.exit(1)

def process_cov_file(cov_file, q, pindex):
    """
    Reads coverage file and returns dict: (start, end, strand): coverage
    """
    coverage = []
    with gzip.open(cov_file, 'rt') as f:
        for line in f:
            coverage.append(float(line.strip()))

    coverage = norm_coverage(coverage, type='trimUQ')

    q.put((pindex, cov_file.split('/')[-1].rsplit('_',1)[0], coverage))

def format_time(seconds):
    """converts seconds into a human readable format"""
    return time.strftime("%H:%M:%S", time.gmtime(seconds))

def read_fasta(fasta):
    ref_seq = []
    with gzip.open(fasta, 'rt') as infile:
        seq=''
        title=''
        for line in infile:
            if line[0]=='>' and title=='':
                title=line[1:].strip()
            elif line[0]=='>' and title!='':
                ref_seq.append((title,seq))
            else:
                seq+=line.strip()
        ref_seq.append((title,seq))
    return ref_seq[0][1]

def create_norm_consensus_coverage(cov_files, strain, ref_len):
    coverages = {}
    sample_count = len(cov_files)
    total_coverage = []
    q = Queue()
    thread_list = []
    pindex = 1
    pindex_count = 0
    p_read = defaults.cores
    if p_read > len(cov_files):
        p_read = len(cov_files)
    start_time = time.time()
    while len(cov_files)>0 or len(thread_list)>0:
        print_ctr = False
        while len(cov_files)>0 and len(thread_list)<p_read:
            thread_list.append((pindex, Process(target=process_cov_file, args=(cov_files.pop(), q, pindex), daemon=True)))
            thread_list[-1][1].start()
            pindex += 1
        while not q.empty():
            print_ctr = True
            pID, name, cov = q.get()
            if not name in coverages: coverages[name] = []
            coverages[name].append(cov)
            index = [i for i, t in thread_list].index(pID)
            thread_list[index][1].join(1)
            thread_list.pop(index)
            pindex_count += 1
        if print_ctr:
            print(F"Est. remaining runtime for strain {strain}: {format_time((time.time()-start_time)/pindex_count*(len(cov_files))+len(thread_list))}", flush=True, end='\r')
    print()
    
    print(F"Collapsing coverages of {strain} into single coverage...", flush=True)
    total_coverage = [sum([sum([cov[p] for cov in l]) for l in coverages.values()])/sample_count for p in range(ref_len)]
    with gzip.open(f"{defaults.temp_folder}mappings/{strain}.avg_cov.cov.gz", 'wt') as outfile:
        for c in total_coverage:
            outfile.write(F"{c}\n")
        
    return total_coverage

def ident_zero_stretches(mask):
    stretches = []
    start = None
    for i, val in enumerate(mask):
        if not val and start is None:
            start = i
        elif val and start is not None:
            if i-1 > start:
                stretches.append((start, i-1))
            start = None
    if start is not None:
        if len(mask)-1 > start:
            stretches.append((start, len(mask)-1))
    return stretches

def identify_regions(total_coverage, annotation, prophages):
    #defining masks
    anno_gene_mask_incl_promoter = np.zeros(len(total_coverage), dtype=bool)
    gene_expression_mask = np.zeros(len(total_coverage), dtype=bool)
    position_expression_mask = np.zeros(len(total_coverage), dtype=bool)
    anno_prophage_mask = np.zeros(len(total_coverage), dtype=bool)
    position_true_zero_mask = np.zeros(len(total_coverage), dtype=bool)
    
    #prophages
    for start, end in prophages:
        anno_prophage_mask[start:end+1] = True
    
    #genes annotation and expression
    for region in annotation:
        if region[2] == '+':
            start = region[0]-defaults.promoter_region
            if start < 0:
                start = 0
            end = region[1]
        else:
            start = region[0]
            end = region[1]+defaults.promoter_region
            if end > len(total_coverage):
                end = len(total_coverage)
        anno_gene_mask_incl_promoter[start:end+1] = True
        val = sum(total_coverage[start:end+1])/(end-start+1)
        if val>defaults.coverage_threshold: #mask region if expressed
            gene_expression_mask[start:end+1] = True
    
    #position expression
    for p in range(len(total_coverage)):
        if total_coverage[p] > defaults.coverage_threshold:
            position_expression_mask[p] = True
    
    #true zero positions
    for p in range(len(total_coverage)):
        if total_coverage[p] > 0: 
            position_true_zero_mask[p] = True
    
    #genes not expression and w/o prophages -> must be annotated gene, not prophage and not expressed
    not_expressed_genes = ident_zero_stretches(np.invert(anno_gene_mask_incl_promoter) + anno_prophage_mask + gene_expression_mask)
    #intergenic regions not expression and w/o prophages -> must not be annotated gene, not prophage and not position expressed
    not_expressed_intergenic = ident_zero_stretches(anno_gene_mask_incl_promoter + anno_prophage_mask + position_expression_mask)
    #intergenic regions expressed but w/o prophages -> must not be annotated gene, not prophage, but position expressed
    expressed_intergenic = ident_zero_stretches(anno_gene_mask_incl_promoter + anno_prophage_mask + np.invert(position_expression_mask))
    #true zero regions no matter if gene or not, but not in prophage
    true_zero_regions = ident_zero_stretches(position_true_zero_mask + anno_prophage_mask)
    
    with open('term_definitions.txt', 'w') as outfile:
        outfile.write('not_expressed_genes: regions annotated as genes incl. promoter region upstream, but not expressed above threshold and not within annotated prophage region\n')
        outfile.write('not_expressed_intergenic: regions not annotated as genes, but not expressed above threshold and not within annotated prophage region\n')
        outfile.write('true_zero_regions: regions with 0 coverage independent of annotation except outside of prophage region\n')
        outfile.write('expressed_intergenic: regions not annotated as genes, but expressed above threshold and not within annotated prophage region\n')
    
    return not_expressed_genes, not_expressed_intergenic, true_zero_regions, expressed_intergenic

def write_output(strain, not_expressed_genes, not_expressed_intergenic, true_zero_regions, expressed_intergenic, consensus_cov):
    print(F"Writing output for {strain}...", flush=True)
    with open(f'{defaults.output_folder}{strain}.regions.tsv', 'w') as outfile:
        outfile.write(F"start(base_1)\tend(base_1_inclusive)\tlength\tregion_type\tavg_norm_coverage\n")
        true_zero_regions = sorted([(e-s, s, e) for s, e in true_zero_regions], reverse=True)
        for l, s, e in true_zero_regions:
            if l>=defaults.min_region_size:
                outfile.write(F"{s+1}\t{e+1}\t{e-s+1}\ttrue_zero_regions\t{sum(consensus_cov[s:e+1])/l}\n")
        
        not_expressed_genes = sorted([(e-s, s, e) for s, e in not_expressed_genes], reverse=True)
        for l, s, e in not_expressed_genes:
            if l>=defaults.min_region_size:
                outfile.write(F"{s+1}\t{e+1}\t{e-s+1}\tnot_expressed_genes\t{sum(consensus_cov[s:e+1])/l}\n")
        
        not_expressed_intergenic = sorted([(e-s, s, e) for s, e in not_expressed_intergenic], reverse=True)
        for l, s, e in not_expressed_intergenic:
            if l>=defaults.min_region_size:
                outfile.write(F"{s+1}\t{e+1}\t{e-s+1}\tnot_expressed_intergenic\t{sum(consensus_cov[s:e+1])/l}\n")
        
        expressed_intergenic = sorted([(e-s, s, e) for s, e in expressed_intergenic], reverse=True)
        for l, s, e in expressed_intergenic:
            if l>=defaults.min_region_size:
                outfile.write(F"{s+1}\t{e+1}\t{e-s+1}\texpressed_intergenic\t{sum(consensus_cov[s:e+1])/l}\n")
            

#main function

def main():
    #parse config file and validate args and parse input table
    input_dataset, anno = parse_config(args.config)
    
    #copying reference files
    input_dataset = copying_ref_files(input_dataset)
    
    #preparing references for mapping
    prepare_references(input_dataset)
    
    #run mapping
    run_mapping(input_dataset)
    
    #extract coverage to file
    coverage_files = extract_coverage_to_file()
    
    #read reference fasta
    references = {strain:[read_fasta(ref), anno[ref]['anno']] for strain, ref in defaults.strains}
    
    #process each strain
    for ref in input_dataset:
        strain = ref.split('/')[-1].rsplit('.',2)[0]
        print(F"Creating normalized consensus coverage for {strain}", flush=True)
        consensus_cov = create_norm_consensus_coverage(coverage_files[strain], strain, len(references[strain][0]))
    
        #ident expressed genes and mask them incl. promoter region, mask prophage if available
        not_expressed_genes, not_expressed_intergenic, true_zero_regions, expressed_intergenic = identify_regions(consensus_cov, references[strain][1], [p for k in anno[ref] if k!='anno' for p in anno[ref][k]])
        
        #write output
        write_output(strain, not_expressed_genes, not_expressed_intergenic, true_zero_regions, expressed_intergenic, consensus_cov)
    
    if not args.keep:
        print(F"Removing temporary files ...", flush=True)
        shutil.rmtree(defaults.temp_folder)
    
    print('Done. Have a nice day!')

main()