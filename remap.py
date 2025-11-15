import pysam
import os
import subprocess
import glob
from numba import njit
import numpy as np
import vacmap_index
from multiprocessing import Process, Lock
def decode_strand_from_flag(flag: int) -> str:
    """
    Decodes the strand information from a pysam.AlignedSegment.flag value.
    Returns "+" for the forward strand, "-" for the reverse strand.
    
    Args:
        flag (int): The flag value from a pysam.AlignedSegment object.

    Returns:
        str: "+" for forward strand, "-" for reverse strand.
    """
    # In SAM format, bit 0x10 (16 decimal) means the read is mapped to the reverse strand
    if flag & 16:
        return "-"
    else:
        return "+"
@njit
def decode_cigar(contig, strand, read_length, cigarstring, refloc, truereadpos = False):
    useref = set(('M', 'D', 'N', '=', 'X'))
    if(truereadpos):
        useread = set(('M', 'I', 'S', '=', 'X', 'H'))
    else:
        useread = set(('M', 'I', 'S', '=', 'X'))
    match = set(('M', '='))
    numset = set(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'))
    op2num = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'X': 9}
    num = 0
    meet_match = False
    readloc = 0
    cigartuple = []
    
    ref_st = -1
    read_st = -1
    
    for c in cigarstring:
        if(c in numset):
            num = num * 10 + ord(c) - ord('0')
        else:
            
            if(ref_st == -1):
                if(c in useref):
                    ref_st = refloc
            if(read_st == -1):
                if(c in ('M', 'I', '=', 'X')):
                    read_st = readloc
                    
                    
            if(c in useref):
                refloc += num

                ref_en = refloc
                

            if(c in useread):
                readloc += num
                if(c in ('M', 'I', '=', 'X')):
                    read_en = readloc
            cigartuple.append((op2num[c], num))
            num = 0
    if(strand == '+'):
        return contig, strand, read_st, read_en, ref_st, ref_en
    else:
        read_length = readloc
        return contig, strand, read_length - read_en, read_length - read_st, ref_st, ref_en
def get_third(x):
    return x[2]
def merge_fastq(files, output, checkunique = False):
    seen = set()
    with open(output, 'w') as out:
        for fname in files:
            with open(fname) as f:
                while True:
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline()
                    plus = f.readline()
                    qual = f.readline()
                    name = header.split()[0]  # FASTQ header, e.g. '@SEQ_ID'
                    if name not in seen:
                        out.write(header)
                        out.write(seq)
                        out.write(plus)
                        out.write(qual)
                        seen.add(name)
                    elif(checkunique == True):
                        print('Duplication detected', name)
def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of a DNA sequence.
    Non-ATCG (case-insensitive) bases are converted to 'N'.

    Args:
        seq (str): DNA sequence.

    Returns:
        str: Reverse complement sequence with non-ATCG bases as 'N'.
    """
    base_map = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c'
    }

    # Build the reverse complement
    rc = []
    for base in reversed(seq):
        rc.append(base_map.get(base, 'N'))
    return ''.join(rc)
def reverse_complement(seq: str) -> str:
    return Seq(seq).reverse_complement()
def removefile(pattern, exclude = ''):
    # Find files matching the pattern
    if(isinstance(pattern, str)):
        files_to_delete = glob.glob(pattern)
    else:
        files_to_delete = pattern
    # Remove each found file
    for file_path in files_to_delete:
        try:
            if(exclude != '' and exclude in file_path.split('/')[-1]):
                continue
            os.remove(file_path)
            #print(f"Removed: {file_path}")
        except OSError as e:
            print(f"Error removing {file_path}: {e}")
import glob
import os
import shutil

def removefile(pattern, exclude=''):
    """
    Remove files and folders matching the pattern.
    If exclude is set, skip any file/folder whose name contains exclude.
    """
    if isinstance(pattern, str):
        files_to_delete = glob.glob(pattern)
    else:
        files_to_delete = pattern

    for file_path in files_to_delete:
        name = os.path.basename(file_path)
        if exclude and exclude in name:
            continue
        try:
            if os.path.isdir(file_path):
                shutil.rmtree(file_path)
                # print(f"Removed folder: {file_path}")
            else:
                os.remove(file_path)
                # print(f"Removed file: {file_path}")
        except OSError as e:
            print(f"Error removing {file_path}: {e}")
def bam2fastq(bamfilepath, savepath, contig, st, en):
    bamfile = pysam.AlignmentFile(bamfilepath, "rb")
    with open(savepath, "w") as ofile:
        usedreadid = set()

        for read in bamfile.fetch(contig, st, en):
            read_name = read.query_name
            if(read_name not in usedreadid):
                seq = read.query_sequence
                qual = read.query_qualities

                # Convert quality scores to ASCII
                if qual is not None:
                    qual = ''.join(chr(q + 33) for q in qual)
                else:
                    qual = 'I' * len(seq)  # Default quality if none available

                # Format FASTQ record
                fastq_record = f"@{read_name}\n{seq}\n+\n{qual}\n"
                ofile.write(fastq_record)
                usedreadid.add(read_name)
    return usedreadid
def gfa_to_fasta(gfa_path_list, fasta_path, id_prefix="noname"):
    """
    Convert a GFA file to a FASTA file.
    Each segment line ('S') in GFA becomes a FASTA entry.

    Args:
        gfa_path (str): Path to the input GFA file.
        fasta_path (str): Path to the output FASTA file.
        id_prefix (str): Optional prefix for segment IDs if GFA 'S' lines lack IDs.
    """
    if(isinstance(gfa_path_list, str)):
        gfa_path_list = [gfa_path_list]
    nonameindex = 1
    count = 0
    with open(fasta_path, 'w') as fasta:
        for gfa_path in gfa_path_list:
            with open(gfa_path, 'r') as gfa:
                for line in gfa:
                    if line.startswith('S'):
                        parts = line.strip().split('\t')
                        # GFA: S\t<name>\t<sequence>\t...
                        if len(parts) >= 3:
                            seg_id = parts[1]
                            seq = parts[2]
                            fasta.write(f">{seg_id}\n{seq}\n")
                            count += 1
    return count
                    
def load_gfa(gfa_path, isfasta = False):
    read2contigpos = dict()
    id2contig = dict()
    asmid2queryname = dict()
    if(isfasta == False):
        with open(gfa_path, 'r') as gfa:
            for line in gfa:
                if line.startswith('A'):
                    parts = line.strip().split('\t')
                    readname = parts[4]
                    if(readname not in read2contigpos):
                        read2contigpos[readname] = parts[1: -2]
                        if(parts[1] not in asmid2queryname):
                            asmid2queryname[parts[1]] = [readname]
                        else:
                            asmid2queryname[parts[1]].append(readname)
                if line.startswith('S'):
                    parts = line.strip().split('\t')
                    # GFA: S\t<name>\t<sequence>\t...
                    if len(parts) >= 3:
                        seg_id = parts[1]
                        seq = parts[2]
                        id2contig[seg_id] = seq
    else:
        for rec in vacmap_index.fastx_read(gfa_path):
            id2contig[rec[0]] = rec[1]
            
    
    if(len(asmid2queryname) == 0):
        for asmid in id2contig:
            asmid2queryname[asmid] = []
    return read2contigpos, id2contig, asmid2queryname

@njit
def cigar_eqx(target, query, cigar):
    ref_pos = 0
    query_pos = 0
    new_cigar = []
    i = 0
    zero = ord('0')
    nine = ord('9')
    newcigar = ''
    while i < len(cigar):
        num = 0
        while i < len(cigar) and zero <= cigar[i] <= nine:
            num = num * 10 + (cigar[i] - zero)
            i += 1
        if i >= len(cigar):
            break
        op = chr(cigar[i])
        i += 1
        if op == 'M':
            if num == 0:
                continue
            current_run = 0
            current_type = None
            for _ in range(num):
                base_t = target[ref_pos]
                base_q = query[query_pos]
                # Uppercase comparison: A=65, a=97; subtract 32 from lowercase
                base_t_upper = base_t - 32 if 97 <= base_t <= 122 else base_t
                base_q_upper = base_q - 32 if 97 <= base_q <= 122 else base_q
                this_type = '=' if base_t_upper == base_q_upper else 'X'
                if this_type == current_type:
                    current_run += 1
                else:
                    if current_type is not None:
                        new_cigar.append(str(current_run))
                        new_cigar.append(current_type)
                    current_type = this_type
                    current_run = 1
                ref_pos += 1
                query_pos += 1
            if current_type is not None:
                new_cigar.append(str(current_run))
                new_cigar.append(current_type)
        else:
            new_cigar.append(str(num))
            new_cigar.append(op)
            if op in 'MDN=X':
                ref_pos += num
            if op in 'MIS=X':
                query_pos += num
    for c in new_cigar:
        newcigar += c
    return newcigar

def bam2fastq_regions(bamfilepath, savepath, regions, keeprawread, verbose):
    
    regionfastq_queryid_list = []
    for region in regions:
        contig, st, en = region[0], int(region[1]), int(region[2])
        if(verbose):
            print(contig+'_'+str(st)+'_'+str(en))
        
        reads_name = savepath+'region_'+contig+'_'+str(st)+'_'+str(en)
        if(keeprawread == False):
            queryset = bam2fastq(bamfilepath, reads_name+'.tmp.fastq', contig, st, en)
            regionfastq_queryid_list.append([reads_name+'.tmp.fastq', queryset])
        else:
            queryset = bam2fastq(bamfilepath, reads_name+'.raw.fastq', contig, st, en)
            regionfastq_queryid_list.append([reads_name+'.raw.fastq', queryset])
    return regionfastq_queryid_list

def local_asm(bamfilepath, savepath, contig, st, en, keeprawread = False):
    
    a = subprocess.run(['hifiasm', '-o', reads_name+'_tmp', '-t1', '-f0', reads_name+'.raw.fastq'], capture_output=True)
    if(a.returncode == 0):
        gfa_path = reads_name+'_tmp.bp.p_utg.gfa'
        fasta_path = reads_name+'.asm.fasta'
        gfa_to_fasta(gfa_path, fasta_path)
        if(keeprawread == False):
            fastq_path = reads_name+'.tmp.fastq'
        else:
            fastq_path = reads_name+'.raw.fastq'
        read2seq = dict()
        for rec in vacmap_index.fastx_read(fastq_path):
            read2seq[rec[0]] = rec[1]
        removefile(reads_name+'*tmp*')
        
        read2contigpos, asmid2contig, asmid2queryname = load_gfa(gfa_path)
        
        for asmid in asmid2queryname:
            asmseq = asmid2contig[asmid]
            asmseq_index_object = mappy.Aligner(seq = asmseq, w = 10, k = 15)  # load index
            if not asmseq_index_object: raise Exception("ERROR: failed to load index")
            for query_name in asmid2queryname[asmid]:
                if(read2contigpos[query_name][2] == '+'):
                    query = read2seq[query_name]
                else:
                    query = reverse_complement(read2seq[query_name])
                check_asmread(asmseq_index_object, asmseq, read2contigpos, query_name, query)
        
    else:
        print('Failed to asm '+ 'region_'+contig+'_'+str(st)+'_'+str(en))
        removefile(reads_name+'*tmp*')
def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of a DNA sequence.
    Non-ATCG (case-insensitive) bases are converted to 'N'.

    Args:
        seq (str): DNA sequence.

    Returns:
        str: Reverse complement sequence with non-ATCG bases as 'N'.
    """
    base_map = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c'
    }

    # Build the reverse complement
    rc = []
    for base in reversed(seq):
        rc.append(base_map.get(base, 'N'))
    return ''.join(rc)
from Bio.Seq import Seq


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())
from collections import defaultdict
from typing import Any, List, Set

class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x: int, y: int) -> None:
        px, py = self.find(x), self.find(y)
        if px != py:
            if self.rank[px] < self.rank[py]:
                self.parent[px] = py
            elif self.rank[px] > self.rank[py]:
                self.parent[py] = px
            else:
                self.parent[py] = px
                self.rank[px] += 1

def group_named_sets(setpack: List[List[Any]]) -> List[List[List[Any]]]:
    """
    Groups a list of [name, set] pairs into groups where sets within a group overlap 
    (directly or indirectly), and groups are disjoint from each other.
    
    Args:
    setpack: A list of [setname, set] where set is a Python set.
    
    Returns:
    A list of lists, each containing the [setname, set] pairs in that group.
    """
    n = len(setpack)
    if n == 0:
        return [], []
    
    uf = UnionFind(n)
    elem_to_indices = defaultdict(list)
    
    for i in range(n):
        name, s = setpack[i]
        if not isinstance(s, set):
            raise ValueError(f"Item at index {i} is not a set.")
        for elem in s:
            elem_to_indices[elem].append(i)
    
    for indices in elem_to_indices.values():
        for j in range(1, len(indices)):
            uf.union(indices[0], indices[j])
    
    root_to_group = defaultdict(list)
    for i in range(n):
        root = uf.find(i)
        root_to_group[root].append(setpack[i])
    
    # Optional: sort groups by the order of their first appearance
    groups = list(root_to_group.values())
    groups.sort(key=lambda g: min(setpack.index(item) for item in g))
    
    # Optional: sort within groups by original order
    for group in groups:
        group.sort(key=lambda item: setpack.index(item))

    groups_size = []
    for group in groups:
        list_of_sets = []
        for pack in group:
            list_of_sets.append(pack[1])
        groups_size.append(len(set().union(*list_of_sets)))
    return groups, groups_size
def check_overlap_and_merge(regionfastq_queryid_list):
    merge_index = 1
    newpaths = []
    merged_regionfastq_queryid_list, groups_size = group_named_sets(regionfastq_queryid_list)
    for item in merged_regionfastq_queryid_list:
        overlaps = [line[0] for line in item]
        if(len(overlaps) > 1):
            output = savepath + 'merged_'+ str(merge_index) + '.fastq'
            merge_index += 1
            merge_fastq(overlaps, output)
            removefile(overlaps)
            newpaths.append(output)
        else:
            newpaths.append(overlaps[0])
    return newpaths, groups_size  


def indexasmid(iloc, tmp_read2contigpos, tmp_asmid2contig, tmp_asmid2queryname, tmp_missed_queryname2alns, read2contigpos, asmid2contig, asmid2queryname, missed_queryname2alns):
    asmid2index = dict()
    for asmid in tmp_asmid2contig:
        asmid2index[asmid] = iloc
        iloc += 1
    for query_name in tmp_read2contigpos:
        assert query_name not in read2contigpos
        tmp_read2contigpos[query_name][0] = asmid2index[tmp_read2contigpos[query_name][0]]
        read2contigpos[query_name] = tmp_read2contigpos[query_name]
    for asmid in tmp_asmid2contig:
        asmid2contig[asmid2index[asmid]] = tmp_asmid2contig[asmid]
        asmid2queryname[asmid2index[asmid]] = tmp_asmid2queryname[asmid]
        
    for query_name in tmp_missed_queryname2alns:
        missed_queryname2alns[query_name] = []
        for aln in tmp_missed_queryname2alns[query_name]:
            aln[0] = asmid2index[aln[0]]
            missed_queryname2alns[query_name].append(aln)
        
    return iloc
def write_fasta_from_dict(seq_dict, fasta_path):
    """
    Writes a FASTA file from a dictionary of {readname: readseq} pairs.

    Parameters:
        seq_dict (dict): Dictionary with read names as keys and sequences as values.
        fasta_path (str): Path to output fasta file.
    """
    with open(fasta_path, 'w') as f:
        for name, seq in seq_dict.items():
            f.write(f">{name}\n{seq}\n")
def align_seqs(alignername, refpath, asmseqs_path, out_alns_path, out_sorted_alns_path, threads, datatype):
    if(alignername == 'vacmap'):
        
        a = subprocess.run(['vacmap', '-mode', 'L', '-ref', refpath, '-read', asmseqs_path, '-t', str(threads), '-o', out_alns_path, '--force', '--eqx'], capture_output=True)
    else:
        if(datatype == 'asm'):
            preset = 'asm5'
        elif(datatype == 'hifi'):
            preset = 'map-hifi'
        elif(datatype == 'ont'):
            preset = 'map-ont'
        a = subprocess.run(['minimap2', '-aYx', preset, '-I', '1000G', refpath, asmseqs_path, '-t', str(threads), '-o', out_alns_path, '--eqx'], capture_output=True)
    if(a.returncode == 0):
        a = subprocess.run(['samtools', 'sort', '-@8', out_alns_path, '-o', out_sorted_alns_path], capture_output=True)
        if(a.returncode == 0):
            a = subprocess.run(['samtools', 'index', '-@8', out_sorted_alns_path], capture_output=True)
            if(a.returncode == 0):
                return True

    print(a.key)

    return False
def load_asm_alns(bamfilepath):

    bamfile = pysam.AlignmentFile(bamfilepath, "rb")
    index2contig = bamfile.header.references
    
    contig2index = dict()
    for iloc in range(len(index2contig)):
        contig2index[index2contig[iloc]] = iloc
    asmid2alns = dict()
    for AlignedSegment in bamfile:
        if(AlignedSegment.is_secondary == True):
            continue
        if(AlignedSegment.is_mapped == False):
            continue
        asmid = AlignedSegment.query_name
        contig = AlignedSegment.reference_name
        strand = decode_strand_from_flag(AlignedSegment.flag)
        read_length = AlignedSegment.query_length
        aln = list(decode_cigar(contig, strand, read_length, AlignedSegment.cigarstring, AlignedSegment.reference_start))
        aln.append(AlignedSegment.cigarstring)
        if(asmid in asmid2alns):
            asmid2alns[asmid].append(aln)
        else:
            asmid2alns[asmid] = [aln]
    for asmid in asmid2alns:
        asmid2alns[asmid].sort(key = get_third)
    bamheader = bamfile.header.to_dict()
    bamheader['PG'].append({'ID': 'Local_asm_guided_realignment', 'PN': 'Local_asm_guided_realignment', 'PP': bamheader['PG'][-1]['PN'],
                           'VN': '1.0', 'CL': 'No'})
    #bamheader = pysam.AlignmentHeader.from_dict(bamheader)
    return asmid2alns, index2contig, contig2index, bamheader
def load_sam(bamfilepath):

    bamfile = pysam.AlignmentFile(bamfilepath, "r")
    asmid2alns = dict()
    for AlignedSegment in bamfile:
        if((AlignedSegment.is_secondary == True) or (AlignedSegment.is_mapped == False)):
            continue
        asmid = AlignedSegment.query_name
        contig = AlignedSegment.reference_name
        strand = decode_strand_from_flag(AlignedSegment.flag)
        read_length = AlignedSegment.query_length
        aln = list(decode_cigar(contig, strand, read_length, AlignedSegment.cigarstring, AlignedSegment.reference_start))
        aln[2] = 0
        aln[3] = read_length
        aln.append(AlignedSegment.cigarstring)
        if(asmid in asmid2alns):
            asmid2alns[asmid].append(aln)
        else:
            asmid2alns[asmid] = [aln]
    for asmid in asmid2alns:
        asmid2alns[asmid].sort(key = get_third)
    return asmid2alns
from cigar import Cigar
def change_alns_strand(alns, query_len):
    newalns = []
    for aln in alns:
        contig, strand, q_st, q_en, r_st, r_en, cigar, mapq = aln
        q_en, q_st = query_len - q_st, query_len - q_en
        if(strand == '+'):
            strand = '-'
        else:
            strand = '+'
        newalns.append((contig, strand, q_st, q_en, r_st, r_en, cigar, mapq))
    return newalns
def change_alns_coo_baseon_strand(alns, query_len):
    newalns = []
    for aln in alns:
        contig, strand, q_st, q_en, r_st, r_en, cigar, mapq = aln

        if(strand == '-'):
            q_en, q_st = query_len - q_st, query_len - q_en

        newalns.append((contig, strand, q_st, q_en, r_st, r_en, cigar, mapq))
    return newalns

@njit
def fill_refpos_forward(onepack, project_arr, project_arr_contig, contig_index, cigar):
    
    refpos = onepack[4]
    readpos = 0
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):

                project_arr[readpos: readpos + num] = np.arange(refpos, refpos + num)
                project_arr_contig[readpos: readpos + num] = contig_index

            if(c in use_query):
                readpos += num
            if(c in use_ref):
                refpos += num
            num = 0

@njit
def fill_refpos_revse(onepack, project_arr, project_arr_contig, contig_index, cigar):
    
    refpos = onepack[4]
    readpos = len(project_arr)
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):
                project_arr[readpos - num: readpos] = np.arange(-refpos - num + 1, -refpos + 1)
                project_arr_contig[readpos - num: readpos] = contig_index
            if(c in use_query):
                readpos -= num
            if(c in use_ref):
                refpos += num
            num = 0      

@njit
def fill_readpos_forward(onepack, project_arr, cigar):
    
    read_project_arr = np.empty(len(project_arr))
    read_project_arr.fill(np.nan)
    

    refpos = onepack[4]
    readpos = onepack[2]
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):
                #print(refpos, refpos + num, readpos, readpos + num)
                read_project_arr[refpos: refpos + num] = np.arange(readpos, readpos + num)
            if(c in use_query):
                readpos += num
            if(c in use_ref):
                refpos += num
            num = 0
    return read_project_arr
@njit
def get_new_alignment_position(read_project_arr, project_arr, project_arr_contig, index2contig, maxreadgap = 50, maxrefgap = 50, maxgap = 50, minaln_size = 40):
    alns = []
    pre_strand = ''
    len_project_arr = len(project_arr)
    for iloc in range(len(project_arr)):
        if(~np.isnan(read_project_arr[iloc]) and ~np.isnan(project_arr[iloc])):
            q_st = int(read_project_arr[iloc])
            pre_contig = project_arr_contig[iloc]
            pre_readpos = q_st
            if(project_arr[iloc] > 0):
                pre_strand = '+'
                r_st = int(project_arr[iloc])
                pre_refpos = r_st
            else:
                pre_strand = '-'
                r_en = 1 - int(project_arr[iloc])
                pre_refpos = r_en - 1
            break
    
    if(pre_strand == ''):
        return alns
    startiloc = iloc + 1
    for iloc in range(startiloc, len_project_arr):
        if(~np.isnan(read_project_arr[iloc]) and ~np.isnan(project_arr[iloc])):
            contig = project_arr_contig[iloc]
            readpos = int(read_project_arr[iloc])
            
            refpos = int(project_arr[iloc])
            if(refpos > 0):
                strand = '+'
            else:
                strand = '-'
                refpos = -refpos
            if(pre_strand == strand):
                if(contig == pre_contig):
                    readgap = readpos - pre_readpos
                    if(strand == '+'):      
                        refgap = refpos - pre_refpos
                    else:
                        refgap = pre_refpos - refpos
                    if(readgap <= maxreadgap):
                        
                        if(refgap <= maxrefgap):
                            gap = readgap - refgap
                            if(gap < maxgap):
                                pre_readpos = readpos
                                pre_refpos = refpos
                                continue
                    elif((readgap <= 100) and (0 <= refgap < 10000)):
                        pre_readpos = readpos
                        pre_refpos = refpos
                        continue
                    elif((0 <= refgap <= 100) and (readgap < 10000)):
                        pre_readpos = readpos
                        pre_refpos = refpos
                        continue
                            
                        
                    
            if(pre_strand == '+'):
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
            else:
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
            q_st = readpos
            if(strand == '+'):
                r_st = refpos
            else:
                r_en = refpos + 1
            pre_contig, pre_strand, pre_readpos, pre_refpos = contig, strand, readpos, refpos
            if((onepack[3] - onepack[2]) > minaln_size):
                alns.append(onepack)
    if(pre_strand == '+'):
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
    else:
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
    if((onepack[3] - onepack[2]) > minaln_size):
        alns.append(onepack)
    return alns
def extend_top(onealn, extend_st, index_object, query):
    #print('top')
    #print(onealn)
    if(onealn[1] == '+'):
        q_st = onealn[2]
        r_st = onealn[4]
        contig = onealn[0]
        if(extend_st == q_st):
            return onealn
        extend_seq_query = query[extend_st: q_st][::-1]
        extend_seq_ref = index_object.seq(contig, max(r_st - q_st + extend_st, 0), r_st)[::-1]
        extend_size = 0
        for iloc in range(q_st - extend_st):
            if(extend_seq_query[iloc] == extend_seq_ref[iloc]):
                extend_size += 1
            else:
                break

        onealn = (onealn[0], onealn[1], onealn[2] - extend_size, onealn[3], onealn[4] - extend_size, onealn[5])
    else:
        q_st = onealn[2]
        r_en = onealn[5]
        contig = onealn[0]
        if(extend_st == q_st):
            return onealn
        extend_seq_query = reverse_complement(query[extend_st: q_st])
        extend_seq_ref = index_object.seq(contig, r_en, r_en + q_st - extend_st)
        extend_size = 0
        for iloc in range(q_st - extend_st):
            if(extend_seq_query[iloc] == extend_seq_ref[iloc]):
                extend_size += 1
            else:
                break

        onealn = (onealn[0], onealn[1], onealn[2] - extend_size, onealn[3], onealn[4], onealn[5] + extend_size)
    #print(extend_seq_query)
    #print(extend_seq_ref)
    #print(onealn)
    #print()
    return onealn
def extend_tail(onealn, extend_en, index_object, query):
    #print('tail')
    #print(onealn)
    if(onealn[1] == '+'):
        q_en = onealn[3]
        r_en = onealn[5]
        contig = onealn[0]
        if(q_en == extend_en):
            return onealn
        extend_seq_query = query[q_en: extend_en]
        extend_seq_ref = index_object.seq(contig, r_en, r_en + extend_en - q_en)
        extend_size = 0
        for iloc in range(extend_en - q_en):
            if(extend_seq_query[iloc] == extend_seq_ref[iloc]):
                extend_size += 1
            else:
                break

        onealn = (onealn[0], onealn[1], onealn[2], onealn[3] + extend_size, onealn[4], onealn[5] + extend_size)
    else:
        q_en = onealn[3]
        r_st = onealn[4]
        contig = onealn[0]
        if(q_en == extend_en):
            return onealn
        extend_seq_query = query[q_en: extend_en]
        extend_seq_ref = reverse_complement(index_object.seq(contig, max(r_st - extend_en + q_en, 0), r_st))
        extend_size = 0
        for iloc in range(extend_en - q_en):
            if(extend_seq_query[iloc] == extend_seq_ref[iloc]):
                extend_size += 1
            else:
                break

        onealn = (onealn[0], onealn[1], onealn[2], onealn[3] + extend_size, onealn[4] - extend_size, onealn[5])
    #print(extend_seq_query)
    #print(extend_seq_ref)
    #print(onealn)
    #print()
    return onealn
def extend_edge(alns, index_object, query):
    iloc = 0
    onealn = alns[iloc]
    extend_st = 0
    onealn = extend_top(onealn, extend_st, index_object, query)
    if(len(alns) > 1):

        for iloc in range(1, len(alns)):
            extend_en = alns[iloc][2]
            if(onealn[3] < extend_en):
                onealn = extend_tail(onealn, extend_en, index_object, query)
            alns[iloc-1] = onealn
            
            onealn = alns[iloc]
            extend_st = alns[iloc-1][3]
            if(onealn[2] > extend_st):
                onealn = extend_top(onealn, extend_st, index_object, query)
            
    extend_en = len(query)  
    if(onealn[3] < extend_en):
        onealn = extend_tail(onealn, extend_en, index_object, query)
    alns[-1] = onealn
    for iloc in range(len(alns)):
        aln = alns[iloc]
        if(aln[1] == '-'):
            tmpquery = reverse_complement(query[aln[2]: aln[3]])
        else:
            tmpquery = query[aln[2]: aln[3]]
        tmpref = index_object.seq(aln[0], aln[4], aln[5])
        try:
            info = vacmap_index.k_cigar(tmpref, tmpquery, bw=50, zdropvalue=-1, eqx = True)
        except:
            print('len(tmpref), len(tmpquery)', len(tmpref), len(tmpquery))
            info = vacmap_index.k_cigar(tmpref, tmpquery, bw=50, zdropvalue=-1, eqx = True)
        if(isinstance(info, type(None))):
            print('len(tmpref), len(tmpquery)', len(tmpref), len(tmpquery))
        if(len(Cigar(info[0])) != len(tmpquery)):
            info = vacmap_index.k_cigar(tmpref, tmpquery, bw=-1, zdropvalue=-1, eqx = True)
        if(len(Cigar(info[0])) != len(tmpquery)):
            print('Failed to compute cigar in extend_edge')
 
        alns[iloc] = (aln[0], aln[1], aln[2], aln[3], aln[4], aln[5], info[0])

def get_asmidsize_list(asmid2contig):
    def get_second(x):
        return x[1]
    asmidsize_list = []
    for asmid in asmid2contig:
        if(asmid in asmid2queryname):
            size = len(asmid2queryname[asmid])
        else:
            size = 0
        if(asmid in missed_asmid2queryname):
            size += len(asmid2queryname[asmid])
        asmidsize_list.append([asmid, size])
    asmidsize_list.sort(key = get_second)
    return asmidsize_list[::-1]


def get_readnames(fastq_path):
    return set([rec[0] for rec in vacmap_index.fastx_read(fastq_path)])
def local_asm_onefile(fastq_path_queue, lock, threads, large_size = 10000):
    while(True):
        fastq_path, fastq_size = fastq_path_queue.get()
        if(isinstance(fastq_path, int)):
            break
        if(fastq_size > large_size):
            lock.acquire()
            a = subprocess.run(['hifiasm', '-o', fastq_path[:-5]+'tmp', '-t'+str(threads), '-f0', fastq_path], capture_output=True)
            lock.release()
        else:

            a = subprocess.run(['hifiasm', '-o', fastq_path[:-5]+'tmp', '-t'+str(1), '-f0', fastq_path], capture_output=True)

        if(a.returncode == 0):

            removefile(fastq_path[:-5]+'tmp*', exclude = 'bp.p_utg.gfa')
            readset = get_readnames(fastq_path)
            gfa_path = fastq_path[:-5]+'tmp.bp.p_utg.gfa'
            tmp_read2contigpos, tmp_asmid2contig, tmp_asmid2queryname = load_gfa(gfa_path)
            drop_readset = set(tmp_read2contigpos.keys())
            readset = list(readset - drop_readset)
            selected_fastq_path = fastq_path[:-5]+'selected.fastq'
            write_selected_fastq(fastq_path, selected_fastq_path, readset)

            asm_fasta_path = fastq_path[:-5]+'.asm.fasta'
            write_fasta_from_dict(tmp_asmid2contig, asm_fasta_path)
            selected_sam_path = selected_fastq_path[:-5] + 'sam'
            if(fastq_size > large_size):
                lock.acquire()
                a = subprocess.run(['minimap2', '-aYx', 'map-ont', '-t', str(threads), '--eqx', asm_fasta_path, selected_fastq_path, '-o', selected_sam_path], capture_output=True)
                lock.release()
            else:
                a = subprocess.run(['minimap2', '-aYx', 'map-ont', '-t', str(1), '--eqx', asm_fasta_path, selected_fastq_path, '-o', selected_sam_path], capture_output=True)

            if(a.returncode == 0):
                continue

        removefile(fastq_path[:-5]+'tmp*')
        print(a)
        continue
def get_alns_for_asmid(asmid, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_object, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig, queryname2alns):
    n_aln = 0
    asmseq = asmid2contig[asmid]
    pack = asmid2alns[str(asmid)]
    project_arr = np.empty(len(asmseq))
    project_arr_contig = np.empty(len(asmseq), dtype = np.int32)
    project_arr.fill(np.nan)

    for onepack in pack:
        if(onepack[1] == '-'):
            fill_refpos_revse(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        else:
            fill_refpos_forward(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
    for query_name in asmid2queryname[int(asmid)]:

        query = all_rawseqs_object[query_name].seq
        if(read2contigpos[query_name][2] == '-'):
            query = reverse_complement(query)
            need_reverse = True
        else:
            need_reverse = False
        aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
        if(isinstance(aln_asm, tuple)):
            read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
            alns = get_new_alignment_position(read_project_arr, project_arr, project_arr_contig, index2contig)  
            extend_edge(alns, index_object, query)
            if(need_reverse == True):
                queryname2alns[query_name] = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
            else:
                queryname2alns[query_name] = change_alns_coo_baseon_strand(alns, len(query))
            n_aln += len(alns)
    
    for query_name in missed_asmid2queryname[int(asmid)]:
        aln_asm = missed_queryname2alns[query_name][0]
        query = all_rawseqs_object[query_name].seq
        if(aln_asm[1] == '-'):
            query = reverse_complement(query)
            need_reverse = True
        else:
            need_reverse = False
        read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
        alns = get_new_alignment_position(read_project_arr, project_arr, project_arr_contig, index2contig)  
        extend_edge(alns, index_object, query)
        if(need_reverse == True):
            queryname2alns[query_name] = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
        else:
            queryname2alns[query_name] = change_alns_coo_baseon_strand(alns, len(query))
        n_aln += len(alns)
    print('asmid =',asmid, 'n_aln =', n_aln)
def queryname2alns_to_bam(out_path, queryname2alns, all_rawseqs_object, header):
    if isinstance(header, dict):
        header = pysam.AlignmentHeader.from_dict(header)
    with pysam.AlignmentFile(out_path, "wb", header=header, threads=16) as outf:
        for queryname in queryname2alns:
            #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
            #      0            1    2   3    4      5         6      7
            #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
            mapinfo = []
            query, qual = all_rawseqs_object[queryname].seq, all_rawseqs_object[queryname].qual
            for line in queryname2alns[queryname]:
                contig, strand, q_st, q_en, r_st, r_en, cigar = line
                mapq = 60
                if(q_st > 0):
                    cigar = str(q_st) + 'S' + cigar
                if((len(query) - q_en) > 0):
                    cigar += str(len(query) - q_en) + 'S'
                mapinfo.append([queryname, contig, strand, q_st, q_en, r_st, r_en, mapq, cigar])

            for onealn in get_bam_dict_str(mapinfo, query, qual):
                aligned = parse_sam_line(onealn, header)
                outf.write(aligned)
def yield_alnstr(onemapinfo, all_rawseqs_object, queryname):
    mapinfo = []
    query, qual = all_rawseqs_object[queryname].seq, all_rawseqs_object[queryname].qual
    for line in onemapinfo:
        contig, strand, q_st, q_en, r_st, r_en, cigar = line
        mapq = 60
        if(q_st > 0):
            cigar = str(q_st) + 'S' + cigar
        if((len(query) - q_en) > 0):
            cigar += str(len(query) - q_en) + 'S'
        mapinfo.append([queryname, contig, strand, q_st, q_en, r_st, r_en, mapq, cigar])
        if(len(Cigar(cigar)) != len(query)):
            print(len(query), len(Cigar(cigar)))
            print(queryname)
            print(contig, strand, q_st, q_en, r_st, r_en)
            print(cigar)
            print()

    for onealn in get_bam_dict_str(mapinfo, query, qual):
        yield onealn

def yield_aln_for_asmid(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig, debug = False):  
    def get_third(x):
        return x[2]
    st = time.time()
    n_aln = 0
    asmseq = asmid2contig[asmid]
    pack = asmid2alns[str(asmid)]
    project_arr = np.empty(len(asmseq))
    project_arr_contig = np.empty(len(asmseq), dtype = np.int32)
    project_arr.fill(np.nan)
    
    pack.sort(key = get_third)
    overlapping = False
    read_end = 0
    for onepack in pack:
        if(onepack[2] < read_end):
            overlapping = True
            break
        read_end = onepack[3]
        if(onepack[1] == '-'):
            fill_refpos_revse(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        else:
            fill_refpos_forward(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
    if(debug == True):
        print('project_arr', time.time() - st)
    if(overlapping == True):
        for alnstr in yield_aln_for_asmid_overlap(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):
            yield alnstr
        return
    
    if(int(asmid) in asmid2queryname):
        for query_name in asmid2queryname[int(asmid)]:

            query = all_rawseqs_object[query_name].seq
            if(read2contigpos[query_name][2] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False

            aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
            if(isinstance(aln_asm, tuple)):
                read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                alns = get_new_alignment_position(read_project_arr, project_arr, project_arr_contig, index2contig) 
                if(len(alns) > 0):
                    extend_edge(alns, index_object, query)
                    if(need_reverse == True):
                        onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                    else:
                        onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                    for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                        yield alnstr
                    n_aln += len(alns)

    if(int(asmid) in missed_asmid2queryname):
        for query_name in missed_asmid2queryname[int(asmid)]:
            aln_asm = missed_queryname2alns[query_name][0]
            query = all_rawseqs_object[query_name].seq
            if(aln_asm[1] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False
            read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
            alns = get_new_alignment_position(read_project_arr, project_arr, project_arr_contig, index2contig) 
            if(len(alns) > 0):
                extend_edge(alns, index_object, query)
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                    yield alnstr
                n_aln += len(alns)
    print('asmid =',asmid, 'n_aln =', n_aln)
def yield_aln_for_asmid_overlap(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  
    #print('overlap')
    n_aln = 0
    asmseq = asmid2contig[asmid]
    pack = asmid2alns[str(asmid)]
    project_arr_list = []
    project_arr_contig_list = []

    for onepack in pack:
        project_arr = np.empty(len(asmseq))
        project_arr_contig = np.empty(len(asmseq), dtype = np.int32)
        project_arr.fill(np.nan)
        if(onepack[1] == '-'):
            fill_refpos_revse(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        else:
            fill_refpos_forward(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        project_arr_list.append(project_arr)
        project_arr_contig_list.append(project_arr_contig)

    if(int(asmid) in asmid2queryname):
        for query_name in asmid2queryname[int(asmid)]:

            query = all_rawseqs_object[query_name].seq
            if(read2contigpos[query_name][2] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False

            aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
            if(isinstance(aln_asm, tuple)):
                alns = []
                for project_arr, project_arr_contig in zip(project_arr_list, project_arr_contig_list):
                    read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                    for aln in get_new_alignment_position(read_project_arr, project_arr, project_arr_contig, index2contig):
                        alns.append(aln)
                if(len(alns) > 0):
                    extend_edge(alns, index_object, query)
                    if(need_reverse == True):
                        onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                    else:
                        onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                    for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                        yield alnstr
                    n_aln += len(alns)

    if(int(asmid) in missed_asmid2queryname):
        for query_name in missed_asmid2queryname[int(asmid)]:
            aln_asm = missed_queryname2alns[query_name][0]
            query = all_rawseqs_object[query_name].seq
            if(aln_asm[1] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False
            alns = []
            for project_arr, project_arr_contig in zip(project_arr_list, project_arr_contig_list):
                read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                for aln in get_new_alignment_position(read_project_arr, project_arr, project_arr_contig, index2contig):
                    alns.append(aln)
            if(len(alns) > 0):
                extend_edge(alns, index_object, query)
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                    yield alnstr
                n_aln += len(alns)
    print('ov asmid =',asmid, 'n_aln =', n_aln)

def get_alns_for_asmid(asmid, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig, debug):
    st = time.time()
    all_rawseqs_object = load_all_rawseqs_object(all_rawseqs_path)
    if(debug == True):
        print('load_all_rawseqs_object', time.time() - st)
    for alnstr in yield_aln_for_asmid(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  

        yield alnstr
        

def mp_get_alns_for_asmid(asmid_queue, alnstr_queue, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):
    all_rawseqs_object = load_all_rawseqs_object(all_rawseqs_path)
    while(True):
        asmid = asmid_queue.get()
        if(isinstance(asmid, str)):
            break
        for alnstr in yield_aln_for_asmid(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  

            alnstr_queue.put(alnstr)
        
    print('Worker stopping')

def mp_write_bam(alnstr_queue, out_path, header, out_path_sorted = ''):
    if(isinstance(header, dict)):
        header = pysam.AlignmentHeader.from_dict(header)
    with pysam.AlignmentFile(out_path, "wb", header=header, threads=16) as outf:
        
        while(True):
            onealn = alnstr_queue.get()
            if(isinstance(onealn, int)):
                break
            aligned = parse_sam_line(onealn, header)
            outf.write(aligned)
    if(out_path_sorted != ''):
        subprocess.run(['samtools', 'sort', '-@4', out_path, '-o', out_path_sorted], check = True)
        subprocess.run(['samtools', 'index', '-@4', out_path_sorted], check = True)

    print('Writing stopping')
def merge_sorted_bam(out_path, bams):
    cmd = ['samtools', 'merge', '-@4', '-f', out_path] + bams
    subprocess.run(cmd, check = True)
    cmd = ['samtools', 'index', '-@4', out_path]
    subprocess.run(cmd, check = True)

@njit
def get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig, maxreadgap = 50, maxrefgap = 50, maxgap = 50, minaln_size = 40):
    alns = []
    alns_anchors = [[(0, 0, 0, 0)]]
    alns_anchors[-1].pop()
    anchor_in_cache = False
    pre_strand = ''
    len_project_arr = len(project_arr)
    for iloc in range(len(project_arr)):
        if((not np.isnan(read_project_arr[iloc])) and (not np.isnan(project_arr[iloc]))):
            q_st = int(read_project_arr[iloc])
            pre_contig = project_arr_contig[iloc]
            pre_readpos = q_st
            if(project_arr[iloc] > 0):
                pre_strand = '+'
                r_st = int(project_arr[iloc])
                pre_refpos = r_st
                anchor = (q_st, q_st + 1, r_st, r_st + 1)
            else:
                pre_strand = '-'
                r_en = 1 - int(project_arr[iloc])
                pre_refpos = r_en - 1
                anchor = (q_st, q_st + 1, r_en - 1, r_en)
            anchor_in_cache = True
            break
    
    if(pre_strand == ''):
        return alns, alns_anchors
    startiloc = iloc + 1
    for iloc in range(startiloc, len_project_arr):
        if((not np.isnan(read_project_arr[iloc])) and (not np.isnan(project_arr[iloc]))):
            contig = project_arr_contig[iloc]
            readpos = int(read_project_arr[iloc])
            
            refpos = int(project_arr[iloc])
            if(refpos > 0):
                strand = '+'
            else:
                strand = '-'
                refpos = -refpos
            if(pre_strand == strand):
                if(contig == pre_contig):
                    in_aln = False
                    readgap = readpos - pre_readpos
                    if(strand == '+'):      
                        refgap = refpos - pre_refpos
                    else:
                        refgap = pre_refpos - refpos
                    if(readgap <= maxreadgap):
                        
                        if(0 < refgap <= maxrefgap):
                            gap = readgap - refgap
                            if(gap < maxgap):
                                in_aln = True
                                
                    elif((readgap <= 100) and (0 < refgap < 10000)):
                        in_aln = True
                    elif((0 < refgap <= 100) and (readgap < 10000)):
                        in_aln = True
                    if(in_aln == True):
                        if(anchor_in_cache == True):
                            if(strand == '+'):
                                if((refpos - pre_refpos) == (readpos - pre_readpos) == 1):
                                    assert anchor[1] == readpos
                                    assert anchor[3] == refpos
                                    anchor = (anchor[0], anchor[1] + 1, anchor[2], anchor[3] + 1)
                                else:
                                    alns_anchors[-1].append(anchor)
                                    anchor = (readpos, readpos + 1, refpos, refpos + 1)
                            else:   
                                if((pre_refpos - refpos) == (readpos - pre_readpos) == 1):
                                    assert anchor[1] == readpos
                                    assert (anchor[2] - 1) == refpos
                                    anchor = (anchor[0], anchor[1] + 1, anchor[2] - 1, anchor[3])
                                else:
                                    alns_anchors[-1].append(anchor)
                                    anchor = (readpos, readpos + 1, refpos, refpos + 1)
                                    
            

                        else:
                            anchor_in_cache = True
                            anchor = (readpos, readpos + 1, refpos, refpos + 1)
                        pre_readpos = readpos
                        pre_refpos = refpos
                        continue
                            
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False               
                        
                    
            if(pre_strand == '+'):
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
            else:
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
            q_st = readpos
            if(strand == '+'):
                r_st = refpos
            else:
                r_en = refpos + 1
            pre_contig, pre_strand, pre_readpos, pre_refpos = contig, strand, readpos, refpos
            if((onepack[3] - onepack[2]) > minaln_size):
                alns.append(onepack)
            else:
                alns_anchors.pop()
            
            anchor = (readpos, readpos + 1, refpos, refpos + 1)
            anchor_in_cache = True
            alns_anchors.append([(0, 0, 0, 0)])
            alns_anchors[-1].pop()

                
        else:
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False
    if(anchor_in_cache == True):
        alns_anchors[-1].append(anchor)
        anchor_in_cache = False        
    if(pre_strand == '+'):
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
    else:
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
    if((onepack[3] - onepack[2]) > minaln_size):
        alns.append(onepack)
    else:
        alns_anchors.pop()
    return alns, alns_anchors


    
def yield_aln_for_asmid_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig, debug = False):  
    def get_third(x):
        return x[2]
    st = time.time()
    n_aln = 0
    try:
        asmseq = asmid2contig[asmid]
        pack = asmid2alns[str(asmid)]
    except:
        return
    project_arr = np.empty(len(asmseq))
    project_arr_contig = np.empty(len(asmseq), dtype = np.int32)
    project_arr.fill(np.nan)
    
    pack.sort(key = get_third)
    overlapping = False
    read_end = 0
    for onepack in pack:
        if(onepack[2] < read_end):
            overlapping = True
            break
        read_end = onepack[3]
        if(onepack[1] == '-'):
            fill_refpos_revse(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        else:
            fill_refpos_forward(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
    if(debug == True):
        print('project_arr', time.time() - st)
    if(overlapping == True):
        for alnstr in yield_aln_for_asmid_overlap_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):
            yield alnstr
        return
    
    if(int(asmid) in asmid2queryname):
        for query_name in asmid2queryname[int(asmid)]:

            query = all_rawseqs_object[query_name].seq
            if(read2contigpos[query_name][2] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False

            aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
            if(isinstance(aln_asm, tuple)):
                read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                alns, alns_anchors = get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig) 
                
                if(len(alns) > 0):
                    extend_edge_anchors(alns, alns_anchors, index_object, query)
                    '''print(query_name)
                    print(alns[0])
                    print(alns_anchors[0])
                    print()'''
                    if(need_reverse == True):
                        onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                    else:
                        onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                    for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                        yield alnstr
                    n_aln += len(alns)

    if(int(asmid) in missed_asmid2queryname):
        for query_name in missed_asmid2queryname[int(asmid)]:
            aln_asm = missed_queryname2alns[query_name][0]
            query = all_rawseqs_object[query_name].seq
            if(aln_asm[1] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False
            read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
            alns, alns_anchors = get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig) 
            
            if(len(alns) > 0):
                extend_edge_anchors(alns, alns_anchors, index_object, query)
                '''print(query_name)
                print(alns[0])
                print(alns_anchors[0])
                print()'''
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                    yield alnstr
                n_aln += len(alns)
    print('asmid =',asmid, 'n_aln =', n_aln)
def yield_aln_for_asmid_overlap_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  
    #print('overlap')
    n_aln = 0
    try:
        asmseq = asmid2contig[asmid]
        pack = asmid2alns[str(asmid)]
    except:
        return
    project_arr_list = []
    project_arr_contig_list = []

    for onepack in pack:
        project_arr = np.empty(len(asmseq))
        project_arr_contig = np.empty(len(asmseq), dtype = np.int32)
        project_arr.fill(np.nan)
        if(onepack[1] == '-'):
            fill_refpos_revse(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        else:
            fill_refpos_forward(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        project_arr_list.append(project_arr)
        project_arr_contig_list.append(project_arr_contig)

    if(int(asmid) in asmid2queryname):
        for query_name in asmid2queryname[int(asmid)]:

            query = all_rawseqs_object[query_name].seq
            if(read2contigpos[query_name][2] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False

            aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
            if(isinstance(aln_asm, tuple)):
                alns, alns_anchors = [], []
                for project_arr, project_arr_contig in zip(project_arr_list, project_arr_contig_list):
                    read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                    tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig)
                    for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):
                        alns.append(aln)
                        alns_anchors.append(aln_anchors)
                if(len(alns) > 0):
                    extend_edge_anchors(alns, alns_anchors, index_object, query)
                    if(need_reverse == True):
                        onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                    else:
                        onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                    for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                        yield alnstr
                    n_aln += len(alns)

    if(int(asmid) in missed_asmid2queryname):
        for query_name in missed_asmid2queryname[int(asmid)]:
            aln_asm = missed_queryname2alns[query_name][0]
            query = all_rawseqs_object[query_name].seq
            if(aln_asm[1] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False
            alns, alns_anchors = [], []
            for project_arr, project_arr_contig in zip(project_arr_list, project_arr_contig_list):
                read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig)
                for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):
                    alns.append(aln)
                    alns_anchors.append(aln_anchors)
            if(len(alns) > 0):
                extend_edge_anchors(alns, alns_anchors, index_object, query)
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                    yield alnstr
                n_aln += len(alns)
    print('ov asmid =',asmid, 'n_aln =', n_aln)
def get_alns_for_asmid_anchors(asmid, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig, debug):
    st = time.time()
    all_rawseqs_object = load_all_rawseqs_object(all_rawseqs_path)
    if(debug == True):
        print('load_all_rawseqs_object', time.time() - st)
    for alnstr in yield_aln_for_asmid_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  

        yield alnstr
def mp_get_alns_for_asmid_anchors(asmid_queue, alnstr_queue, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):
    all_rawseqs_object = load_all_rawseqs_object(all_rawseqs_path)
    while(True):
        asmid = asmid_queue.get()
        if(isinstance(asmid, str)):
            break
        for alnstr in yield_aln_for_asmid_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  

            alnstr_queue.put(alnstr)
        
    print('Worker stopping')
def yield_aln_for_asmid_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig, debug = False):  
    def get_third(x):
        return x[2]
    st = time.time()
    n_aln = 0
    try:
        asmseq = asmid2contig[asmid]
        pack = asmid2alns[str(asmid)]
    except:
        return
    project_arr = np.empty(len(asmseq))
    project_arr_contig = np.empty(len(asmseq), dtype = np.int32)
    project_arr.fill(np.nan)
    
    pack.sort(key = get_third)
    overlapping = False
    read_end = 0
    for onepack in pack:
        if(onepack[2] < read_end):
            overlapping = True
            break
        read_end = onepack[3]
        if(onepack[1] == '-'):
            fill_refpos_revse(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        else:
            fill_refpos_forward(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
    if(debug == True):
        print('project_arr', time.time() - st)
    if(overlapping == True):
        for alnstr in yield_aln_for_asmid_overlap_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):
            yield alnstr
        return
    
    if(int(asmid) in asmid2queryname):
        for query_name in asmid2queryname[int(asmid)]:

            query = all_rawseqs_object[query_name].seq
            if(read2contigpos[query_name][2] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False

            aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
            if(isinstance(aln_asm, tuple)):
                read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                alns, alns_anchors = get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig) 
                
                if(len(alns) > 0):
                    extend_edge_anchors(alns, alns_anchors, index_object, query)
                    '''print(query_name)
                    print(alns[0])
                    print(alns_anchors[0])
                    print()'''
                    if(need_reverse == True):
                        onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                    else:
                        onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                    for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                        yield alnstr
                    n_aln += len(alns)

    if(int(asmid) in missed_asmid2queryname):
        for query_name in missed_asmid2queryname[int(asmid)]:
            aln_asm = missed_queryname2alns[query_name][0]
            query = all_rawseqs_object[query_name].seq
            if(aln_asm[1] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False
            read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
            alns, alns_anchors = get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig) 
            
            if(len(alns) > 0):
                extend_edge_anchors(alns, alns_anchors, index_object, query)
                '''print(query_name)
                print(alns[0])
                print(alns_anchors[0])
                print()'''
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                    yield alnstr
                n_aln += len(alns)
    print('asmid =',asmid, 'n_aln =', n_aln)
def yield_aln_for_asmid_overlap_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  
    #print('overlap')
    n_aln = 0
    try:
        asmseq = asmid2contig[asmid]
        pack = asmid2alns[str(asmid)]
    except:
        return
    project_arr_list = []
    project_arr_contig_list = []

    for onepack in pack:
        project_arr = np.empty(len(asmseq))
        project_arr_contig = np.empty(len(asmseq), dtype = np.int32)
        project_arr.fill(np.nan)
        if(onepack[1] == '-'):
            fill_refpos_revse(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        else:
            fill_refpos_forward(tuple(onepack), project_arr, project_arr_contig, contig2index[(onepack[0])], onepack[6])
        project_arr_list.append(project_arr)
        project_arr_contig_list.append(project_arr_contig)

    if(int(asmid) in asmid2queryname):
        for query_name in asmid2queryname[int(asmid)]:

            query = all_rawseqs_object[query_name].seq
            if(read2contigpos[query_name][2] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False

            aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
            if(isinstance(aln_asm, tuple)):
                alns, alns_anchors = [], []
                for project_arr, project_arr_contig in zip(project_arr_list, project_arr_contig_list):
                    read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                    tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig)
                    for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):
                        alns.append(aln)
                        alns_anchors.append(aln_anchors)
                if(len(alns) > 0):
                    extend_edge_anchors(alns, alns_anchors, index_object, query)
                    if(need_reverse == True):
                        onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                    else:
                        onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                    for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                        yield alnstr
                    n_aln += len(alns)

    if(int(asmid) in missed_asmid2queryname):
        for query_name in missed_asmid2queryname[int(asmid)]:
            aln_asm = missed_queryname2alns[query_name][0]
            query = all_rawseqs_object[query_name].seq
            if(aln_asm[1] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False
            alns, alns_anchors = [], []
            for project_arr, project_arr_contig in zip(project_arr_list, project_arr_contig_list):
                read_project_arr = fill_readpos_forward(tuple(aln_asm), project_arr, aln_asm[6])
                tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig)
                for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):
                    alns.append(aln)
                    alns_anchors.append(aln_anchors)
            if(len(alns) > 0):
                extend_edge_anchors(alns, alns_anchors, index_object, query)
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                    yield alnstr
                n_aln += len(alns)
    print('ov asmid =',asmid, 'n_aln =', n_aln)
def extend_edge_anchors(alns, alns_anchors, mapqs, index_object, query, extend = False):
    if(extend == True):
        iloc = 0
        onealn = alns[iloc]
        extend_st = 0
        onealn = extend_top(onealn, extend_st, index_object, query)
        if(len(alns) > 1):

            for iloc in range(1, len(alns)):
                extend_en = alns[iloc][2]
                if(onealn[3] < extend_en):
                    onealn = extend_tail(onealn, extend_en, index_object, query)
                alns[iloc-1] = onealn

                onealn = alns[iloc]
                extend_st = alns[iloc-1][3]
                if(onealn[2] > extend_st):
                    onealn = extend_top(onealn, extend_st, index_object, query)

        extend_en = len(query)  
        if(onealn[3] < extend_en):
            onealn = extend_tail(onealn, extend_en, index_object, query)
        alns[-1] = onealn
    for iloc in range(len(alns)):
        aln = alns[iloc]
        anchors = alns_anchors[iloc]
        
        try:
            tmpcigar = get_cigar_anchors(index_object, aln, anchors, query)
        except:
            print(aln)
            print(anchors)
            print(query)
            tmpcigar = get_cigar_anchors(index_object, aln, anchors, query)
        assert len(Cigar(tmpcigar)) == (aln[3] - aln[2])
  
        alns[iloc] = (aln[0], aln[1], aln[2], aln[3], aln[4], aln[5], tmpcigar, mapqs[iloc])
@njit
def get_new_alignment_position_anchors(read_project_arr, project_arr, project_arr_contig, index2contig, maxreadgap = 50, maxrefgap = 50, maxgap = 50, minaln_size = 40):
    alns = []
    alns_anchors = [[(0, 0, 0, 0)]]
    alns_anchors[-1].pop()
    anchor_in_cache = False
    pre_strand = ''
    len_project_arr = len(project_arr)
    for iloc in range(len(project_arr)):
        if((not np.isnan(read_project_arr[iloc])) and (not np.isnan(project_arr[iloc]))):
            q_st = int(read_project_arr[iloc])
            pre_contig = project_arr_contig[iloc]
            pre_readpos = q_st
            if(project_arr[iloc] > 0):
                pre_strand = '+'
                r_st = int(project_arr[iloc])
                pre_refpos = r_st
                anchor = (q_st, q_st + 1, r_st, r_st + 1)
            else:
                pre_strand = '-'
                r_en = 1 - int(project_arr[iloc])
                pre_refpos = r_en - 1
                anchor = (q_st, q_st + 1, r_en - 1, r_en)
            anchor_in_cache = True
            break
    
    if(pre_strand == ''):
        return alns, alns_anchors
    startiloc = iloc + 1
    for iloc in range(startiloc, len_project_arr):
        if((not np.isnan(read_project_arr[iloc])) and (not np.isnan(project_arr[iloc]))):
            contig = project_arr_contig[iloc]
            readpos = int(read_project_arr[iloc])
            
            refpos = int(project_arr[iloc])
            if(refpos > 0):
                strand = '+'
            else:
                strand = '-'
                refpos = -refpos
            if(pre_strand == strand):
                if(contig == pre_contig):
                    in_aln = False
                    readgap = readpos - pre_readpos
                    if(strand == '+'):      
                        refgap = refpos - pre_refpos
                    else:
                        refgap = pre_refpos - refpos
                    if(readgap <= maxreadgap):
                        
                        if(0 < refgap <= maxrefgap):
                            gap = readgap - refgap
                            if(gap < maxgap):
                                in_aln = True
                                
                    elif((readgap <= 100) and (0 < refgap < 10000)):
                        in_aln = True
                    elif((0 < refgap <= 100) and (readgap < 10000)):
                        in_aln = True
                    if(in_aln == True):
                        if(anchor_in_cache == True):
                            if(strand == '+'):
                                if((refpos - pre_refpos) == (readpos - pre_readpos) == 1):
                                    assert anchor[1] == readpos
                                    assert anchor[3] == refpos
                                    anchor = (anchor[0], anchor[1] + 1, anchor[2], anchor[3] + 1)
                                else:
                                    alns_anchors[-1].append(anchor)
                                    anchor = (readpos, readpos + 1, refpos, refpos + 1)
                            else:   
                                if((pre_refpos - refpos) == (readpos - pre_readpos) == 1):
                                    assert anchor[1] == readpos
                                    assert (anchor[2] - 1) == refpos
                                    anchor = (anchor[0], anchor[1] + 1, anchor[2] - 1, anchor[3])
                                else:
                                    alns_anchors[-1].append(anchor)
                                    anchor = (readpos, readpos + 1, refpos, refpos + 1)
                                    
            

                        else:
                            anchor_in_cache = True
                            anchor = (readpos, readpos + 1, refpos, refpos + 1)
                        pre_readpos = readpos
                        pre_refpos = refpos
                        continue
                            
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False               
                        
                    
            if(pre_strand == '+'):
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
            else:
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
            q_st = readpos
            if(strand == '+'):
                r_st = refpos
            else:
                r_en = refpos + 1
            pre_contig, pre_strand, pre_readpos, pre_refpos = contig, strand, readpos, refpos
            if((onepack[3] - onepack[2]) > minaln_size):
                alns.append(onepack)
            else:
                alns_anchors.pop()
            
            anchor = (readpos, readpos + 1, refpos, refpos + 1)
            anchor_in_cache = True
            alns_anchors.append([(0, 0, 0, 0)])
            alns_anchors[-1].pop()

                
        else:
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False
    if(anchor_in_cache == True):
        alns_anchors[-1].append(anchor)
        anchor_in_cache = False        
    if(pre_strand == '+'):
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
    else:
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
    if((onepack[3] - onepack[2]) > minaln_size):
        alns.append(onepack)
    else:
        alns_anchors.pop()
    return alns, alns_anchors
def check_asmread(asmseq, read2contigpos, query_name, query):
    largesize = 20000
    r_st = int(read2contigpos[query_name][1])
    r_en = r_st + int(read2contigpos[query_name][5])
    tmp_asmseq = asmseq[r_st: r_en]
    if((len(tmp_asmseq) < largesize) and (len(query) < largesize)):
        info = vacmap_index.k_cigar(tmp_asmseq, query, bw = 30, eqx = True)
        if(len(Cigar(info[0])) == len(query)):
            return (read2contigpos[query_name][0], '+', 0, len(query), r_st, r_en, info[0])

    r_bias_st, r_bias_en, cigar = compute_large_cigar(tmp_asmseq, query)

    if(r_bias_st != -1):

        assert len(query) == len(Cigar(cigar))


        return (read2contigpos[query_name][0], '+', 0, len(query), r_st + r_bias_st, r_st + r_bias_en, cigar)
    else:
        return []
def compute_large_cigar(target, query):
    minimap2_aligner = mappy.Aligner(seq = target)
    for rec in minimap2_aligner.map(query):
        if(rec.q_st > 0):
            if((len(query) - rec.q_en) > 0):
                cigar =  (str(rec.q_st)+'S'+ rec.cigar_str + str(len(query) - rec.q_en) + 'S')
            else:
                cigar =  (str(rec.q_st)+'S'+ rec.cigar_str)
        else:
            if((len(query) - rec.q_en) > 0):
                cigar =  rec.cigar_str + str(len(query) - rec.q_en) + 'S'
            else:
                cigar = rec.cigar_str
        #print(cigar)
        return rec.r_st, rec.r_en, cigar_eqx(target[rec.r_st: rec.r_en].encode(), query.encode(), cigar.encode())
    return -1, -1, ''
def get_cigar_anchors(index_object, aln, anchors, query):
    cigar = []
    contig = aln[0]
    strand = aln[1]
    if strand == '+':
        if anchors[0][0] > aln[2]:
            anchors[0] = (aln[2], anchors[0][1], aln[4], anchors[0][3])
        if aln[3] > anchors[-1][1]:
            anchors[-1] = (anchors[-1][0], aln[3], anchors[-1][2], aln[5])
        pre_read_pos = anchors[0][0]
        pre_ref_pos = anchors[0][2]
        for anchor in anchors:
            assert pre_read_pos <= anchor[0]
            tmpquery = query[pre_read_pos: anchor[0]]
            assert pre_ref_pos <= anchor[2]
            if pre_ref_pos < anchor[2]:
                tmpref = index_object.seq(contig, pre_ref_pos, anchor[2])
            else:
                tmpref = ''
            if len(tmpquery) == 0:
                if len(tmpref) > 0:  # deletion
                    cigar.append(str(len(tmpref)) + 'D')
            else:  # len(tmpquery) > 0
                if len(tmpref) == 0:  # insertion
                    cigar.append(str(len(tmpquery)) + 'I')
                else:
                    if(len(tmpquery) == len(tmpref) and (len(tmpquery) < 2) and (tmpquery != tmpref)):
                        cigar.append(str(len(tmpquery)) + 'X')
                    else:
                        info = vacmap_index.k_cigar(tmpref, tmpquery, bw=-1, zdropvalue=-1, eqx=True)
                        assert len(Cigar(info[0])) == len(tmpquery)
                        cigar.append(info[0])
            cigar.append(str(anchor[1] - anchor[0]) + '=')
            pre_read_pos = anchor[1]
            pre_ref_pos = anchor[3]
        return ''.join(cigar)
    else:
        if anchors[0][0] > aln[2]:
            anchors[0] = (aln[2], anchors[0][1], anchors[0][2], aln[5])
        if aln[3] > anchors[-1][1]:
            anchors[-1] = (anchors[-1][0], aln[3], aln[4], anchors[-1][3])
        pre_read_pos = anchors[0][0]
        pre_ref_pos = anchors[0][3]
        for anchor in anchors:
            assert pre_read_pos <= anchor[0]
            tmpquery = reverse_complement(query[pre_read_pos: anchor[0]])
            assert pre_ref_pos >= anchor[3]
            if pre_ref_pos > anchor[3]:
                tmpref = index_object.seq(contig, anchor[3], pre_ref_pos)
            else:
                tmpref = ''
            if len(tmpquery) == 0:
                if len(tmpref) > 0:  # deletion
                    cigar.append(str(len(tmpref)) + 'D')
            else:  # len(tmpquery) > 0
                if len(tmpref) == 0:  # insertion
                    cigar.append(str(len(tmpquery)) + 'I')
                else:
                    if(len(tmpquery) == len(tmpref) and (len(tmpquery) < 2) and (tmpquery != tmpref)):
                        cigar.append(str(len(tmpquery)) + 'X')
                    else:
                        info = vacmap_index.k_cigar(tmpref, tmpquery, bw=-1, zdropvalue=-1, eqx=True)
                        assert len(Cigar(info[0])) == len(tmpquery)
                        cigar.append(info[0])
            cigar.append(str(anchor[1] - anchor[0]) + '=')
            pre_read_pos = anchor[1]
            pre_ref_pos = anchor[2]
        return ''.join(cigar[::-1])
from Bio.Seq import Seq
def P_alignmentstring(infodict):
                #  0     1     2     3     4     5     6      7      8     9     10
                #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL
    infolist = ['*',   '4',   '*',  '0', '255', '*',  '*',   '0',  '0',   '*',   '*']
    name2iloc = {
        'QNAME': 0,
        'FLAG': 1,
        'RNAME': 2,
        'POS': 3,
        'MAPQ': 4,
        'CIGAR': 5,
        'RNEXT': 6,
        'PNEXT': 7,
        'TLEN': 8,
        'SEQ': 9,
        'QUAL': 10
    }
    def p_other_tag(tag, value):
        if(type(value) == int):
            code = 'i'
        elif(type(value) == float):
            code = 'f'
        elif(type(value) == str):
            code = 'Z'
        else:
            code = 'Z'
        return tag+':'+code+':'+str(value)

    for key in infodict:
        if(key in name2iloc):
            infolist[name2iloc[key]] = infodict[key]
        else:
            infolist.append(p_other_tag(key, infodict[key]))
    return '\t'.join(infolist)
def compute_NM_tag(query, target):
    return edlib.align(query = query, target = target, task = 'distance')['editDistance']
@njit
def mergecigar_(cigarstring):

    numberrange = (ord('0'), ord('9'))
    
    oplist = []
    num = 0
    preop = '0'
    prenum = 0
    for ch in cigarstring:
        c = ord(ch)
        if(numberrange[0] <= c and numberrange[1] >= c):
            num = num * 10 + c - numberrange[0]
        else:
            if(preop == ch):
                prenum = prenum + num
                oplist[-2] = str(prenum)
                num = 0
            else:
                prenum = num
                oplist.append(str(num))
                oplist.append(ch)
                preop = ch
                num = 0
    return oplist
@njit
def compute_NM_tag_eqx(cigarstring):
    nm = 0
    numberrange = (ord('0'), ord('9'))
    num = 0

    for ch in cigarstring:
        c = ord(ch)
        if(numberrange[0] <= c and numberrange[1] >= c):
            num = num * 10 + c - numberrange[0]
        else:
            if(ch in ('X', 'D', 'I')):
                nm += num

            num = 0
    return nm
def mergecigar_n(cigarstring):  
    oplist = mergecigar_(cigarstring)
    n_cigar = len(oplist)
    return ''.join(oplist), n_cigar
def get_bam_dict_str(mapinfo, query, qual):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    option = {'H': False, 'fakecigar': False}
    hardclip = option['H']
    
    rc_query = str(Seq(query).reverse_complement())

    #mapinfo.sort(key = sortbycontig)
    iloc2nm = dict()
    iloc2n_cigar = dict()
    tmpiloc = -1
    fakecigar = option['fakecigar']
    if(fakecigar == True):
        iloc2fakecigar = dict()
        if(hardclip == True):
            clipsyb = 'H'
        else:
            clipsyb = 'S'

    for item in mapinfo:
        item[-1], n_cigar = mergecigar_n(item[-1])
        tmpiloc += 1
        nm = compute_NM_tag_eqx(item[-1])
        iloc2nm[tmpiloc] = nm
        iloc2n_cigar[tmpiloc] = n_cigar
        if(fakecigar == True):
            if(item[3] > 0):
                onefaketop = str(item[3]) + clipsyb
            else:
                onefaketop = ''
            if((len(query) - item[4]) > 0):
                onefaketail = str(len(query) - item[4]) + clipsyb
            else:
                onefaketail = ''
            tmpdiff = item[4] - item[3] - item[6] + item[5]
            if(tmpdiff > 0):
                onefakebody = str(item[6] - item[5]) + 'M' + str(tmpdiff) + 'I'
            elif(tmpdiff < 0):
                onefakebody = str(item[4] - item[3]) + 'M' + str(abs(tmpdiff)) + 'D'
            else:
                onefakebody = str(item[4] - item[3]) + 'M'
            iloc2fakecigar[tmpiloc] =  ''.join((onefaketop, onefakebody, onefaketail))
    
    if((qual != None) and (len(qual) == len(query))):
        query_qualities = qual
        rc_query_qualities = query_qualities[::-1]
    a_list = []
    primary_iloc = 0
    #if(len(mapinfo) > 1):
        #if(mapinfo[0][7] == 1 and mapinfo[1][7] != 1):
            #primary_iloc = 1

    #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL
    for iloc in range(len(mapinfo)):
        
        bam_dict = dict()
        if('rg-id' in option):
            bam_dict['RG'] = option['rg-id']
        primary = mapinfo[iloc]
        bam_dict['QNAME'] = primary[0]
        bam_dict['RNAME'] = primary[1]
        if(iloc == primary_iloc):
            base_value = 0
        else:
            base_value = 2048
        if(primary[2] == '+'):
            bam_dict['FLAG'] = str(base_value)

        else:
            bam_dict['FLAG'] = str(16 + base_value)


        bam_dict['POS'] = str(primary[5] + 1)# SAM Format


        bam_dict['CIGAR'] = primary[8]

        if(len(mapinfo) > 1):
            salist = []
            tmpiloc = -1
            for item in mapinfo:
                tmpiloc += 1
                if(tmpiloc == iloc):
                    continue
                mq = mapinfo[tmpiloc][7]
                if(mq != 0):
                    mq = 60
                else:
                    mq = 1
                nm = iloc2nm[tmpiloc]
                if(fakecigar == False):
                    salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', item[8], ',', str(mq), ',', str(nm)+';')))
                else:
                    salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', iloc2fakecigar[tmpiloc], ',', str(mq), ',', str(nm)+';')))


            bam_dict['SA'] = ''.join(salist)
        mq = mapinfo[iloc][7]
    
        item = primary

        bam_dict['MAPQ'] = str(mq)

        if(item[2] == '+'):
            if(hardclip == False):
                bam_dict['SEQ'] = query
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = query_qualities
            else:
                bam_dict['SEQ'] = query[item[3]: item[4]]
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = query_qualities[item[3]: item[4]]
                

            
        else:
            if(hardclip == False):
                bam_dict['SEQ'] = rc_query
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = rc_query_qualities
            else:
                bam_dict['SEQ'] = rc_query[item[3]: item[4]]
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = rc_query_qualities[item[3]: item[4]]

            
        bam_dict['NM'] = iloc2nm[iloc]


        a_list.append(P_alignmentstring(bam_dict))
    return a_list

def parse_sam_line(sam_line, header):
    """
    Parse a SAM alignment line (str) into a pysam.AlignedSegment object.
    This function is a pure Python replacement for pysam.AlignedSegment.fromstring.

    Args:
        sam_line (str): One alignment line in SAM format.
        header (dict): pysam header dictionary.

    Returns:
        pysam.AlignedSegment
    """
    fields = sam_line.strip().split('\t')
    if len(fields) < 11:
        raise ValueError(f"Invalid SAM line (too few fields): {sam_line}")

    # Core fields
    qname = fields[0]
    flag = int(fields[1])
    rname = fields[2]
    pos = int(fields[3]) - 1  # pysam uses 0-based positions
    mapq = int(fields[4])
    cigar = fields[5]
    rnext = fields[6]
    pnext = int(fields[7]) - 1 if fields[7] != "0" else -1
    tlen = int(fields[8])
    seq = fields[9]
    qual = fields[10]

    seg = pysam.AlignedSegment(header)
    seg.query_name = qname
    seg.flag = flag
    seg.reference_id = seg.header.get_tid(rname) if rname != '*' else -1
    seg.reference_start = pos if rname != '*' else -1
    seg.mapping_quality = mapq
    seg.cigarstring = cigar if cigar != '*' else None
    seg.next_reference_id = seg.header.get_tid(rnext) if rnext != '*' else -1
    seg.next_reference_start = pnext if rnext != '*' else -1
    seg.template_length = tlen
    seg.query_sequence = seq if seq != '*' else None
    seg.query_qualities = pysam.qualitystring_to_array(qual) if qual != '*' else None

    # Optional fields
    for opt in fields[11:]:
        if ":" in opt:
            tag, type_, value = opt.split(":", 2)
            # Try to convert value to int/float if possible
            if type_ == 'i':
                value = int(value)
            elif type_ == 'f':
                value = float(value)
            seg.set_tag(tag, value, value_type=type_)

    return seg
def local_asm_onefile(fastq_path_queue, asm_method, datatype, lock, threads, large_size = 10000):
    if(asm_method == 'hifiasm'):
        while(True):
            fastq_path, fastq_size = fastq_path_queue.get()
            if(isinstance(fastq_path, int)):
                break
            hifipath = '/home/dhy/hifiasm/hifiasm'
            if(fastq_size > large_size):
                lock.acquire()
                if(datatype == 'hifi'):
                    a = subprocess.run([hifipath, '-o', fastq_path[:-5]+'tmp', '-t'+str(threads), '-f0', fastq_path], capture_output=True)
                else:
                    a = subprocess.run([hifipath, '--ont', '-o', fastq_path[:-5]+'tmp', '-t'+str(threads), '-f0', fastq_path], capture_output=True)

                lock.release()
            else:

                if(datatype == 'hifi'):
                    a = subprocess.run([hifipath, '-o', fastq_path[:-5]+'tmp', '-t'+str(threads), '-f0', fastq_path], capture_output=True)
                else:
                    a = subprocess.run([hifipath, '--ont', '-o', fastq_path[:-5]+'tmp', '-t'+str(threads), '-f0', fastq_path], capture_output=True)

            if(a.returncode == 0):
                #print(a)
                removefile(fastq_path[:-5]+'tmp*', exclude = 'bp.p_utg.gfa')
                readset = get_readnames(fastq_path)
                gfa_path = fastq_path[:-5]+'tmp.bp.p_utg.gfa'
                tmp_read2contigpos, tmp_asmid2contig, tmp_asmid2queryname = load_gfa(gfa_path)
                drop_readset = set(tmp_read2contigpos.keys())
                readset = list(readset - drop_readset)
                selected_fastq_path = fastq_path[:-5]+'selected.fastq'
                write_selected_fastq(fastq_path, selected_fastq_path, readset)

                asm_fasta_path = fastq_path[:-5]+'asm.fasta'
                write_fasta_from_dict(tmp_asmid2contig, asm_fasta_path)
                selected_sam_path = selected_fastq_path[:-5] + 'sam'
                if(fastq_size > large_size):
                    lock.acquire()
                    a = subprocess.run(['minimap2', '-aYx', 'map-ont', '-t', str(threads), '--eqx', asm_fasta_path, selected_fastq_path, '-o', selected_sam_path], capture_output=True)
                    lock.release()
                else:
                    a = subprocess.run(['minimap2', '-aYx', 'map-ont', '-t', str(1), '--eqx', asm_fasta_path, selected_fastq_path, '-o', selected_sam_path], capture_output=True)

                if(a.returncode == 0):
                    continue

            removefile(fastq_path[:-5]+'tmp*')
            print(a)
            continue
    elif(asm_method == "flye"):
        while(True):
            fastq_path, fastq_size = fastq_path_queue.get()
            if(isinstance(fastq_path, int)):
                break
            flyepath = 'flye'
            lock.acquire()
            lock.release()
            if(fastq_size > large_size):
                lock.acquire()

                a = subprocess.run([flyepath, '--nano-raw', fastq_path, '--keep-haplotypes', '-o', fastq_path[:-5]+'tmp', '-t', str(threads)], capture_output=True)

                lock.release()
            else:

                a = subprocess.run([flyepath, '--nano-raw', fastq_path, '--keep-haplotypes', '-o', fastq_path[:-5]+'tmp', '-t', str(1)], capture_output=True)
            
            if(a.returncode == 0):
                #print(a)
                #removefile(fastq_path[:-5]+'tmp*', exclude = 'bp.p_utg.gfa')
                readset = get_readnames(fastq_path)
                gfa_path = fastq_path[:-5]+'tmp/assembly_graph.gfa'
                a = subprocess.run(['cp', gfa_path, fastq_path[:-5]+'tmp.bp.p_utg.gfa'], capture_output=True)
                if(a.returncode == 0):
                    tmp_read2contigpos, tmp_asmid2contig, tmp_asmid2queryname = load_gfa(gfa_path)
                    drop_readset = set(tmp_read2contigpos.keys())
                    readset = list(readset - drop_readset)
                    selected_fastq_path = fastq_path[:-5]+'selected.fastq'
                    write_selected_fastq(fastq_path, selected_fastq_path, readset)

                    asm_fasta_path = fastq_path[:-5]+'asm.fasta'
                    write_fasta_from_dict(tmp_asmid2contig, asm_fasta_path)
                    selected_sam_path = selected_fastq_path[:-5] + 'sam'
                    lock.acquire()
                    lock.release()
                    if(fastq_size > large_size):
                        lock.acquire()
                        a = subprocess.run(['minimap2', '-aYx', 'map-ont', '-t', str(threads), '--eqx', asm_fasta_path, selected_fastq_path, '-o', selected_sam_path], capture_output=True)
                        lock.release()
                    else:
                        a = subprocess.run(['minimap2', '-aYx', 'map-ont', '-t', str(1), '--eqx', asm_fasta_path, selected_fastq_path, '-o', selected_sam_path], capture_output=True)
                    removefile(fastq_path[:-5]+'tmp')
                    if(a.returncode == 0):

                        continue

            removefile(fastq_path[:-5]+'tmp')
            print(a)
            continue
    else:
        
        while(True):
            fastq_path, fastq_size = fastq_path_queue.get()
            if(isinstance(fastq_path, int)):
                break
            wtdbg2path = '/home/dhy/wtdbg2/wtdbg2.pl'
            lock.acquire()
            lock.release()
            if(fastq_size > large_size):
                lock.acquire()

                a = subprocess.run([wtdbg2path, '-x', 'ont',  '-o', 
                                    fastq_path[:-5]+'tmp', '-t', str(threads),
                                   fastq_path], capture_output=True)

                lock.release()
            else:

                a = subprocess.run([wtdbg2path, '-x', 'ont',  '-o', 
                                    fastq_path[:-5]+'tmp', '-t', str(threads),
                                   fastq_path], capture_output=True)
            
            if(a.returncode == 0):
                #print(a)
                removefile(fastq_path[:-5]+'tmp*', exclude = 'cns.fa')
                readset = get_readnames(fastq_path)
                gfa_path = fastq_path[:-5]+'tmp.cns.fa'
               
                tmp_read2contigpos, tmp_asmid2contig, tmp_asmid2queryname = load_gfa(gfa_path, isfasta = True)
                drop_readset = set(tmp_read2contigpos.keys())
                readset = list(readset - drop_readset)
                selected_fastq_path = fastq_path[:-5]+'selected.fastq'
                write_selected_fastq(fastq_path, selected_fastq_path, readset)

                asm_fasta_path = fastq_path[:-5]+'asm.fasta'
                write_fasta_from_dict(tmp_asmid2contig, asm_fasta_path)
                selected_sam_path = selected_fastq_path[:-5] + 'sam'
                lock.acquire()
                lock.release()
                if(fastq_size > large_size):
                    lock.acquire()
                    a = subprocess.run(['minimap2', '-aYx', 'map-ont', '-t', str(threads), '--eqx', asm_fasta_path, selected_fastq_path, '-o', selected_sam_path], capture_output=True)
                    lock.release()
                else:
                    a = subprocess.run(['minimap2', '-aYx', 'map-ont', '-t', str(1), '--eqx', asm_fasta_path, selected_fastq_path, '-o', selected_sam_path], capture_output=True)

                if(a.returncode == 0):

                    continue

            removefile(fastq_path[:-5]+'tmp*')
            print(a)
            continue
@njit
def check_largeindel(cigarstring, size):
    nm = 0
    numberrange = (ord('0'), ord('9'))
    num = 0

    for ch in cigarstring:
        c = ord(ch)
        if(numberrange[0] <= c and numberrange[1] >= c):
            num = num * 10 + c - numberrange[0]
        else:
            if(ch in ('D', 'I')):
                if(num >= size):
                    return True
            num = 0
    return False
@njit
def fill_refpos_forward_bias(onepack, contig_index, cigar):
    
    project_arr = np.empty(onepack[3] - onepack[2])
    project_arr_contig = np.empty(onepack[3] - onepack[2], dtype = np.int32)
    project_arr.fill(np.nan)
    
    refpos = onepack[4]
    readpos = 0
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):

                project_arr[readpos - onepack[2]: readpos + num - onepack[2]] = np.arange(refpos, refpos + num)
                project_arr_contig[readpos - onepack[2]: readpos + num - onepack[2]] = contig_index

            if(c in use_query):
                readpos += num
            if(c in use_ref):
                refpos += num
            num = 0
    return onepack[2], project_arr, project_arr_contig

@njit
def fill_refpos_revse_bias(len_asmseq, onepack, contig_index, cigar):
    
    project_arr = np.empty(onepack[3] - onepack[2])
    project_arr_contig = np.empty(onepack[3] - onepack[2], dtype = np.int32)
    project_arr.fill(np.nan)
    
    refpos = onepack[4]
    readpos = len_asmseq
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):
                project_arr[readpos - num - onepack[2]: readpos - onepack[2]] = np.arange(-refpos - num + 1, -refpos + 1)
                project_arr_contig[readpos - num - onepack[2]: readpos - onepack[2]] = contig_index
            if(c in use_query):
                readpos -= num
            if(c in use_ref):
                refpos += num
            num = 0 
    return onepack[2], project_arr, project_arr_contig

@njit
def fill_readpos_forward_bias(onepack, cigar):
    
    read_project_arr = np.empty(onepack[5] - onepack[4])
    read_project_arr.fill(np.nan)
    

    refpos = onepack[4]
    readpos = onepack[2]
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):
                #print(refpos, refpos + num, readpos, readpos + num)
                read_project_arr[refpos - onepack[4]: refpos + num - onepack[4]] = np.arange(readpos, readpos + num)
            if(c in use_query):
                readpos += num
            if(c in use_ref):
                refpos += num
            num = 0
    return onepack[4], read_project_arr
@njit
def get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, project_arr_contig, index2contig, maxreadgap = 50, maxrefgap = 50, maxgap = 50, minaln_size = 40):
    alns = []
    alns_anchors = [[(0, 0, 0, 0)]]
    alns_anchors[-1].pop()
    anchor_in_cache = False
    pre_strand = ''
    len_project_arr = min(read_slice_st + len(read_project_arr), ref_slice_st + len(project_arr))
    for iloc in range(max(read_slice_st, ref_slice_st), len_project_arr):
        if((not np.isnan(read_project_arr[iloc - read_slice_st])) and (not np.isnan(project_arr[iloc - ref_slice_st]))):
            q_st = int(read_project_arr[iloc - read_slice_st])
            pre_contig = project_arr_contig[iloc - ref_slice_st]
            pre_readpos = q_st
            if(project_arr[iloc - ref_slice_st] > 0):
                pre_strand = '+'
                r_st = int(project_arr[iloc - ref_slice_st])
                pre_refpos = r_st
                anchor = (q_st, q_st + 1, r_st, r_st + 1)
            else:
                pre_strand = '-'
                r_en = 1 - int(project_arr[iloc - ref_slice_st])
                pre_refpos = r_en - 1
                anchor = (q_st, q_st + 1, r_en - 1, r_en)
            anchor_in_cache = True
            break
    
    if(pre_strand == ''):
        return alns, alns_anchors
    startiloc = iloc + 1
    for iloc in range(startiloc, len_project_arr):
        if((not np.isnan(read_project_arr[iloc - read_slice_st])) and (not np.isnan(project_arr[iloc - ref_slice_st]))):
            contig = project_arr_contig[iloc - ref_slice_st]
            readpos = int(read_project_arr[iloc - read_slice_st])
            
            refpos = int(project_arr[iloc - ref_slice_st])
            if(refpos > 0):
                strand = '+'
            else:
                strand = '-'
                refpos = -refpos
            if(pre_strand == strand):
                if(contig == pre_contig):
                    in_aln = False
                    readgap = readpos - pre_readpos
                    if(strand == '+'):      
                        refgap = refpos - pre_refpos
                    else:
                        refgap = pre_refpos - refpos
                    if(readgap <= maxreadgap):
                        
                        if(0 < refgap <= maxrefgap):
                            gap = readgap - refgap
                            if(gap < maxgap):
                                in_aln = True
                                
                    elif((readgap <= 100) and (0 < refgap < 10000)):
                        in_aln = True
                    elif((0 < refgap <= 100) and (readgap < 10000)):
                        in_aln = True
                    if(in_aln == True):
                        if(anchor_in_cache == True):
                            if(strand == '+'):
                                if((refpos - pre_refpos) == (readpos - pre_readpos) == 1):
                                    assert anchor[1] == readpos
                                    assert anchor[3] == refpos
                                    anchor = (anchor[0], anchor[1] + 1, anchor[2], anchor[3] + 1)
                                else:
                                    alns_anchors[-1].append(anchor)
                                    anchor = (readpos, readpos + 1, refpos, refpos + 1)
                            else:   
                                if((pre_refpos - refpos) == (readpos - pre_readpos) == 1):
                                    assert anchor[1] == readpos
                                    assert (anchor[2] - 1) == refpos
                                    anchor = (anchor[0], anchor[1] + 1, anchor[2] - 1, anchor[3])
                                else:
                                    alns_anchors[-1].append(anchor)
                                    anchor = (readpos, readpos + 1, refpos, refpos + 1)
                                    
            

                        else:
                            anchor_in_cache = True
                            anchor = (readpos, readpos + 1, refpos, refpos + 1)
                        pre_readpos = readpos
                        pre_refpos = refpos
                        continue
                            
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False               
                        
                    
            if(pre_strand == '+'):
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
            else:
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
            q_st = readpos
            if(strand == '+'):
                r_st = refpos
            else:
                r_en = refpos + 1
            pre_contig, pre_strand, pre_readpos, pre_refpos = contig, strand, readpos, refpos
            if((onepack[3] - onepack[2]) > minaln_size):
                alns.append(onepack)
            else:
                alns_anchors.pop()
            
            anchor = (readpos, readpos + 1, refpos, refpos + 1)
            anchor_in_cache = True
            alns_anchors.append([(0, 0, 0, 0)])
            alns_anchors[-1].pop()

                
        else:
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False
    if(anchor_in_cache == True):
        alns_anchors[-1].append(anchor)
        anchor_in_cache = False        
    if(pre_strand == '+'):
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
    else:
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
    if((onepack[3] - onepack[2]) > minaln_size):
        alns.append(onepack)
    else:
        alns_anchors.pop()
    return alns, alns_anchors
import bisect


def yield_aln_for_asmid_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  
    #print('overlap')
    def get_overlaped_chains1(read_slice_st, read_project_arr, chains):
        for ref_slice_st, project_arr, project_arr_contig in chains:
            max_st = max(read_slice_st, ref_slice_st)
            min_en = min(read_slice_st + len(read_project_arr), ref_slice_st + len(project_arr))
            if(min_en > max_st):
                yield ref_slice_st, project_arr, project_arr_contig
        
    n_aln = 0
    try:
        asmseq = asmid2contig[asmid]
        pack = asmid2alns[str(asmid)]
    except:
        return
    chains = []


    for onepack in pack:

        if(onepack[1] == '-'):
            slice_st, project_arr, project_arr_contig = fill_refpos_revse_bias(len(asmseq), tuple(onepack), contig2index[(onepack[0])], onepack[6])
        else:
            slice_st, project_arr, project_arr_contig = fill_refpos_forward_bias(tuple(onepack), contig2index[(onepack[0])], onepack[6])
        chains.append((slice_st, project_arr, project_arr_contig))
        
    max_chain_len = max(len(p) for _, p, _ in chains) if chains else 0
    
    if(int(asmid) in asmid2queryname):
        for query_name in asmid2queryname[int(asmid)]:

            query = all_rawseqs_object[query_name].seq
            if(read2contigpos[query_name][2] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False

            aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
            if(isinstance(aln_asm, tuple)):
                alns, alns_anchors = [], []
                
                read_slice_st, read_project_arr = fill_readpos_forward_bias(tuple(aln_asm), aln_asm[6])
                for ref_slice_st, project_arr, project_arr_contig in get_overlaped_chains(read_slice_st, read_project_arr, chains, max_chain_len):
                    tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, project_arr_contig, index2contig)

                    for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):
                        alns.append(aln)
                        alns_anchors.append(aln_anchors)
                if(len(alns) > 0):
                    extend_edge_anchors(alns, alns_anchors, index_object, query)
                    if(need_reverse == True):
                        onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                    else:
                        onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                    for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                        yield alnstr
                    n_aln += len(alns)

    if(int(asmid) in missed_asmid2queryname):
        for query_name in missed_asmid2queryname[int(asmid)]:
            aln_asm = missed_queryname2alns[query_name][0]
            query = all_rawseqs_object[query_name].seq
            if(aln_asm[1] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False
            alns, alns_anchors = [], []
            read_slice_st, read_project_arr = fill_readpos_forward_bias(tuple(aln_asm), aln_asm[6])
            for ref_slice_st, project_arr, project_arr_contig in get_overlaped_chains(read_slice_st, read_project_arr, chains, max_chain_len):
                tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, project_arr_contig, index2contig)

                for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):
                    alns.append(aln)
                    alns_anchors.append(aln_anchors)
            if(len(alns) > 0):
                extend_edge_anchors(alns, alns_anchors, index_object, query)
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                    yield alnstr
                n_aln += len(alns)
    print('ov asmid =',asmid, 'n_aln =', n_aln)

@njit
def fill_refpos_forward_bias(onepack, contig_index, cigar):
    max_int32 = np.iinfo(np.int32).max
    project_arr = np.empty(onepack[3] - onepack[2], dtype = np.int32)
    project_arr_contig = np.empty(onepack[3] - onepack[2], dtype = np.uint16)
    project_arr.fill(max_int32)
    
    refpos = onepack[4]
    readpos = 0
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):

                project_arr[readpos - onepack[2]: readpos + num - onepack[2]] = np.arange(refpos, refpos + num)
                project_arr_contig[readpos - onepack[2]: readpos + num - onepack[2]] = contig_index

            if(c in use_query):
                readpos += num
            if(c in use_ref):
                refpos += num
            num = 0
    return onepack[2], project_arr, project_arr_contig

@njit
def fill_refpos_revse_bias(len_asmseq, onepack, contig_index, cigar):
    

    
    max_int32 = np.iinfo(np.int32).max
    project_arr = np.empty(onepack[3] - onepack[2], dtype = np.int32)
    project_arr_contig = np.empty(onepack[3] - onepack[2], dtype = np.uint16)
    project_arr.fill(max_int32)
    
    refpos = onepack[4]
    readpos = len_asmseq
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):
                project_arr[readpos - num - onepack[2]: readpos - onepack[2]] = np.arange(-refpos - num + 1, -refpos + 1)
                project_arr_contig[readpos - num - onepack[2]: readpos - onepack[2]] = contig_index
            if(c in use_query):
                readpos -= num
            if(c in use_ref):
                refpos += num
            num = 0 
    return onepack[2], project_arr, project_arr_contig

@njit
def fill_readpos_forward_bias(onepack, cigar):
    
    max_int32 = np.iinfo(np.int32).max
    read_project_arr = np.empty(onepack[5] - onepack[4], dtype = np.int32)
    read_project_arr.fill(max_int32)
    

    refpos = onepack[4]
    readpos = onepack[2]
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):
                #print(refpos, refpos + num, readpos, readpos + num)
                read_project_arr[refpos - onepack[4]: refpos + num - onepack[4]] = np.arange(readpos, readpos + num)
            if(c in use_query):
                readpos += num
            if(c in use_ref):
                refpos += num
            num = 0
    return onepack[4], read_project_arr
@njit
def get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, project_arr_contig, index2contig, maxreadgap = 50, maxrefgap = 50, maxgap = 50, minaln_size = 40):
    max_int32 = np.iinfo(np.int32).max
    alns = []
    alns_anchors = [[(0, 0, 0, 0)]]
    alns_anchors[-1].pop()
    anchor_in_cache = False
    pre_strand = ''
    len_project_arr = min(read_slice_st + len(read_project_arr), ref_slice_st + len(project_arr))
    for iloc in range(max(read_slice_st, ref_slice_st), len_project_arr):
        if((read_project_arr[iloc - read_slice_st] != max_int32) and (project_arr[iloc - ref_slice_st] != max_int32)):
            q_st = (read_project_arr[iloc - read_slice_st])
            pre_contig = project_arr_contig[iloc - ref_slice_st]
            pre_readpos = q_st
            if(project_arr[iloc - ref_slice_st] > 0):
                pre_strand = '+'
                r_st = (project_arr[iloc - ref_slice_st])
                pre_refpos = r_st
                anchor = (q_st, q_st + 1, r_st, r_st + 1)
            else:
                pre_strand = '-'
                r_en = 1 - (project_arr[iloc - ref_slice_st])
                pre_refpos = r_en - 1
                anchor = (q_st, q_st + 1, r_en - 1, r_en)
            anchor_in_cache = True
            break
    
    if(pre_strand == ''):
        return alns, alns_anchors
    startiloc = iloc + 1
    for iloc in range(startiloc, len_project_arr):
        if((read_project_arr[iloc - read_slice_st] != max_int32) and (project_arr[iloc - ref_slice_st] != max_int32)):
            contig = project_arr_contig[iloc - ref_slice_st]
            readpos = (read_project_arr[iloc - read_slice_st])
            
            refpos = int(project_arr[iloc - ref_slice_st])
            if(refpos > 0):
                strand = '+'
            else:
                strand = '-'
                refpos = -refpos
            if(pre_strand == strand):
                if(contig == pre_contig):
                    in_aln = False
                    readgap = readpos - pre_readpos
                    if(strand == '+'):      
                        refgap = refpos - pre_refpos
                    else:
                        refgap = pre_refpos - refpos
                    
                                
                    if((readgap <= 100) and (0 < refgap < 10000)):
                        in_aln = True
                    elif((0 < refgap <= 100) and (readgap < 10000)):
                        in_aln = True
                    if(in_aln == True):
                        if(anchor_in_cache == True):
                            if(strand == '+'):
                                if((refpos - pre_refpos) == (readpos - pre_readpos) == 1):
                                    assert anchor[1] == readpos
                                    assert anchor[3] == refpos
                                    anchor = (anchor[0], anchor[1] + 1, anchor[2], anchor[3] + 1)
                                else:
                                    alns_anchors[-1].append(anchor)
                                    anchor = (readpos, readpos + 1, refpos, refpos + 1)
                            else:   
                                if((pre_refpos - refpos) == (readpos - pre_readpos) == 1):
                                    assert anchor[1] == readpos
                                    assert (anchor[2] - 1) == refpos
                                    anchor = (anchor[0], anchor[1] + 1, anchor[2] - 1, anchor[3])
                                else:
                                    alns_anchors[-1].append(anchor)
                                    anchor = (readpos, readpos + 1, refpos, refpos + 1)
                                    
            

                        else:
                            anchor_in_cache = True
                            anchor = (readpos, readpos + 1, refpos, refpos + 1)
                        pre_readpos = readpos
                        pre_refpos = refpos
                        continue
                            
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False               
                        
                    
            if(pre_strand == '+'):
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
            else:
                onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
            q_st = readpos
            if(strand == '+'):
                r_st = refpos
            else:
                r_en = refpos + 1
            pre_contig, pre_strand, pre_readpos, pre_refpos = contig, strand, readpos, refpos
            if((onepack[3] - onepack[2]) > minaln_size):
                alns.append(onepack)
            else:
                alns_anchors.pop()
            
            anchor = (readpos, readpos + 1, refpos, refpos + 1)
            anchor_in_cache = True
            alns_anchors.append([(0, 0, 0, 0)])
            alns_anchors[-1].pop()

                
        else:
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False
    if(anchor_in_cache == True):
        alns_anchors[-1].append(anchor)
        anchor_in_cache = False        
    if(pre_strand == '+'):
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
    else:
        onepack = (index2contig[pre_contig], pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
    if((onepack[3] - onepack[2]) > minaln_size):
        alns.append(onepack)
    else:
        alns_anchors.pop()
    return alns, alns_anchors
import bisect
def getfirst(x):
    return x[0]
def get_overlaped_chains(read_slice_st, read_project_arr, chains, max_chain_len):
    """
    Binary search version to yield overlapped chains.
    Assumes chains is a list of tuples (ref_slice_st, project_arr, project_arr_contig),
    sorted by ref_slice_st.
    """
    if not chains:
        return

    qs = read_slice_st
    qe = read_slice_st + len(read_project_arr)
    
    # Compute max chain length for conservative left bound
    #max_chain_len = max(len(p) for _, p, _ in chains) if chains else 0
    
    # Key function for bisect: ref_slice_st
    def key(c):
        return c[0]
    
    # Left: first chain where start >= qs - max_chain_len
    left_key = qs - max_chain_len
    left = bisect.bisect_left(chains, left_key, key=key)
    
    # Right: first chain where start >= qe
    right = bisect.bisect_left(chains, qe, key=key)
    
    # Check candidates in [left, right)
    for i in range(max(0, left), right):
        ref_slice_st, project_arr, project_arr_contig, mapq = chains[i]
        ce_st = ref_slice_st
        ce_en = ref_slice_st + len(project_arr)
        max_st = max(qs, ce_st)
        min_en = min(qe, ce_en)
        if min_en > max_st:
            yield ref_slice_st, project_arr, project_arr_contig, mapq
def yield_aln_for_asmid_anchors(asmid, all_rawseqs_object, asmid2contig, asmid2alns, contig2index, asmid2queryname, all_rawseqs_path, read2contigpos, missed_asmid2queryname, missed_queryname2alns, index_object, index2contig):  
    #print('overlap')
    def get_overlaped_chains1(read_slice_st, read_project_arr, chains):
        for ref_slice_st, project_arr, project_arr_contig in chains:
            max_st = max(read_slice_st, ref_slice_st)
            min_en = min(read_slice_st + len(read_project_arr), ref_slice_st + len(project_arr))
            if(min_en > max_st):
                yield ref_slice_st, project_arr, project_arr_contig
        
    n_aln = 0
    try:
        asmseq = asmid2contig[asmid]
        pack = asmid2alns[str(asmid)]
    except:
        return
    chains = []


    for onepack in pack:

        if(onepack[1] == '-'):
            slice_st, project_arr, project_arr_contig = fill_refpos_revse_bias(len(asmseq), tuple(onepack), contig2index[(onepack[0])], onepack[6])
        else:
            slice_st, project_arr, project_arr_contig = fill_refpos_forward_bias(tuple(onepack), contig2index[(onepack[0])], onepack[6])
        chains.append((slice_st, project_arr, project_arr_contig))
    chains.sort(key = getfirst)    
    max_chain_len = max(len(p) for _, p, _ in chains) if chains else 0
    
    if(int(asmid) in asmid2queryname):
        for query_name in asmid2queryname[int(asmid)]:

            query = all_rawseqs_object[query_name].seq
            if(read2contigpos[query_name][2] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False

            aln_asm = check_asmread(asmseq, read2contigpos, query_name, query)
            if(isinstance(aln_asm, tuple)):
                alns, alns_anchors = [], []
                
                read_slice_st, read_project_arr = fill_readpos_forward_bias(tuple(aln_asm), aln_asm[6])
                for ref_slice_st, project_arr, project_arr_contig in get_overlaped_chains(read_slice_st, read_project_arr, chains, max_chain_len):
                    tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, project_arr_contig, index2contig)

                    for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):
                        alns.append(aln)
                        alns_anchors.append(aln_anchors)
                if(len(alns) > 0):
                    extend_edge_anchors(alns, alns_anchors, index_object, query)
                    if(need_reverse == True):
                        onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                    else:
                        onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                    for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                        yield alnstr
                    n_aln += len(alns)

    if(int(asmid) in missed_asmid2queryname):
        for query_name in missed_asmid2queryname[int(asmid)]:
            aln_asm = missed_queryname2alns[query_name][0]
            query = all_rawseqs_object[query_name].seq
            if(aln_asm[1] == '-'):
                query = reverse_complement(query)
                need_reverse = True
            else:
                need_reverse = False
            alns, alns_anchors = [], []
            read_slice_st, read_project_arr = fill_readpos_forward_bias(tuple(aln_asm), aln_asm[6])
            for ref_slice_st, project_arr, project_arr_contig in get_overlaped_chains(read_slice_st, read_project_arr, chains, max_chain_len):
                tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, project_arr_contig, index2contig)

                for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):
                    alns.append(aln)
                    alns_anchors.append(aln_anchors)
            if(len(alns) > 0):
                extend_edge_anchors(alns, alns_anchors, index_object, query)
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                for alnstr in yield_alnstr(onemapinfo, all_rawseqs_object, query_name):
                    yield alnstr
                n_aln += len(alns)
    print('ov asmid =',asmid, 'n_aln =', n_aln)
def compute_one_chain(aln_queue, cooked_queue, contig2index, convert_eqx):
    size = 0
    count = 0
    asmid2chains = dict()
    while(True):
        contig_length, (contig, strand, cigarstring, reference_start, query_sequence, query_name, mapq) = aln_queue.get()
        if(contig_length):
            aln = list(decode_cigar(contig, strand, contig_length, cigarstring, reference_start, truereadpos = True))

            if(convert_eqx == True):

                cigarstring = cigar_eqx(index_object.seq(contig, aln[4], aln[5]).encode(), query_sequence.encode(), cigarstring.encode())

            cigarstring = cigarstring.replace('H', 'S')

            aln.append(cigarstring)
            
            onepack = aln
            if(onepack[1] == '-'):
                slice_st, project_arr, project_arr_contig = fill_refpos_revse_bias(contig_length, tuple(onepack), contig2index[(onepack[0])], onepack[6])
            else:
                slice_st, project_arr, project_arr_contig = fill_refpos_forward_bias(tuple(onepack), contig2index[(onepack[0])], onepack[6])

            try:
                asmid2chains[query_name].append((slice_st, project_arr, project_arr_contig, mapq))
            except:
                asmid2chains[query_name] = [(slice_st, project_arr, project_arr_contig, mapq)]

            size += len(project_arr)
            count += 1
            if(count % 100 == 0):
                print(count, size)

        else:
            break
    cooked_queue.put(asmid2chains)
    print('compute_one_chain Stoping')

def yield_alnstr_dirct(onemapinfo, query, qual, queryname):
    mapinfo = []

    for line in onemapinfo:
        contig, strand, q_st, q_en, r_st, r_en, cigar, mapq = line
        if(q_st > 0):
            cigar = str(q_st) + 'S' + cigar
        if((len(query) - q_en) > 0):
            cigar += str(len(query) - q_en) + 'S'
        mapinfo.append([queryname, contig, strand, q_st, q_en, r_st, r_en, mapq, cigar])
        if(len(Cigar(cigar)) != len(query)):
            print(len(query), len(Cigar(cigar)))
            print(queryname)
            print(contig, strand, q_st, q_en, r_st, r_en)
            print(cigar)
            print()

    for onealn in get_bam_dict_str(mapinfo, query, qual):
        yield onealn

def qualities_to_fastq_str(qualities):
    """
    Convert list of Phred qualities (ints) to FASTQ string.
    """
    return "".join(chr(q + 33) for q in qualities)

#mark
def lift_one_aln(raw_queue, cooked_queue, asmid2chains, asmid2max_chain_len, index2contig, index_object, missedpath, nofilter = False):
    string_cache = []
    string_cache_max = 1000
    if(nofilter == False):
        fq = open(missedpath, 'w')
    while(True):
        data_rev = raw_queue.get()
        if(data_rev is None):
            break
        for aln_asm, t_query, t_qual, query_name, asmid, query in data_rev:
            if(asmid not in asmid2chains):
                continue
            if(nofilter == False and check_largeindel(aln_asm[6], size = 30) == True):
                t_qual = qualities_to_fastq_str(t_qual) 
                fq.write(f"@{query_name}\n")
                fq.write(f"{t_query}\n")
                fq.write("+\n")
                fq.write(f"{t_qual}\n")
                continue
            t_qual = qualities_to_fastq_str(t_qual) 
            if(aln_asm[1] == '-'):

                need_reverse = True
            else:
                need_reverse = False
            alns, alns_anchors, mapqs = [], [], []
            #print(query_name)
            read_slice_st, read_project_arr = fill_readpos_forward_bias(tuple(aln_asm), aln_asm[6])
            for ref_slice_st, project_arr, project_arr_contig, mapq in get_overlaped_chains(read_slice_st, read_project_arr, asmid2chains[asmid], asmid2max_chain_len[asmid]):
                #print('hit')
                tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, project_arr_contig, index2contig)
                #print(project_arr)
                #print(read_project_arr)
                #print(tmp_alns)
                #print()

                for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):

                    alns.append(aln)
                    alns_anchors.append(aln_anchors)
                    mapqs.append(mapq)
            if(len(alns) > 0):
                extend_edge_anchors(alns, alns_anchors, mapqs, index_object, query)
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                #for line in onemapinfo:
                    #print(line)
                    #if(line[1] == '+'):
                        #print((t_query)[line[2]: line[2]+20])
                    #else:
                        #print(reverse_complement(t_query)[line[2]: line[2]+20])
                    #print(index_object.seq(line[0], line[4], line[4]+20))
                    #print()
                #print(len(onemapinfo))

                for alnstr in yield_alnstr_dirct(onemapinfo, t_query, t_qual, query_name):
                    if(len(string_cache) > string_cache_max):
                        alnstring = '\n'.join(string_cache)
                        alnstring = alnstring.encode('utf-8') + b'\n'
                        cooked_queue.put(alnstring)
                        string_cache = []
                    string_cache.append(alnstr)
                #print(tmpcount)
                #print()
    if(len(string_cache) > 0):
        alnstring = '\n'.join(string_cache)
        alnstring = alnstring.encode('utf-8') + b'\n'
        cooked_queue.put(alnstring)
    if(nofilter == False):
        fq.close()
    #break
def mp_write_bam_sorted(alnstr_queue, header, out_path_sorted):
    # Create output dir if needed
    #os.makedirs(os.path.dirname(out_path_sorted), exist_ok=True)
    
    if isinstance(header, dict):
        header = str(pysam.AlignmentHeader.from_dict(header))
    assert isinstance(header, str), "Header must be str or dict"
    
    proc = subprocess.Popen(
        ['samtools', 'sort', '-@', '8', '-', '-o', out_path_sorted, '-O', 'BAM'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False
    )
    
    try:
        proc.stdin.write(header.rstrip().encode('utf-8') + b'\n')
        
        while(True):
            onealn = alnstr_queue.get()
            if(isinstance(onealn, int)):
                break
            proc.stdin.write(onealn)
        
        # NO manual close here -- communicate() will flush and close stdin
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, proc.args, stdout, stderr)
        
        # Index only if sort succeeds
        #subprocess.run(['samtools', 'index', '-@', '4', out_path_sorted], check=True)
    
    
    except Exception as e:
        proc.terminate()
        raise RuntimeError(f"Error in mp_write_bam_sorted: {e}")
    
    finally:
        if proc.stdin and not proc.stdin.closed:
            proc.stdin.close()

def get_alns_for_lift(bamfilepath, missedpath, nofilter = False):
    if(bamfilepath.endswith('sam')):
        bamfile = pysam.AlignmentFile(bamfilepath,'r', threads = 16)
    else:
        bamfile = pysam.AlignmentFile(bamfilepath,'rb', threads = 16)

    if(nofilter == False):
        countunmap = 0
        usedset = set()
        with open(missedpath, 'w') as fq:
            for AlignedSegment in bamfile:
                
                if((AlignedSegment.is_secondary == True) or (AlignedSegment.is_supplementary == True)):
                    continue
                query_name = AlignedSegment.query_name

                query = AlignedSegment.query_sequence
                if(AlignedSegment.is_mapped == False):
                    countunmap += 1

                    qual = qualities_to_fastq_str(AlignedSegment.query_qualities)
                    fq.write(f"@{query_name}\n")
                    fq.write(f"{query}\n")
                    fq.write("+\n")
                    fq.write(f"{qual}\n")
                asmid = AlignedSegment.reference_name
                strand = decode_strand_from_flag(AlignedSegment.flag)
                read_length = AlignedSegment.query_length
                aln_asm = [asmid, strand, 0, read_length, AlignedSegment.reference_start, AlignedSegment.reference_end, AlignedSegment.cigarstring]


                t_query = AlignedSegment.get_forward_sequence()


                if(strand == '+'):
                    t_qual = AlignedSegment.query_qualities
                else:
                    t_qual = AlignedSegment.query_qualities[::-1]
                if((AlignedSegment.has_tag('SA') == True) or ((read_length - AlignedSegment.query_alignment_length) > 100)):
                    t_qual = qualities_to_fastq_str(t_qual) 
                    fq.write(f"@{query_name}\n")
                    fq.write(f"{t_query}\n")
                    fq.write("+\n")
                    fq.write(f"{t_qual}\n")
                else:
                    yield aln_asm, t_query, t_qual, query_name, asmid, query
                #break
        print('countunmap =', countunmap)
    else:
        for AlignedSegment in bamfile:
            if((AlignedSegment.is_secondary == True) or (AlignedSegment.is_mapped == False)):
                continue


            asmid = AlignedSegment.reference_name
            strand = decode_strand_from_flag(AlignedSegment.flag)
            read_length = AlignedSegment.query_length
            aln_asm = [asmid, strand, 0, read_length, AlignedSegment.reference_start, AlignedSegment.reference_end, AlignedSegment.cigarstring]

            query = AlignedSegment.query_sequence
            query_name = AlignedSegment.query_name
            t_query = AlignedSegment.get_forward_sequence()


            if(strand == '+'):
                t_qual = AlignedSegment.query_qualities
            else:
                t_qual = AlignedSegment.query_qualities[::-1]
            
            yield aln_asm, t_query, t_qual, query_name, asmid, query
            #break
@njit
def fill_refpos_forward_bias(onepack, cigar):
    max_int32 = np.iinfo(np.int32).max
    project_arr = np.empty(onepack[3] - onepack[2], dtype = np.int32)
    project_arr.fill(max_int32)
    
    refpos = onepack[4]
    readpos = 0
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):

                project_arr[readpos - onepack[2]: readpos + num - onepack[2]] = np.arange(refpos, refpos + num)

            if(c in use_query):
                readpos += num
            if(c in use_ref):
                refpos += num
            num = 0
    return onepack[2], project_arr

@njit
def fill_refpos_revse_bias(len_asmseq, onepack, cigar):
    

    
    max_int32 = np.iinfo(np.int32).max
    project_arr = np.empty(onepack[3] - onepack[2], dtype = np.int32)
    project_arr.fill(max_int32)
    
    refpos = onepack[4]
    readpos = len_asmseq
    num = 0
    zero = ord('0')
    nine = ord('9')
    use_query = ord('M'), ord('I'), ord('S'), ord('X'), ord('=')
    inssyb = ord('I')
    delsyb = ord('D')
    use_ref = ord('M'), ord('D'), ord('X'), ord('=')
    match = ord('M'), ord('X'), ord('=')
    equal = ord('=')
    for ch in cigar:
        c = ord(ch)
        if(c >= zero and nine >= c):
            num = num * 10 + c - zero
        else:         
            if(c == equal):
                project_arr[readpos - num - onepack[2]: readpos - onepack[2]] = np.arange(-refpos - num + 1, -refpos + 1)
            if(c in use_query):
                readpos -= num
            if(c in use_ref):
                refpos += num
            num = 0 
    return onepack[2], project_arr
@njit
def get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, contig, maxreadgap = 50, maxrefgap = 50, maxgap = 50, minaln_size = 40):
    max_int32 = np.iinfo(np.int32).max
    alns = []
    alns_anchors = [[(0, 0, 0, 0)]]
    alns_anchors[-1].pop()
    anchor_in_cache = False
    pre_strand = ''
    len_project_arr = min(read_slice_st + len(read_project_arr), ref_slice_st + len(project_arr))
    for iloc in range(max(read_slice_st, ref_slice_st), len_project_arr):
        if((read_project_arr[iloc - read_slice_st] != max_int32) and (project_arr[iloc - ref_slice_st] != max_int32)):
            q_st = (read_project_arr[iloc - read_slice_st])
            
            pre_readpos = q_st
            if(project_arr[iloc - ref_slice_st] > 0):
                pre_strand = '+'
                r_st = (project_arr[iloc - ref_slice_st])
                pre_refpos = r_st
                anchor = (q_st, q_st + 1, r_st, r_st + 1)
            else:
                pre_strand = '-'
                r_en = 1 - (project_arr[iloc - ref_slice_st])
                pre_refpos = r_en - 1
                anchor = (q_st, q_st + 1, r_en - 1, r_en)
            anchor_in_cache = True
            break
    
    if(pre_strand == ''):
        return alns, alns_anchors
    startiloc = iloc + 1
    for iloc in range(startiloc, len_project_arr):
        if((read_project_arr[iloc - read_slice_st] != max_int32) and (project_arr[iloc - ref_slice_st] != max_int32)):
            
            readpos = (read_project_arr[iloc - read_slice_st])
            
            refpos = int(project_arr[iloc - ref_slice_st])
            if(refpos > 0):
                strand = '+'
            else:
                strand = '-'
                refpos = -refpos
            if(pre_strand == strand):
      
                in_aln = False
                readgap = readpos - pre_readpos
                if(strand == '+'):      
                    refgap = refpos - pre_refpos
                else:
                    refgap = pre_refpos - refpos


                if((readgap <= 100) and (0 < refgap < 10000)):
                    in_aln = True
                elif((0 < refgap <= 100) and (readgap < 10000)):
                    in_aln = True
                if(in_aln == True):
                    if(anchor_in_cache == True):
                        if(strand == '+'):
                            if((refpos - pre_refpos) == (readpos - pre_readpos) == 1):
                                assert anchor[1] == readpos
                                assert anchor[3] == refpos
                                anchor = (anchor[0], anchor[1] + 1, anchor[2], anchor[3] + 1)
                            else:
                                alns_anchors[-1].append(anchor)
                                anchor = (readpos, readpos + 1, refpos, refpos + 1)
                        else:   
                            if((pre_refpos - refpos) == (readpos - pre_readpos) == 1):
                                assert anchor[1] == readpos
                                assert (anchor[2] - 1) == refpos
                                anchor = (anchor[0], anchor[1] + 1, anchor[2] - 1, anchor[3])
                            else:
                                alns_anchors[-1].append(anchor)
                                anchor = (readpos, readpos + 1, refpos, refpos + 1)



                    else:
                        anchor_in_cache = True
                        anchor = (readpos, readpos + 1, refpos, refpos + 1)
                    pre_readpos = readpos
                    pre_refpos = refpos
                    continue
                            
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False               
                        
                    
            if(pre_strand == '+'):
                onepack = (contig, pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
            else:
                onepack = (contig, pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
            q_st = readpos
            if(strand == '+'):
                r_st = refpos
            else:
                r_en = refpos + 1
            pre_strand, pre_readpos, pre_refpos =  strand, readpos, refpos
            if((onepack[3] - onepack[2]) > minaln_size):
                alns.append(onepack)
            else:
                alns_anchors.pop()
            
            anchor = (readpos, readpos + 1, refpos, refpos + 1)
            anchor_in_cache = True
            alns_anchors.append([(0, 0, 0, 0)])
            alns_anchors[-1].pop()

                
        else:
            if(anchor_in_cache == True):
                alns_anchors[-1].append(anchor)
                anchor_in_cache = False
    if(anchor_in_cache == True):
        alns_anchors[-1].append(anchor)
        anchor_in_cache = False        
    if(pre_strand == '+'):
        onepack = (contig, pre_strand, q_st, pre_readpos + 1, r_st, pre_refpos + 1)
    else:
        onepack = (contig, pre_strand, q_st, pre_readpos + 1, pre_refpos, r_en)
    if((onepack[3] - onepack[2]) > minaln_size):
        alns.append(onepack)
    else:
        alns_anchors.pop()
    return alns, alns_anchors
def compute_one_chain(aln_queue, cooked_queue, contig2index, convert_eqx):
    size = 0
    count = 0
    asmid2chains = dict()
    while(True):
        contig_length, (contig, strand, cigarstring, reference_start, query_sequence, query_name, mapq) = aln_queue.get()
        if(contig_length):
            aln = list(decode_cigar(contig, strand, contig_length, cigarstring, reference_start, truereadpos = True))

            if(convert_eqx == True):

                cigarstring = cigar_eqx(index_object.seq(contig, aln[4], aln[5]).encode(), query_sequence.encode(), cigarstring.encode())

            cigarstring = cigarstring.replace('H', 'S')

            aln.append(cigarstring)
            
            onepack = aln
            if(onepack[1] == '-'):
                slice_st, project_arr = fill_refpos_revse_bias(contig_length, tuple(onepack), onepack[6])
            else:
                slice_st, project_arr = fill_refpos_forward_bias(tuple(onepack), onepack[6])

            try:
                asmid2chains[query_name].append((slice_st, project_arr, contig, mapq))
            except:
                asmid2chains[query_name] = [(slice_st, project_arr, contig, mapq)]

            size += len(project_arr)
            count += 1
            #if(count % 100 == 0):
                #print(count, size)

        else:
            break
    cooked_queue.put(asmid2chains)
    #print('compute_one_chain Stoping')
def lift_one_aln(raw_queue, cooked_queue, asmid2chains, asmid2max_chain_len, index2contig, index_object, missedpath, nofilter = False):
    string_cache = []
    string_cache_max = 1000
    if(nofilter == False):
        fq = open(missedpath, 'w')
    while(True):
        data_rev = raw_queue.get()
        if(data_rev is None):
            break
        for aln_asm, t_query, t_qual, query_name, asmid, query in data_rev:
            if(asmid not in asmid2chains):
                continue
            if(nofilter == False and check_largeindel(aln_asm[6], size = 30) == True):
                t_qual = qualities_to_fastq_str(t_qual) 
                fq.write(f"@{query_name}\n")
                fq.write(f"{t_query}\n")
                fq.write("+\n")
                fq.write(f"{t_qual}\n")
                continue
            t_qual = qualities_to_fastq_str(t_qual) 
            if(aln_asm[1] == '-'):

                need_reverse = True
            else:
                need_reverse = False
            alns, alns_anchors, mapqs = [], [], []
            #print(query_name)
            read_slice_st, read_project_arr = fill_readpos_forward_bias(tuple(aln_asm), aln_asm[6])
            for ref_slice_st, project_arr, contig, mapq in get_overlaped_chains(read_slice_st, read_project_arr, asmid2chains[asmid], asmid2max_chain_len[asmid]):
                #print('hit')
                tmp_alns, tmp_alns_anchors = get_new_alignment_position_anchors_bias(read_slice_st, ref_slice_st, read_project_arr, project_arr, contig)
                #print(project_arr)
                #print(read_project_arr)
                #print(tmp_alns)
                #print()

                for aln, aln_anchors in zip(tmp_alns, tmp_alns_anchors):

                    alns.append(aln)
                    alns_anchors.append(aln_anchors)
                    mapqs.append(mapq)
            if(len(alns) > 0):
                extend_edge_anchors(alns, alns_anchors, mapqs, index_object, query)
                if(need_reverse == True):
                    onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                else:
                    onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
                #for line in onemapinfo:
                    #print(line)
                    #if(line[1] == '+'):
                        #print((t_query)[line[2]: line[2]+20])
                    #else:
                        #print(reverse_complement(t_query)[line[2]: line[2]+20])
                    #print(index_object.seq(line[0], line[4], line[4]+20))
                    #print()
                #print(len(onemapinfo))

                for alnstr in yield_alnstr_dirct(onemapinfo, t_query, t_qual, query_name):
                    if(len(string_cache) > string_cache_max):
                        alnstring = '\n'.join(string_cache)
                        alnstring = alnstring.encode('utf-8') + b'\n'
                        cooked_queue.put(alnstring)
                        string_cache = []
                    string_cache.append(alnstr)
                #print(tmpcount)
                #print()
    if(len(string_cache) > 0):
        alnstring = '\n'.join(string_cache)
        alnstring = alnstring.encode('utf-8') + b'\n'
        cooked_queue.put(alnstring)
    if(nofilter == False):
        fq.close()
import numba
import multiprocessing
import time

def remap_func(asm2ref_path, read2asm_path, workdir, threads, datatype, refpath, index_object, header, nofilter = False, realign_unlifted = True):
    def run_minimap_to_bam(minimap_args, output_path):
        minimap_cmd = aligner_env + minimap_args
        samtools_cmd = aligner_env + ['samtools', 'view', '-@', '8', '-b', '-o', output_path, '-']
        with subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE) as minimap_proc:
            try:
                subprocess.run(samtools_cmd, stdin=minimap_proc.stdout, check=True)
            finally:
                if minimap_proc.stdout is not None:
                    minimap_proc.stdout.close()
                minimap_return = minimap_proc.wait()
        if minimap_return != 0:
            raise subprocess.CalledProcessError(minimap_return, minimap_cmd)
        return subprocess.CompletedProcess(minimap_cmd, 0)

    aligner_env_name = 'minimap2_env'
    aligner_env = ["conda", "run", "-n", aligner_env_name]
    aligner_env = []

    minimap_presets = {
        'ont': 'map-ont',
        'hifi': 'map-hifi',
        'hqlr': 'lr:hq',
    }

    if asm2ref_path != '':
        print('Build Chain')
        convert_eqx = False
        bamfilepath = asm2ref_path
        if bamfilepath.endswith('sam'):
            bamfile = pysam.AlignmentFile(bamfilepath, 'r', check_sq=False)
        else:
            bamfile = pysam.AlignmentFile(bamfilepath, 'rb', check_sq=False)

        contig2index = dict()
        index2contig = bamfile.header.references
        iloc = -1
        for contigname, contiglength in zip(bamfile.references, bamfile.lengths):
            iloc += 1
            contig2index[contigname] = iloc

        index2contig = numba.typed.List(index2contig)

        aln_queue, cooked_queue = multiprocessing.Manager().Queue(maxsize=200), multiprocessing.Manager().Queue()
        workers = []
        for iloc in range(threads):
            worker = multiprocessing.Process(target=compute_one_chain, args=(aln_queue, cooked_queue, contig2index, convert_eqx))
            worker.start()
            workers.append(worker)
        count = 0
        st_time = time.time()
        for AlignedSegment in bamfile:
            if AlignedSegment.is_secondary:
                continue
            if not AlignedSegment.is_mapped:
                continue

            contig = AlignedSegment.reference_name
            strand = decode_strand_from_flag(AlignedSegment.flag)

            read_length = AlignedSegment.infer_read_length()
            cigarstring = AlignedSegment.cigarstring

            tmpdata = (
                contig,
                strand,
                cigarstring,
                AlignedSegment.reference_start,
                AlignedSegment.query_sequence,
                AlignedSegment.query_name,
                AlignedSegment.mapping_quality,
            )

            aln_queue.put((read_length, tmpdata))

            count += 1
            if count % 10000 == 0:
                print(time.time() - st_time)
                st_time = time.time()

        for iloc in range(threads):
            aln_queue.put((None, (0, 0, 0, 0, 0, 0, 0)))

        for worker in workers:
            worker.join()
        asmid2chains = dict()
        for iloc in range(threads):
            tmp_asmid2chains = cooked_queue.get()
            for asmid in tmp_asmid2chains:
                if asmid in asmid2chains:
                    asmid2chains[asmid] += tmp_asmid2chains[asmid]
                else:
                    asmid2chains[asmid] = tmp_asmid2chains[asmid]
        asmid2max_chain_len = dict()
        for asmid in asmid2chains:
            chains = asmid2chains[asmid]
            max_chain_len = max(len(p) for _, p, _, _ in chains) if chains else 0
            asmid2max_chain_len[asmid] = max_chain_len
            asmid2chains[asmid].sort(key=getfirst)

        print('Build Chain completed')
        print('Liftover')
        bamfilepath = read2asm_path
        st_time = time.time()

        sorted_lifed_read_alns_path = workdir + 'read2ref.liftover.sorted.bam'

        missedreadprefix = workdir + 'filtered_read/'
        try:
            os.mkdir(missedreadprefix)
        except:
            shutil.rmtree(missedreadprefix)
            os.mkdir(missedreadprefix)
        missediloc = 0

        raw_queue, cooked_queue = multiprocessing.Manager().Queue(maxsize=200), multiprocessing.Manager().Queue(maxsize=200)
        workers = []
        for iloc in range(threads):
            missediloc += 1
            worker = multiprocessing.Process(
                target=lift_one_aln,
                args=(
                    raw_queue,
                    cooked_queue,
                    asmid2chains,
                    asmid2max_chain_len,
                    index2contig,
                    index_object,
                    missedreadprefix + str(missediloc) + '.fastq',
                    nofilter,
                ),
            )
            worker.start()
            workers.append(worker)
        write_process = multiprocessing.Process(
            target=mp_write_bam_sorted,
            args=(cooked_queue, header, sorted_lifed_read_alns_path)
        )
        write_process.start()
        cache = []
        cache_size = 100
        missediloc += 1
        for datatosend in get_alns_for_lift(bamfilepath, missedreadprefix + str(missediloc) + '.fastq', nofilter):
            if len(cache) > cache_size:
                raw_queue.put(cache)
                cache = []
            cache.append(datatosend)
        if len(cache) > 0:
            raw_queue.put(cache)
        for work_p in workers:
            raw_queue.put(None)
        for work_p in workers:
            work_p.join()
        cooked_queue.put(0)
        write_process.join()
        time.time() - st_time
        print('Remapping complete')
        if realign_unlifted:
            missed_reads_fastqs = glob.glob(missedreadprefix + '*.fastq')
            if not missed_reads_fastqs:
                print('No missed reads to realign.')
                return sorted_lifed_read_alns_path, ''
            if datatype not in minimap_presets:
                print('not support datatype ' + datatype)
                return
            missedread2ref_unsorted = workdir + 'missedread2ref.unsorted.bam'
            missedread2refpath_bam = workdir + 'missedread2ref.sorted.bam'
            minimap_args = [
                'minimap2',
                '-I', '1000G',
                '-aYx', minimap_presets[datatype],
                '-t', str(threads),
                '--eqx',
                refpath,
            ] + missed_reads_fastqs
            run_minimap_to_bam(minimap_args, missedread2ref_unsorted)

            sort_cmd = aligner_env + ['samtools', 'sort', '-@', str(threads), missedread2ref_unsorted, '-o', missedread2refpath_bam]
            sort_proc = subprocess.run(sort_cmd, capture_output=True)
            if sort_proc.returncode != 0:
                print(sort_proc)
            else:
                if os.path.exists(missedread2ref_unsorted):
                    os.remove(missedread2ref_unsorted)
                return sorted_lifed_read_alns_path, missedread2refpath_bam
        else:
            return sorted_lifed_read_alns_path, missedreadprefix + '*.fastq'

    else:
        inputreadpath = read2asm_path
        if datatype not in minimap_presets:
            print('not support datatype ' + datatype)
            return
        input_reads = glob.glob(inputreadpath)
        if not input_reads:
            print('No input reads found for remapping.')
            return
        missedread2ref_unsorted = workdir + 'missedread2ref.unsorted.bam'
        missedread2ref_bam = workdir + 'missedread2ref.bam'
        minimap_args = [
            'minimap2',
            '-I', '1000G',
            '-aYx', minimap_presets[datatype],
            '-t', str(threads),
            '--eqx',
            refpath,
        ] + input_reads
        run_minimap_to_bam(minimap_args, missedread2ref_unsorted)

        sort_cmd = aligner_env + ['samtools', 'sort', '-@', str(threads), missedread2ref_unsorted, '-o', missedread2ref_bam]
        sort_proc = subprocess.run(sort_cmd, capture_output=True)
        if sort_proc.returncode != 0:
            print(sort_proc)
        else:
            if os.path.exists(missedread2ref_unsorted):
                os.remove(missedread2ref_unsorted)
            return '', missedread2ref_bam

    raise RuntimeError("Remapping failed: minimap2 alignment did not complete successfully.")


def get_liftover_data(inputreadpath, refpath, workdir, datatype, threads, asmdata='', self_align=True, rutg = True):
    assert datatype in ('ont', 'hqlr', 'hifi')
    aligner_env_name = 'minimap2_env'
    aligner_env = ["conda", "run", "-n", aligner_env_name]
    aligner_env = []
    asmfasta = asmdata
    def run_minimap_to_bam(minimap_args, output_path):
        minimap_cmd = aligner_env + minimap_args
        samtools_cmd = aligner_env + ['samtools', 'view', '-@', '8', '-b', '-o', output_path, '-']
        with subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE) as minimap_proc:
            try:
                subprocess.run(samtools_cmd, stdin=minimap_proc.stdout, check=True)
            finally:
                if minimap_proc.stdout is not None:
                    minimap_proc.stdout.close()
                minimap_return = minimap_proc.wait()
        if minimap_return != 0:
            raise subprocess.CalledProcessError(minimap_return, minimap_cmd)
        return subprocess.CompletedProcess(minimap_cmd, 0)

    if not os.path.isdir(workdir):
        os.mkdir(workdir)
    if not os.path.isdir(workdir + 'ASM_data'):
        os.mkdir(workdir + 'ASM_data')
    asmout = workdir + 'ASM_data/asmdata'

    if (asmfasta == '') and (not os.path.isfile(asmout + '.bp.r_utg.gfa')):
        print('Start assembly')
        hifipath = 'hifiasm'
        if datatype == 'hqlr':
            subprocess.run([hifipath, '--ont', '-o', asmout, '-t' + str(threads)] + glob.glob(inputreadpath), check=True)
        elif datatype == 'hifi':
            subprocess.run([hifipath, '-o', asmout, '-t' + str(threads)] + glob.glob(inputreadpath), check=True)
        else:
            assert datatype in ('hqlr', 'hifi')
    


    if asmfasta == '':
        if(rutg == True):
            asmfasta = asmout + '.bp.r_utg.fasta'
            count = gfa_to_fasta(asmout + '.bp.r_utg.gfa', asmfasta)
        else:#asmdata.bp.hap1.p_ctg.gfa
            asmfasta = asmout + '.bp.hap12.p_ctg.fasta'
            count = gfa_to_fasta([asmout + '.bp.hap1.p_ctg.gfa', asmout + '.bp.hap2.p_ctg.gfa'], asmfasta)
        if count == 0:
            asmfasta = ''
    elif asmfasta[-3:] == 'gfa':
        print('Assembly provided')
        count = gfa_to_fasta(asmfasta, asmout + '.bp.r_utg.fasta')
        asmfasta = asmout + '.bp.provided.fasta'
        if count == 0:
            asmfasta = ''
    if asmfasta != '':
        asm2ref_path = workdir + 'asm2ref.bam'
        read2asm_path = workdir + 'read2asm.bam'
        print('Aligning assembly to reference')
        a = run_minimap_to_bam(
            ['minimap2', '-I', '1000G', '-ax', 'asm5', '--cs', '-r2k', '-t', str(threads), '--eqx', refpath, asmfasta],
            asm2ref_path
        )
        if a.returncode == 0:
            print('Aligning read to assembly')
            input_reads = glob.glob(inputreadpath)
            if datatype == 'ont':
                a = run_minimap_to_bam(
                    ['minimap2', '-I', '1000G', '-aYx', 'map-ont', '-t', str(threads), '--eqx', asmfasta] + input_reads,
                    read2asm_path
                )
            elif datatype in ('hifi',):
                if self_align:
                    a = run_minimap_to_bam(
                        ['minimap2', '-I', '1000G', '-aYx', 'lr:hqae', '-t', str(threads), '--eqx', asmfasta] + input_reads,
                        read2asm_path
                    )
                else:
                    print('Liftover')
                    a = run_minimap_to_bam(
                        ['minimap2', '-I', '1000G', '-aYx', 'map-hifi', '-t', str(threads), '--eqx', asmfasta] + input_reads,
                        read2asm_path
                    )
            elif datatype in ('hqlr',):
                if self_align:
                    a = run_minimap_to_bam(
                        ['minimap2', '-I', '1000G', '-aYx', 'lr:hqae', '-t', str(threads), '--eqx', asmfasta] + input_reads,
                        read2asm_path
                    )
                else:
                    print('Liftover')
                    a = run_minimap_to_bam(
                        ['minimap2', '-I', '1000G', '-aYx', 'lr:hq', '-t', str(threads), '--eqx', asmfasta] + input_reads,
                        read2asm_path
                    )
    else:
        return '', inputreadpath


    if a.returncode == 0:
        return asm2ref_path, read2asm_path
    else:
        print(a)
        assert a.returncode == 0

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional, Tuple




def samtools_merge(
    output: Path,
    inputs: list[str],
    threads: int,
) -> None:
    cmd = [
        "samtools",
        "merge",
        "-f",
        "-@",
        str(threads),
        "-o",
        str(output),
    ] + inputs
    subprocess.run(cmd, check=True)


def samtools_index(bam_path: Path, threads: int) -> None:
    cmd = [
        "samtools",
        "index",
        "-@",
        str(threads),
        str(bam_path),
    ]
    subprocess.run(cmd, check=True)


def run_pipeline(args: argparse.Namespace) -> None:
    asmdata = args.asmdata or ""
    if(args.workdir[-1] != '/'):
        args.workdir += '/'
    # Determine hybrid mode from flags.
    if args.hybrid is None:
        hybrid = False
    else:
        hybrid = args.hybrid
        if asmdata == "":
            raise RuntimeError("--hybrid requires --asmdata.")

    if asmdata and not Path(asmdata).exists():
        raise FileNotFoundError(f"{asmdata} not exist.")

    if args.datatype == "ont" and not asmdata:
        raise RuntimeError("Old ONT data requires --asmdata.")

    rutg = not args.germline
    self_align = not hybrid
    realign_unlifted = not hybrid

    asm2ref_path, read2asm_path = get_liftover_data(
        args.inputreadpath,
        args.refpath,
        args.workdir,
        args.datatype,
        args.threads,
        asmdata=asmdata if hybrid else "",
        rutg=rutg,
        self_align=self_align,
    )

    index_object = vacmap_index.Aligner(args.refpath, w=20, k=15)
    if not index_object:
        raise RuntimeError("ERROR: failed to load index")

    header = {"HD": {"VN": "1.0"}, "SQ": []}
    for name_bytes, *_ in index_object.seq_offset:
        seq_name = name_bytes.decode()
        header["SQ"].append(
            {"LN": len(index_object.seq(seq_name)), "SN": seq_name}
        )
    #asm2ref_path, read2asm_path, workdir, threads, datatype, refpath, nofilter = False, realign_unlifted = True
    firstbam, missedreadpath = remap_func(
        asm2ref_path,
        read2asm_path,
        args.workdir,
        args.threads,
        args.datatype,
        args.refpath,
        index_object,
        header,
        nofilter=False,
        realign_unlifted=realign_unlifted
    )

    secondbam = ""
    if hybrid:
        remap_workdir = Path(missedreadpath).resolve().parent / "REMAP_data"
        remap_workdir.mkdir(parents=True, exist_ok=True)

        asm2ref_path, read2asm_path = get_liftover_data(
            missedreadpath,
            args.refpath,
            str(remap_workdir),
            args.datatype,
            args.threads,
            asmdata="",
            rutg=rutg,
            self_align=True,
        )
        secondbam, missedreadpath = remap_func(
            asm2ref_path,
            read2asm_path,
            str(remap_workdir),
            args.threads,
            args.datatype,
            args.refpath,
            index_object,
            header,
            nofilter=False,
            realign_unlifted=True,
        )

    if not firstbam:
        print("No Realigned BAM generated.", file=sys.stderr)
        return

    if not missedreadpath.endswith("bam"):
        raise RuntimeError("Error: missed read path is not a BAM file.")

    samtools_threads = min(args.threads, 8)
    lifted_bam = Path(f"{args.outputprefix}lifted.bam").resolve()
    all_bam = Path(f"{args.outputprefix}all.bam").resolve()
    lifted_bam.parent.mkdir(parents=True, exist_ok=True)

    if secondbam:
        print("Merging Realigned BAMs...", flush=True)
        samtools_merge(
            lifted_bam,
            [firstbam, secondbam],
            samtools_threads,
        )
        samtools_index(lifted_bam, samtools_threads)

        print("Merging Realigned and missed BAM...", flush=True)
        samtools_merge(
            all_bam,
            [str(lifted_bam), missedreadpath],
            samtools_threads,
        )
        samtools_index(all_bam, samtools_threads)
    else:
        print("Renaming Realigned BAM...", flush=True)
        Path(firstbam).replace(lifted_bam)
        samtools_index(lifted_bam, samtools_threads)

        print("Merging Realigned and missed BAM...", flush=True)
        samtools_merge(
            all_bam,
            [str(lifted_bam), missedreadpath],
            samtools_threads,
        )
        samtools_index(all_bam, samtools_threads)

    print("Pipeline completed successfully.", flush=True)


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Liftover pipeline command-line tool."
    )
    parser.add_argument("--inputreadpath", required=True, help="Glob for read files.")
    parser.add_argument("--refpath", required=True, help="Reference FASTA.")
    parser.add_argument("--workdir", required=True, help="Working directory.")
    parser.add_argument(
        "--datatype",
        required=True,
        choices=("ont", "hqlr", "hifi"),
        help="Read type.",
    )
    parser.add_argument(
        "--threads",
        required=True,
        type=int,
        help="Number of threads.",
    )
    parser.add_argument(
        "--asmdata",
        default="",
        help="Optional assembly file (FASTA or GFA).",
    )
    parser.add_argument(
        "--outputprefix",
        required=True,
        help="Prefix for output BAMs, e.g. /path/to/sample_.",
    )
    parser.add_argument(
        "--germline",
        action="store_true",
        help="Use if working with germline samples.",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--hybrid",
        dest="hybrid",
        action="store_true",
        help="Force hybrid mode.",
    )
    group.add_argument(
        "--non-hybrid",
        dest="hybrid",
        action="store_false",
        help="Force non-hybrid mode.",
    )
    parser.set_defaults(hybrid=None)
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    try:
        args = parse_args(argv)
        run_pipeline(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
