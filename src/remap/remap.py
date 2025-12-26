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
                    


from Bio.Seq import Seq


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())
from collections import defaultdict
from typing import Any, List, Set

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



def extend_edge_anchors(alns, alns_anchors, mapqs, index_object, query, minanchor, extend = False):
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
        

        tmpcigar = get_cigar_anchors(index_object, aln, anchors, query, minanchor)

  
        alns[iloc] = (aln[0], aln[1], aln[2], aln[3], aln[4], aln[5], tmpcigar, mapqs[iloc])
@njit
def checkindel_count(cigar: bytes) -> int:
    DEL = ord('D')
    INS = ord('I')

    total = 0
    num = 0
    
    hasdel = 0
    hasins = 0
    
    for i in range(len(cigar)):
        b = cigar[i]
        if 48 <= b <= 57:  # '0'..'9'
            num = num * 10 + (b - 48)
        else:
            if(b == DEL and num > 5):
                hasdel += 1
            elif(b == INS and num > 5):
                hasins += 1
            num = 0
    return hasdel, hasins
    
def get_cigar_anchors(index_object, aln, anchors, query, minanchor = 0):
    cigar = []
    contig = aln[0]
    strand = aln[1]
    if strand == '+':
        if anchors[0][0] > aln[2]:
            anchors[0] = (aln[2], anchors[0][1], aln[4], anchors[0][3])
        if aln[3] > anchors[-1][1]:
            anchors[-1] = (anchors[-1][0], aln[3], anchors[-1][2], aln[5])
        iloc = 0
        indelsize = 30
        while((iloc + 2) < len(anchors)):
            midsize = anchors[iloc + 1][1] - anchors[iloc + 1][0]
            readgap_1 = anchors[iloc + 1][0] - anchors[iloc][1]
            refgap_1 = anchors[iloc + 1][2] - anchors[iloc][3]
            gap_1 = readgap_1 - refgap_1
            
            readgap_2 = anchors[iloc + 2][0] - anchors[iloc + 1][1]
            refgap_2 = anchors[iloc + 2][2] - anchors[iloc + 1][3]
            gap_2 = readgap_2 - refgap_2
            if(gap_1 > indelsize):
                if(gap_2 < 0):
                    if(minanchor > midsize):
                        anchors.pop(iloc + 1)
                        continue
                    else:
                        iloc += 1
                        continue
                else:
                    iloc += 1
                    continue
            elif(gap_1 < -indelsize):
                if(gap_2 > 0):
                    if(minanchor > midsize):
                        anchors.pop(iloc + 1)
                        continue
                    else:
                        iloc += 1
                        continue
                else:
                    iloc += 1
                    continue
            elif(gap_2 > indelsize):
                if(gap_1 < 0):
                    if(minanchor > midsize):
                        anchors.pop(iloc + 1)
                        continue
                    else:
                        iloc += 1
                        continue
                else:
                    iloc += 1
                    continue
            elif(gap_2 < -indelsize):
                if(gap_1 > 0):
                    if(minanchor > midsize):
                        anchors.pop(iloc + 1)
                        continue
                    else:
                        iloc += 1
                        continue
                else:
                    iloc += 1
                    continue
            else:
                iloc += 1
                continue
        pre_read_pos = anchors[0][0]
        pre_ref_pos = anchors[0][2]
        for iloc in range(len(anchors)):
            anchor = anchors[iloc]

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
        iloc = 0
        indelsize = 30
        while((iloc + 2) < len(anchors)):
            midsize = anchors[iloc + 1][1] - anchors[iloc + 1][0]
            readgap_1 = anchors[iloc + 1][0] - anchors[iloc][1]
            refgap_1 = anchors[iloc][2] - anchors[iloc + 1][3]
            gap_1 = readgap_1 - refgap_1
            
            readgap_2 = anchors[iloc + 2][0] - anchors[iloc + 1][1]
            refgap_2 = anchors[iloc + 1][2] - anchors[iloc + 2][3]
            gap_2 = readgap_2 - refgap_2
            if(gap_1 > indelsize):
                if(gap_2 < 0):
                    if(minanchor > midsize):
                        anchors.pop(iloc + 1)
                        continue
                    else:
                        iloc += 1
                        continue
                else:
                    iloc += 1
                    continue
            elif(gap_1 < -indelsize):
                if(gap_2 > 0):
                    if(minanchor > midsize):
                        anchors.pop(iloc + 1)
                        continue
                    else:
                        iloc += 1
                        continue
                else:
                    iloc += 1
                    continue
            elif(gap_2 > indelsize):
                if(gap_1 < 0):
                    if(minanchor > midsize):
                        anchors.pop(iloc + 1)
                        continue
                    else:
                        iloc += 1
                        continue
                else:
                    iloc += 1
                    continue
            elif(gap_2 < -indelsize):
                if(gap_1 > 0):
                    if(minanchor > midsize):
                        anchors.pop(iloc + 1)
                        continue
                    else:
                        iloc += 1
                        continue
                else:
                    iloc += 1
                    continue
            else:
                iloc += 1
                continue
        pre_read_pos = anchors[0][0]
        pre_ref_pos = anchors[0][3]
        for iloc in range(len(anchors)):
            anchor = anchors[iloc]
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
        
def get_cigar_anchors(index_object, aln, anchors, query, minanchor = 0):
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
        for iloc in range(len(anchors)):
            anchor = anchors[iloc]

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
        for iloc in range(len(anchors)):
            anchor = anchors[iloc]
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

def mergecigar_n(cigarstring):  
    oplist = mergecigar_(cigarstring)
    n_cigar = len(oplist)
    return ''.join(oplist), n_cigar

def P_alignmentstring_comments(infodict, comments):
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
    def checkcomment(comments, infolist, tags):
        if(isinstance(comments, str)):
            comments_info = comments.split('\t')
            for onecomment in comments_info:
                info = onecomment.split(':')
                if((len(info) == 3) and (len(info[0]) == 2) and (info[0] not in tags) and (info[1] in ('A', 'i', 'f', 'Z', 'H', 'B'))):
                    infolist.append(onecomment)
                    tags.add(info[0])
    tags = set(('QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL', 'SA','NM','MD','cs'))
    for key in infodict:
        if(key in name2iloc):
            infolist[name2iloc[key]] = infodict[key]
        else:
            infolist.append(p_other_tag(key, infodict[key]))
            tags.add(key)
    checkcomment(comments, infolist, tags)
    return '\t'.join(infolist)

def get_bam_dict_str(mapinfo, query, qual, comments, index_object, option):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    

    
    hardclip = option['H']
    
    rc_query = str(Seq(query).reverse_complement())
    

    iloc2nm = dict()
    iloc2md = dict()
    iloc2cs = dict()
    iloc2n_cigar = dict()
    tmpiloc = -1
    fakecigar = option['fakecigar']
    if(fakecigar == True):
        iloc2fakecigar = dict()
        if(hardclip == True):
            clipsyb = 'H'
        else:
            clipsyb = 'S'
    want_css_short, want_css_long, want_md = False, False, False
    if('shortcs' in option):
        if(option['shortcs'] == False):

            want_css_long = True
        else:
            want_css_short = True

    if('md' in option and option['md'] == True):
        want_md = True
    if('cigar2cg' in option):
        cigar2cg = True
    
    for item in mapinfo:
        item[-1], n_cigar = mergecigar_n(item[-1])
        tmpiloc += 1
        if(item[2] == '+'):
            if(option['H'] == True):
                tmp_query = query[item[3]: item[4]]
            else:
                tmp_query = query
            tmp_target = index_object.seq(item[1], item[5], item[6])
        else:
            if(option['H'] == True):
                tmp_query = rc_query[item[3]: item[4]]
            else:
                tmp_query = rc_query
            tmp_target = index_object.seq(item[1], item[5], item[6])
        try:
            mdstring, css, csl, nm = vacmap_index.k_md_cs(item[-1], tmp_target, tmp_query, want_md, want_css_short, want_css_long)
        except:
            print(item[-1])
            print()
            print(tmp_target)
            print()
            print(tmp_query)
            print()
            print(want_md, want_css_short, want_css_long)
            print()
            raise
        iloc2n_cigar[tmpiloc] = n_cigar
        iloc2nm[tmpiloc] = nm
        if(want_md == True):
            iloc2md[tmpiloc] = mdstring
        if(want_css_short == True):
            iloc2cs[tmpiloc] = css
        elif(want_css_long == True):
            iloc2cs[tmpiloc] = csl

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

        if(iloc2n_cigar[iloc] > 65535):
            if(cigar2cg == True):
                bam_dict['CG'] = primary[8]
                logging.info('Write long CIGAR to CG tag.')
            else:
                bam_dict['CIGAR'] = primary[8]
        else:
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
        if(mq != 0):
            mq = 60
        else:
            mq = 1
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
        if(want_md == True):
            bam_dict['MD'] = iloc2md[tmpiloc]
        if(want_css_short == True):
            bam_dict['cs'] = iloc2cs[tmpiloc]
        elif(want_css_long == True):
            bam_dict['cs'] = iloc2cs[tmpiloc]

        a_list.append(P_alignmentstring_comments(bam_dict, comments))
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

import bisect





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



def yield_alnstr_dirct(onemapinfo, query, qual, queryname, comments, index_object, option):
    mapinfo = []
    if(option['H'] == True):
        clipsyb = 'H'
    else:
        clipsyb = 'S'
    for line in onemapinfo:
        contig, strand, q_st, q_en, r_st, r_en, cigar, mapq = line
        if(q_st > 0):
            cigar = str(q_st) + clipsyb + cigar
        if((len(query) - q_en) > 0):
            cigar += str(len(query) - q_en) + clipsyb
        mapinfo.append([queryname, contig, strand, q_st, q_en, r_st, r_en, mapq, cigar])
        
    for onealn in get_bam_dict_str(mapinfo, query, qual, comments, index_object, option):
        yield onealn

def qualities_to_fastq_str(qualities):
    """
    Convert list of Phred qualities (ints) to FASTQ string.
    """
    return "".join(chr(q + 33) for q in qualities)


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
    
    

    if proc.stdin and not proc.stdin.closed:
        proc.stdin.close()


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
def lift_one_aln(stat_queue, cooked_queue, asmid2chains, asmid2max_chain_len,  index_object, missedpath, workerindex, bamfilepath, threads, option, sameref):
    string_cache = []
    string_cache_max = 1000

    fq = open(missedpath, 'w')
        

    batchsize = 100
    allbatchsize = threads * batchsize
    if(bamfilepath.endswith('sam')):
        bamfile = pysam.AlignmentFile(bamfilepath,'r', threads = 1)
    else:
        bamfile = pysam.AlignmentFile(bamfilepath,'rb', threads = 1)

    readcount_st, readcount_en = (workerindex - 1) * batchsize, workerindex * batchsize
    readcount = 0
    mapunmap = [0, 0]
    for AlignedSegment in bamfile:
        readcount += 1
        if(readcount_st <= readcount % allbatchsize < readcount_en):
            try:
                if((AlignedSegment.is_secondary == True) or (AlignedSegment.is_supplementary == True)):
                    continue
                query_name = AlignedSegment.query_name

                query = AlignedSegment.query_sequence
                asmid = AlignedSegment.reference_name
                strand = decode_strand_from_flag(AlignedSegment.flag)
                read_length = AlignedSegment.query_length
                aln_asm = [asmid, strand, 0, read_length, AlignedSegment.reference_start, AlignedSegment.reference_end, AlignedSegment.cigarstring]


                t_query = AlignedSegment.get_forward_sequence()

                quals = AlignedSegment.query_qualities
                if quals is None:
                    quals = [0] * AlignedSegment.query_length  # default values

                t_qual = quals if strand == '+' else quals[::-1]

                t_qual = qualities_to_fastq_str(t_qual)
                if(AlignedSegment.is_mapped == False):

                    fq.write(f"@{query_name}\n{t_query}\n+\n{t_qual}\n")
                    mapunmap[1] += 1
                    continue
                
                if((AlignedSegment.has_tag('SA') == True) or ((read_length - AlignedSegment.query_alignment_length) > option['largeclip'])):

                    fq.write(f"@{query_name}\n{t_query}\n+\n{t_qual}\n")
                    mapunmap[1] += 1
                    continue
                    
                if((sameref == False) and (asmid not in asmid2chains)):

                    fq.write(f"@{query_name}\n{t_query}\n+\n{t_qual}\n")
                    mapunmap[1] += 1
                    continue
                    
                if(check_largeindel(aln_asm[6], size = option['largeindel']) == True):

                    fq.write(f"@{query_name}\n{t_query}\n+\n{t_qual}\n")
                    mapunmap[1] += 1
                    continue
                if(sameref == True):
                    if(len(string_cache) > string_cache_max):
                        alnstring = '\n'.join(string_cache)
                        alnstring = alnstring.encode('utf-8') + b'\n'
                        cooked_queue.put(alnstring)
                        string_cache = []
                    string_cache.append(AlignedSegment.to_string())
                else:
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

                        extend_edge_anchors(alns, alns_anchors, mapqs, index_object, query, minanchor = option['minanchor'])
                        if(need_reverse == True):
                            onemapinfo = change_alns_coo_baseon_strand(change_alns_strand(alns, len(query)), len(query))
                        else:
                            onemapinfo = change_alns_coo_baseon_strand(alns, len(query))
      
                        comments = 'XR:Z:true'
                        for alnstr in yield_alnstr_dirct(onemapinfo, t_query, t_qual, query_name, comments, index_object, option):
                            if(len(string_cache) > string_cache_max):
                                alnstring = '\n'.join(string_cache)
                                alnstring = alnstring.encode('utf-8') + b'\n'
                                cooked_queue.put(alnstring)
                                string_cache = []
                            string_cache.append(alnstr)
                        mapunmap[0] += 1
                    else:
                        fq.write(f"@{query_name}\n{t_query}\n+\n{t_qual}\n")
                        mapunmap[1] += 1
                        continue
                    
            except:
                query_name = AlignedSegment.query_name
                t_query = AlignedSegment.get_forward_sequence()
                if(query_name):
                    if(t_query):
                        quals = AlignedSegment.query_qualities
                        if quals is None:
                            quals = [0] * AlignedSegment.query_length  # default values
                        t_qual = quals if strand == '+' else quals[::-1]
                        t_qual = qualities_to_fastq_str(t_qual)
                        fq.write(f"@{query_name}\n{t_query}\n+\n{t_qual}\n")
                mapunmap[1] += 1
                continue

    stat_queue.put(mapunmap)
    if(len(string_cache) > 0):
        alnstring = '\n'.join(string_cache)
        alnstring = alnstring.encode('utf-8') + b'\n'
        cooked_queue.put(alnstring)
    fq.close()
import numba
import multiprocessing
import time

def remap_func(asm2ref_path, read2asm_path, workdir, threads, datatype, refpath, index_object, header, option, sameref = False, realign_unlifted = True):
    def getminimap2option(option):
        minimap2option = []
        if(option['H'] == False):
            minimap2option.append('-Y')
        if(option['md'] == True):
            minimap2option.append('--MD')
        if('shortcs' in option):
            if(option['shortcs'] == True):
                minimap2option.append('--cs')
            else:
                minimap2option.append('--cs=long')
        return minimap2option
    def run_minimap_to_bam(minimap_args, output_path, sort = False):
        minimap_cmd = aligner_env + minimap_args
        if(sort == True):
            samtools_cmd = aligner_env + ['samtools', 'sort', '-@', '8', '-o', output_path, '-']
        else:
            samtools_cmd = aligner_env + ['samtools', 'view', '-@', '8', '-b', '-o', output_path, '-']
        with subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE) as minimap_proc:

            subprocess.run(samtools_cmd, stdin=minimap_proc.stdout, check=True)

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
        sorted_lifed_read_alns_path = workdir + 'read2ref.remapped.sorted.bam'
        missedreadprefix = workdir + 'filtered_read/'
        os.makedirs(missedreadprefix, exist_ok=True)
        if(os.path.isfile(sorted_lifed_read_alns_path) == False):
            if(sameref == False): 
                print('Loading ASM to REF Alignments to memory')
                convert_eqx = False
                bamfilepath = asm2ref_path
                if bamfilepath.endswith('sam'):
                    bamfile = pysam.AlignmentFile(bamfilepath, 'r', check_sq=False)
                else:
                    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', check_sq=False)

                contig2index = dict()

                iloc = -1
                for contigname, contiglength in zip(bamfile.references, bamfile.lengths):
                    iloc += 1
                    contig2index[contigname] = iloc



                aln_queue, cooked_queue = multiprocessing.Manager().Queue(maxsize=200), multiprocessing.Manager().Queue(maxsize=2000)
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

            else:

                asmid2chains = dict()
                asmid2max_chain_len = dict()
            stat_queue = multiprocessing.Manager().Queue()
            print('Remapping')



            bamfilepath = read2asm_path

            sorted_lifed_read_alns_path = workdir + 'read2ref.liftover.sorted.bam'

            

            missediloc = 0

            cooked_queue =  multiprocessing.Manager().Queue(maxsize=200)
            workers = []
            for iloc in range(threads):
                missediloc += 1
                worker = multiprocessing.Process(
                    target=lift_one_aln,
                    args=(
                        stat_queue,
                        cooked_queue,
                        asmid2chains,
                        asmid2max_chain_len,
                        index_object,
                        missedreadprefix + str(missediloc) + '.fastq',
                        missediloc,
                        bamfilepath,
                        threads, 
                        option,
                        sameref,
                    ),
                )
                worker.start()
                workers.append(worker)
            write_process = multiprocessing.Process(
                target=mp_write_bam_sorted,
                args=(cooked_queue, header, sorted_lifed_read_alns_path)
            )
            write_process.start()
            mapunmap = [0, 0]
            for work_p in workers:
                tmpmapunmap = stat_queue.get()
                mapunmap[0] += tmpmapunmap[0]
                mapunmap[1] += tmpmapunmap[1]
                work_p.join()
            cooked_queue.put(0)
            write_process.join()

            print('Remapped reads:', str(mapunmap[0])+'/'+str(mapunmap[0] + mapunmap[1]))
            print('Remapping completed')
        

        #print(f"Remapping Wall time (real elapsed): {format_seconds(wall_time_sec)} ({wall_time_sec:.4f} seconds)")
        #print(f"Remapping CPU time (processor used): {format_seconds(cpu_time_sec)} ({cpu_time_sec:.4f} seconds)")
        
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
            minimap2option = getminimap2option(option)
            minimap_args = [
                'minimap2',
                '-I', '1000G',
                '-ax', minimap_presets[datatype]] + minimap2option + [
                '-t', str(threads),
                '--eqx',
                refpath,
            ] + missed_reads_fastqs
            run_minimap_to_bam(minimap_args, missedread2refpath_bam, sort = True)
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
        
        minimap2option = getminimap2option(option)
        minimap_args = [
            'minimap2',
            '-I', '1000G',
            '-ax', minimap_presets[datatype]] + minimap2option + [
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


def get_liftover_data(inputreadpath, refpath, workdir, datatype, threads, asmdata='', option = dict(), self_align=True, rutg = True):
    assert datatype in ('ont', 'hqlr', 'hifi')
    aligner_env_name = 'minimap2_env'
    aligner_env = ["conda", "run", "-n", aligner_env_name]
    aligner_env = []
    asmfasta = asmdata
    def run_minimap_to_bam(minimap_args, output_path):
        minimap_cmd = aligner_env + minimap_args
        samtools_cmd = aligner_env + ['samtools', 'view', '-@', '8', '-b', '-o', output_path, '-']
        with subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE) as minimap_proc:

            subprocess.run(samtools_cmd, stdin=minimap_proc.stdout, check=True)

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

    if (asmfasta == ''):
        if(not os.path.isfile(asmout + '.bp.r_utg.gfa')):
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
            if(not os.path.isfile(asmfasta)):
                count = gfa_to_fasta([asmout + '.bp.hap1.p_ctg.gfa', asmout + '.bp.hap2.p_ctg.gfa'], asmfasta)
            else:
                count = 1
        if count == 0:
            asmfasta = ''
    elif asmfasta[-3:] == 'gfa':
        print('Assembly provided')
        count = gfa_to_fasta(asmfasta, asmout + '.bp.r_utg.fasta')
        asmfasta = asmout + '.bp.provided.fasta'
        if count == 0:
            asmfasta = ''
    if(asmfasta != ''):
        if(('asm2ref'not in option) or (option['asm2ref'] == '')):
            asm2ref_path = workdir + 'asm2ref.bam'
        else:
            asm2ref_path = option['asm2ref']
            assert os.path.isfile(asm2ref_path) == True
        if(('read2asm'not in option) or (option['read2asm'] == '')):
            read2asm_path = workdir + 'read2asm.bam'
        else:
            read2asm_path = option['read2asm']
            assert os.path.isfile(read2asm_path) == True
        if(asmfasta != refpath):
            if(os.path.isfile(asm2ref_path) == False):
                print('Aligning assembly to reference')
                a = run_minimap_to_bam(
                    ['minimap2', '-I', '1000G', '-ax', 'asm5', '--cs', '-r2k', '-t', str(threads), '--eqx', refpath, asmfasta],
                    asm2ref_path
                )
                assert a.returncode == 0
        if(os.path.isfile(read2asm_path) == False):
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
                    a = run_minimap_to_bam(
                        ['minimap2', '-I', '1000G', '-aYx', 'lr:hq', '-t', str(threads), '--eqx', asmfasta] + input_reads,
                        read2asm_path
                    )
            assert a.returncode == 0
        return asm2ref_path, read2asm_path
    else:
        return '', inputreadpath




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

class FastaReader:
    

    def __init__(self, filepath):
        
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"FASTA file not found at: {filepath}")

        self.sequences = {}

        self.name_len = []

        for rec in vacmap_index.fastx_read(filepath):

            self.sequences[rec[0]] = rec[1].upper()#.encode('latin-1')
            self.name_len.append((rec[0], len(rec[1])))

    def seq(self, contig_name, start, end):
        return self.sequences[contig_name][start:end]
    
RG_ARG_MAP: dict[str, str] = {
    "rg-id": "ID",
    "rg-sm": "SM",
    "rg-lb": "LB",
    "rg-pl": "PL",
    "rg-ds": "DS",
    "rg-dt": "DT",
    "rg-pu": "PU",
    "rg-pi": "PI",
    "rg-pg": "PG",
    "rg-cn": "CN",
    "rg-fo": "FO",
    "rg-ks": "KS",
    "rg-pm": "PM",
    "rg-bc": "BC",
}


def collect_rg_metadata(args: argparse.Namespace, option) -> dict[str, str]:
    """Build a SAM @RG dictionary from CLI args."""
    rg_values: dict[str, str] = {}
    for attr, tag in RG_ARG_MAP.items():
        value = getattr(args, attr, None)
        if value is None:
            continue
        rg_values[tag] = str(value)
        option[attr] = str(value)

    if rg_values and "ID" not in rg_values:
        raise ValueError("The --rg-id option is required when any other --rg-* option is supplied.")
    return rg_values
def format_seconds(seconds):
    """Converts a floating point number of seconds into hh:mm:ss.ms format."""
    # Split the seconds into minutes/seconds/milliseconds
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    # Format the components as strings
    # Use '{:02.4f}'.format(seconds) to keep 4 decimal places for seconds component
    return f"{int(hours):02}:{int(minutes):02}:{seconds:07.4f}"
def run_pipeline(args: argparse.Namespace) -> None:
    asmdata = args.asmdata or ""
    if(args.workdir[-1] != '/'):
        args.workdir += '/'
    # Determine hybrid mode from flags.
    sameref = False
    if args.localasm == True:
        if(asmdata == ''):
            asmdata = args.refpath
            sameref = True
        self_align = False
        realign_unlifted = False
    else:
        self_align = True
        realign_unlifted = True

    if asmdata and not Path(asmdata).exists():
        raise FileNotFoundError(f"{asmdata} not exist.")

    if args.datatype == "ont" and not asmdata:
        raise RuntimeError("Old ONT data requires --asmdata.")
    if args.datatype == "ont" and args.localasm:
        raise RuntimeError("Old ONT data not supported with local assembly")
    if(args.datatype == 'hifi' and args.minanchor == 0):
        args.minanchor = 5
        
    rutg = args.rutg
    option = {'H': args.H, 'fakecigar': True, 'largeindel': args.largeindel,
             'largeclip': args.largeclip, 'asm2ref': args.asm2ref,
             'read2asm': args.read2asm, 'md':args.md, 'minanchor': args.minanchor}
    
    if(args.csshort == True):
        option['shortcs'] = True
    elif(args.cslong == True):
        option['shortcs'] = False
    

    asm2ref_path, read2asm_path = get_liftover_data(
        args.inputreadpath,
        args.refpath,
        args.workdir,
        args.datatype,
        args.threads,
        asmdata=asmdata,
        option=option,
        rutg=rutg,
        self_align=self_align,
    )

    index_object = FastaReader(args.refpath)
    if not index_object:
        raise RuntimeError("ERROR: failed to load index")


    header = {"HD": {"VN": "1.0"}, "SQ": [], "PG": [{
        "ID": "remap",
        "PN": "remap",
        "VN": "1.0",
        "CL": " ".join(sys.argv)
    }]}
    rg_metadata = collect_rg_metadata(args, option)

    if rg_metadata:
        header["RG"] = [rg_metadata]
    for seq_name, seq_len in index_object.name_len:

        header["SQ"].append(
            {"LN": seq_len, "SN": seq_name}
        )
    

    
    firstbam, missedreadpath = remap_func(
        asm2ref_path,
        read2asm_path,
        args.workdir,
        args.threads,
        args.datatype,
        args.refpath,
        index_object,
        header,
        option,
        sameref=sameref,
        realign_unlifted=realign_unlifted
    )
    

    
    secondbam = ""
    if(args.localasm == True):
        remap_workdir = Path(missedreadpath).resolve().parent / "REMAP_data"
        remap_workdir.mkdir(parents=True, exist_ok=True)

        asm2ref_path, read2asm_path = get_liftover_data(
            missedreadpath,
            args.refpath,
            str(remap_workdir)+'/',
            args.datatype,
            args.threads,
            asmdata="",
            option=dict(),
            rutg=rutg,
            self_align=True,
        )
        secondbam, missedreadpath = remap_func(
            asm2ref_path,
            read2asm_path,
            str(remap_workdir)+'/',
            args.threads,
            args.datatype,
            args.refpath,
            index_object,
            header,
            option,
            sameref=False,
            realign_unlifted=True,
        )

    if not firstbam:
        print("No Realigned BAM generated.", file=sys.stderr)
        return

    if not missedreadpath.endswith("bam"):
        raise RuntimeError("Error: missed read path is not a BAM file.")

    samtools_threads = min(args.threads, 8)
    lifted_bam = Path(f"{args.outputprefix}remapped.bam").resolve()
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
    parser.add_argument('-i', "--inputreadpath", dest="inputreadpath", required=True, help="Glob for read files.")
    parser.add_argument('-r', "--refpath", dest="refpath", required=True, help="Reference FASTA.")
    parser.add_argument('-w', "--workdir", dest="workdir", required=True, help="Working directory.")
    parser.add_argument(
        '-o',
        "--outputprefix",
        dest="outputprefix",
        required=True,
        help="Prefix for output BAMs, e.g. /path/to/sample_.",
    )
    parser.add_argument(
        '-d',
        "--datatype",
        dest="datatype",
        required=True,
        choices=("ont", "hqlr", "hifi"),
        help="Read type R9.x: ont, R10.x: hqlr, PacBio HiFi: hifi.",
    )
    parser.add_argument(
        '-t',
        "--threads",
        dest="threads",
        default=8,
        type=int,
        help="Number of threads.",
    )
    parser.add_argument(
        '-a',
        "--asmdata",
        default="",
        help="Optional assembly file (FASTA or GFA).",
    )
    parser.add_argument(
        "--asm2ref",
        dest="asm2ref",
        default="",
        help="Optional assembly to reference alignment (BAM or SAM).",
    )
    parser.add_argument(
        "--read2asm",
        dest="read2asm",
        default="",
        help="Optional read to assembly alignment (BAM or SAM).",
    )

    
    parser.add_argument(
        "--rutg",
        dest='rutg',
        default=False,
        action="store_true",
        help="Use raw unitig graph.",
    )
    parser.add_argument(
        "--minanchor",
        default=0,
        type=int,
        help="Minimum size (base pair) of an anchor.",
    )
    parser.add_argument(
        "--largeindel",
        default=30,
        type=int,
        help="Direct mapping read contains large indels in CIGAR string",
    )
    parser.add_argument(
        "--largeclip",
        default=100,
        type=int,
        help="Direct mapping read with large clip",
    )

    parser.add_argument(
        "--localasm",
        dest="localasm",
        default=False,
        action="store_true",
        help="Only assemble reads contain SVs.",
    )
    parser.add_argument(
        "--cs",
        dest="csshort",
        default=False,
        action="store_true",
        help="Add short cs tag.",
    )
    parser.add_argument(
        "--cs=long",
        dest="cslong",
        default=False,
        action="store_true",
        help="Add long cs tag.",
    )
    parser.add_argument(
        "--MD",
        dest="md",
        default=False,
        action="store_true",
        help="Add MD tag.",
    )

    parser.add_argument(
        '--H',
        "--hardclip",
        dest="H",
        default=False,
        action="store_true",
        help="Use hardclip instead of softclip.",
    )
    rg_group = parser.add_argument_group("Read-group (RG) metadata")
    rg_group.add_argument("--rg-id", dest="rg-id", metavar="STRING",
                          help="Adds RG:Z:<string> to all alignments (required when using other --rg-* flags).")
    rg_group.add_argument("--rg-sm", dest="rg-sm", metavar="STRING", help="RG header: Sample (SM).")
    rg_group.add_argument("--rg-lb", dest="rg-lb", metavar="STRING", help="RG header: Library (LB).")
    rg_group.add_argument("--rg-pl", dest="rg-pl", metavar="STRING", help="RG header: Platform (PL).")
    rg_group.add_argument("--rg-ds", dest="rg-ds", metavar="STRING", help="RG header: Description (DS).")
    rg_group.add_argument("--rg-dt", dest="rg-dt", metavar="YYYY-MM-DD", help="RG header: Date (DT).")
    rg_group.add_argument("--rg-pu", dest="rg-pu", metavar="STRING", help="RG header: Platform unit (PU).")
    rg_group.add_argument("--rg-pi", dest="rg-pi", metavar="INT", help="RG header: Median insert size (PI).")
    rg_group.add_argument("--rg-pg", dest="rg-pg", metavar="STRING", help="RG header: Programs (PG).")
    rg_group.add_argument("--rg-cn", dest="rg-cn", metavar="STRING", help="RG header: Sequencing center (CN).")
    rg_group.add_argument("--rg-fo", dest="rg-fo", metavar="STRING", help="RG header: Flow order (FO).")
    rg_group.add_argument("--rg-ks", dest="rg-ks", metavar="STRING", help="RG header: Key sequence (KS).")
    rg_group.add_argument("--rg-pm", dest="rg-pm", metavar="STRING",
                          help="RG header: Platform model (PM).")
    rg_group.add_argument("--rg-bc", dest="rg-bc", metavar="STRING",
                          help="RG header: Barcode sequence (BC).")

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:

    args = parse_args(argv)
    run_pipeline(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
