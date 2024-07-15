import pickle
import sys
import os
from itertools import islice

## assumes sorted score for sam file

def proc_sam_reads_dict(l1, l2, mappings, vals, inds):
    # print(inds[0])
    l1=l1.split("\n")[0].split("\t")
    l2=l2.split("\n")[0].split("\t")

    prev_read = vals[0]
    score_flag = vals[1]
    i = vals[2] ## i denotes the current pointer in the mappings
    read_name=l1[0]
    inds[0] += 1

    flag1 = int(l1[1])
    flag2 = int(l2[1])
    
    is_paired = (flag1 & 1 == 1) and (flag2 & 1 == 1)
    is_proper_pair = (flag1 & 2 == 2) and (flag2 & 2 == 2)
    if read_name!=l2[0]:
        if not is_proper_pair:
            inds[0] -= 1
            return [read_name, 0, i] 
        sys.exit("incompatible read names")
    if score_flag == -1 and prev_read == read_name:
        return [read_name, -1, i]
    if is_paired and is_proper_pair:
        score1 = int(l1[11].split(":")[2].strip())
        score2 = int(l2[11].split(":")[2].strip())
        score = score1+score2
        if read_name != prev_read :
            mappings[read_name]=[score,1] ## new read, thus create an entry, score_flag 0 means check, -1 means prev_read same as cur_read but score is inferior
            return [read_name, 0, i+1]
        else: ## same read
            if score > mappings[read_name][0]:
                mappings[read_name]=[score,1]
                return [read_name, 0, i]
            elif mappings[read_name][0] == score: ## append s
                mappings[read_name][1] += 1
                return [read_name, 0, i]
            else:
                return [read_name, -1, i]
    return [read_name, 0, i]

def proc_sam_reads_list(l1, l2, mappings, vals, inds):
    # print(inds[0])
    l1=l1.split("\n")[0].split("\t")
    l2=l2.split("\n")[0].split("\t")

    prev_read = vals[0]
    best_score = vals[1]
    score_flag = vals[2]
    i = vals[3] ## i denotes the current pointer in the mappings
    read_name=l1[0]
    inds[0] += 1

    flag1 = int(l1[1])
    flag2 = int(l2[1])
    
    is_paired = (flag1 & 1 == 1) and (flag2 & 1 == 1)
    is_proper_pair = (flag1 & 2 == 2) and (flag2 & 2 == 2)
    if read_name!=l2[0]:
        if not is_proper_pair:
            inds[0] -= 1
            return [read_name,best_score, 0, i] 
        sys.exit("incompatible read names")
    if score_flag == -1 and prev_read == read_name:
        return [read_name,-1, -1, i]
    if is_paired and is_proper_pair:
        score1 = int(l1[11].split(":")[2].strip())
        score2 = int(l2[11].split(":")[2].strip())
        score = score1+score2
        if read_name != prev_read :
            mappings.append(1) ## new read, thus create an entry, score_flag 0 means check, -1 means prev_read same as cur_read but score is inferior
            return [read_name, score, 0, i+1]
        else: ## same read
            if score > best_score:
                mappings[i]=1
                return [read_name, score, 0, i]
            elif best_score == score: ## append s
                mappings[i]+=1
                return [read_name, score, 0, i]
            else:
                return [read_name, best_score, -1, i]
    return [read_name, best_score, 0, i]
    
def best_mapping_sam(sam_file, pickle_file, keep_recname=False, num_chunks=100): 
    check_header = True
    first = True
    read_name=""
    inds=[0]
    mappings = []
    vals = [read_name, -1, 0, -1]
    
    if keep_recname:
        vals = [read_name, 0, -1]    
        mappings = {}
    
    with open(sam_file, 'r') as f:
        while True:
            inds[0]=0
            chunk = list(islice(f, num_chunks))
            l=len(chunk)
            if not chunk:
                break
            if check_header:
                while inds[0] < l:
                    if chunk[inds[0]].startswith("@"):
                        if inds[0]==0:
                            if chunk[inds[0]+l-1].startswith("@"):    
                                inds[0] += l-1
                            else:
                                check_header=False
                    else:    
                        if first:
                            check_header=False
                            first=False
                            if inds[0]%2 == 1:
                                cc = list(islice(f, 1))
                                if cc:
                                    chunk.extend(cc)
                                    l += 1
                        ind=inds[0]
                        if keep_recname:
                            vals = proc_sam_reads_dict(chunk[inds[0]], chunk[inds[0]+1], mappings, vals, inds)
                        else:
                            vals = proc_sam_reads_list(chunk[inds[0]], chunk[inds[0]+1], mappings, vals, inds)
                        if inds[0]==ind:
                            cc = list(islice(f, 1))
                            if cc:
                                chunk.extend(cc)
                                l += 1
                        
                    inds[0] += 1
            else:
                while inds[0] < l:
                    ind=inds[0]
                    if keep_recname:
                        vals = proc_sam_reads_dict(chunk[inds[0]], chunk[inds[0]+1], mappings, vals, inds)
                    else:
                        vals = proc_sam_reads_list(chunk[inds[0]], chunk[inds[0]+1], mappings, vals, inds)
                    if inds[0]==ind:
                        cc = list(islice(f, 1))
                        if cc:
                            chunk.extend(cc)
                            l += 1
                    inds[0] += 1
    with open(pickle_file, "wb") as ff:
        pickle.dump(mappings, ff)
    return mappings

def count_multimapping(bed_file, pickle_file="", num_chunks=1000000):
    skip = 0
    map_count = dict()
    ind = 0
    with open(bed_file, "r") as f:
        while True:
            chunk = list(islice(f, num_chunks))
            if not chunk:
                break
            while ind < len(chunk):
                l_chunk = chunk[ind%num_chunks]
                l_chunk = l_chunk.split("\n")[0].split("\t")
                count = int(l_chunk[4])
                if count not in map_count.keys():
                    map_count[count] = 0
                map_count[count] += 1
                ind = ind+count
            ind = ind%num_chunks
    if pickle_file != "":
        with open(pickle_file, "wb") as f:
            pickle.dump(map_count, f)
    return map_count
        

## Assuming barcode is unique
def extract_mappings_from_bed(bed_file):
    bed_mappings = {}
    with open(bed_file, "r") as f:
        for line in f:
            fields = line.split('\t')
            ref = fields[0]
            pos = int(fields[1])
            barcode = fields[3]
            if barcode not in bed_mappings.keys():
                bed_mappings[barcode] = []
            bed_mappings[barcode].append((ref, pos))
    return bed_mappings

def comp_bed_bed(bed_file1, bed_file2, barcode_map_pickle, mismatch_thresh = 10):
    barcode_map = pickle.load(open(barcode_map_pickle, 'rb'))[0]
    bed_map1 = extract_mappings_from_bed(bed_file1)
    bed_map2 = extract_mappings_from_bed(bed_file2)
    miss_barcodes_bed1 = []
    miss_barcodes_bed2 = []
    missing_both = []
    mismatch_query = {}
    matches = []
    count_matches = 0

    for b in barcode_map.keys():    
        barcode = barcode_map[b]
        key_in_bed1 = barcode in bed_map1.keys()
        key_in_bed2 = barcode in bed_map2.keys()
        if not key_in_bed1 and not key_in_bed2:
            missing_both.append(barcode)
        elif not key_in_bed1:
            miss_barcodes_bed1.append(barcode)
        elif not key_in_bed2:
            miss_barcodes_bed2.append(barcode)
        else:
            bed_vals1 = bed_map1[barcode] ## A list of mapping
            bed_vals1 = sorted(bed_vals1, key=lambda x:int(x[1]))
            
            bed_vals2 = bed_map2[barcode] ## A list of mapping
            bed_vals2 = sorted(bed_vals2, key=lambda x:int(x[1]))
            
            pos_bed1=-1
            pos_bed2=-1
            miss = True
            i=0
            j=0
            while i < len(bed_vals1) and j < len(bed_vals2):
                bed_val1 = bed_vals1[i]
                bed_val2 = bed_vals2[j]
                pos_bed1 = int(bed_val1[1])
                pos_bed2 = int(bed_val2[1])
                if abs(pos_bed1 - pos_bed2) <= mismatch_thresh:
                    matches.append(barcode)
                    miss = False
                    count_matches += 1
                    break
                elif pos_bed1 < pos_bed2:
                    i += 1
                else:
                    j += 1
            if miss:
                mismatch_query[barcode] = [pos_bed1, pos_bed2]

    print("count matches are", count_matches)
    return (mismatch_query, miss_barcodes_bed1, miss_barcodes_bed2, missing_both, matches)

def comp_samp_map(sam_file, bed_file, barcode_map_pickle, mismatch_thresh = 10):
    barcode_map = pickle.load(open(barcode_map_pickle, 'rb'))[0]
    sam_map = extract_mappings_from_sam(sam_file)
    bed_map = extract_mappings_from_bed(bed_file)
    miss_barcodes_bed = []
    miss_barcodes_sam = []
    missing_both = []
    mismatch_query = {}
    matches = []
    count_matches = 0
    
    for b in barcode_map.keys():    
        barcode = barcode_map[b]
        key_in_sam = b in sam_map.keys()
        key_in_bed = barcode in bed_map.keys()
        if not key_in_bed and not key_in_sam:
            missing_both.append(barcode)
        elif not key_in_bed:
            miss_barcodes_bed.append(barcode)
        elif not key_in_sam:
            miss_barcodes_sam.append(barcode)
        else:
            sam_val_read1 = sam_map[b]['read1']
            sam_val_read2 = sam_map[b]['read2']
            bed_vals = bed_map[barcode] ## A list of mapping
            bed_vals = sorted(bed_vals, key=lambda x:int(x[1]))
            pos_sam=-1
            pos_bed=-1
            miss = True
            for bed_val in bed_vals:
                if bed_val[0] == sam_val_read1[0] or bed_val[0] == sam_val_read2[0]:
                    pos_sam1 = int(sam_val_read1[1])
                    pos_sam2 = int(sam_val_read2[1])
                    pos_bed = int(bed_val[1])
                    pos_sam = pos_sam1 if abs(pos_sam1-pos_bed) < abs(pos_sam2-pos_bed) else pos_sam2
                    if abs(pos_sam - pos_bed) <= mismatch_thresh:
                        count_matches += 1
                        matches.append(barcode)
                        miss = False
                        break
                    if (pos_bed-pos_sam) > 100:
                        break
            if miss:
                mismatch_query[barcode] = [pos_sam, pos_bed]
                    
    print("count matches are", count_matches)
    return (mismatch_query, miss_barcodes_bed, miss_barcodes_sam, missing_both, matches)

def extract_mappings_from_sam(sam_file):
    mappings = {}
    
    with open(sam_file, 'r') as f:
        for line in f:
            if not line.startswith('@'):  # Skip header lines
                fields = line.split('\t')
                read_name = fields[0]
                flag= int(fields[1])
                ref = fields[2]
                pos = int(fields[3])
                is_paired = flag & 1 == 1
                is_proper_pair = flag & 2 == 2

                if is_paired and is_proper_pair:
                    if read_name not in mappings:
                        mappings[read_name] = {'read1': None, 'read2': None}
                    
                    read_type = "read1"
                    if flag & 128 == 128:  # Read 2
                        read_type = "read2"
                    
                    mappings[read_name][read_type] = (ref, pos)

    return mappings

def main():
    file1 = sys.argv[1] ## could be sam or bed
    file2 = sys.argv[2]
    bar_pickle = sys.argv[3]
    file_type = sys.argv[4]
    if file_type == "sam":
        out = comp_samp_map(file1, file2, bar_pickle)
        print(len(out[0]), len(out[1]), len(out[2]), len(out[3]))
    # pick_path = sys.argv[2]
    # pick_prefix = sys.argv[3]

    # mappings = extract_mappings_from_sam(sam_file)
    # with open(os.path.join(pick_path, pick_prefix) +".pi", "wb") as f:
    #     pickle.dump(mappings, f)

if __name__ == "__main__":
    main()