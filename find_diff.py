from itertools import islice
import pickle
import sys
import os

### Pseudocode
### Check if piscem_id is not keys of pisc_comm_but_multihits given a multihit
### If yes - increment the piscem iterator
### If not a match, 
###     If piscem position smaller, increment the iterator for piscem, else for chromap
### Else
###     
### Match is defined is 0 <= start_pos_chromap - start_pos_piscem <= 6 and barcodes map
### If Match:
###     If multihits:
###         check if the key exists in pisc_diff_hits
###             if yes - remove it from there and put
###         add it to the pisc_comm_but_multihits
###    
###     increment counter for common
###     move iterator for both chromap and piscem
### If not a match:
###     if pisc < chrom_pos
###         put the piscem value in pisc_diff_hits
###         pisc++
###     else
###         put the chromap value in chrom_unique
###         chrom++


def is_match(chrom_l, pisc_l):
    if chrom_l[3] == pisc_l[3]:
        diff = int(chrom_l[1]) - int(pisc_l[1])
        if diff > 2 and diff < 6:
            return True
    return False

def add_from_file_from_spot(start_ind, fmap, nlines, map_diff, map_comm=None,dedup=None):
    with open(fmap, "r") as f1:
        chunk = list(islice(f1, start_ind))
        while True:
            chunk = list(islice(f1, nlines))
            if not chunk:
                break
            add_miss(0, chunk, nlines, map_diff, map_comm, dedup)


def add_miss(start_ind, chunk, nlines, map_diff, map_comm=None, dedup=None):
    ind = start_ind
    while ind < len(chunk):
        # print(ind)
        l_chunk = chunk[ind%nlines]
        l_chunk = l_chunk.split("\n")[0].split("\t")
        pos = int(l_chunk[1])
        if map_comm is not None:
            rec_id = int(l_chunk[-1])
            n_hits = int(l_chunk[-2])
            
            if (rec_id not in map_comm.keys()) and (rec_id not in dedup.keys()):
                map_diff[rec_id] = (pos, n_hits)
        else:
            # print("natch", pos)
            map_diff.append((pos, l_chunk[-2]))
        ind += 1
    return ind

def proc_chunks(start_ind, flag_chunks, chunks, pisc_diff_hits, pisc_comm_hits, chrom_unique, dedup, nlines, counter):
    len_chr = len(chunks[0])
    len_pisc = len(chunks[1])
    inds = [start_ind[0]%nlines, start_ind[1]%nlines]
    prev_pos = -1
    prev_barcode = ""
    prev_rec=-1
    while (inds[0] < len_chr and inds[1] < len_pisc):
        l_chr = chunks[0][inds[0]] ## chromap chunk
        l_chr = l_chr.split("\n")[0].split("\t")
        l_pisc = chunks[1][inds[1]] ## piscem chunk
        l_pisc = l_pisc.split("\n")[0].split("\t")
        rec_id = int(l_pisc[-1])
        n_hits = int(l_pisc[-2])
        pos_pisc = int(l_pisc[1])
        pos_chr = int(l_chr[1])
        bar_pis = l_pisc[3]

        # print('rec', rec_id)
        if n_hits >= 2:
            if rec_id in pisc_comm_hits.keys() or rec_id in dedup.keys():
                inds[1] += 1
                start_ind[1] += 1
                continue
        
        match = is_match(l_chr, l_pisc)
        if not match:
            if (pos_pisc - prev_pos) <= 2 and bar_pis == prev_barcode: ## this would have been a valid hit but already accounted for
                if rec_id != prev_rec:
                    dedup[rec_id] = n_hits
                    start_ind[1] += 1
                    inds[1] += 1
                    continue
        if match:
            # print("match", pos_chr)
            if n_hits >= 2:
                pisc_diff_hits.pop(rec_id, None)
                ## Only a single match for piscem - since chromap outputs only a single hit
                ## Thus once the match is found, we wont find any other match
            pisc_comm_hits[rec_id] = (pos_pisc, n_hits)
            prev_pos = pos_pisc
            prev_barcode = bar_pis
            prev_rec = rec_id
            counter += 1
            inds[0] += 1
            inds[1] += 1
            start_ind[0] += 1
            start_ind[1] += 1
        else:
            if pos_pisc < pos_chr:
                pisc_diff_hits[rec_id] = (pos_pisc, n_hits)
                inds[1] += 1
                start_ind[1] += 1
            else:
                # print("unmatch", pos_chr)
                chrom_unique.append((pos_chr, l_chr[-2]))
                inds[0] += 1
                start_ind[0] += 1
    # print(inds)
    if inds[0] == len_chr:
        # start_ind[0:2] = [0, inds[1]]
        flag_chunks[0:2] = [True,False]
    else:
        # start_ind[0:2] = [inds[0], 0]
        flag_chunks[0:2] = [False, True]
    # print(flag_chunks)
    
    return counter

def find_diff(fchromap, fpiscem, pisc_diff_hits, pisc_comm_hits, chrom_unique, dedup, nlines=1000, counter=0):
    total_lines=0
    start_ind = [0,0] ### Starting indexes in chunks
    
    flag_chunks = [True, True] ### whether to read next chunk or not from the files
    with open(fchromap, "r") as f1, open(fpiscem, "r") as f2:
        while True:
            if flag_chunks[0]:
                chunk_chr = list(islice(f1, nlines))
            if flag_chunks[1]:
                chunk_pisc = list(islice(f2, nlines))
            if not chunk_chr:
                break
            if not chunk_pisc:
                break
            counter = proc_chunks(start_ind, flag_chunks, [chunk_chr, chunk_pisc], 
                pisc_diff_hits, pisc_comm_hits, chrom_unique, dedup, nlines, counter)
            print(start_ind)
    # print(len(chrom_unique))
    # print(chrom_unique)
    if chunk_pisc:
        # print("yo")
        start_ind[1] = add_miss(start_ind[1], chunk_pisc, nlines, pisc_diff_hits, pisc_comm_hits, dedup)
        add_from_file_from_spot(start_ind[1], fpiscem, nlines, pisc_diff_hits, pisc_comm_hits, dedup)
    if chunk_chr:
        # print("sup")
        start_ind[0] = add_miss(start_ind[0], chunk_chr, nlines, chrom_unique)
        # print(start_ind)
        add_from_file_from_spot(start_ind[0], fchromap, nlines, chrom_unique)
    # print(len(chunk_chr), len(chunk_pisc))
    return counter
    ### Still to implement remaining chunk part, writing to file

def main():
    # fchromap="sorted_chromap.bed"
    # fpiscem="sorted_piscem_all.bed"
    fchromap = sys.argv[1]
    fpiscem = sys.argv[2]
    pick_path = sys.argv[3]
    pick_prefix = sys.argv[4]
    chunk1=[]
    chunk2=[]
    ## If either barcode does not exist in chromap or the position does not exist in piscem
    pisc_diff_hits={} ## Dictionary with key rec id and value tuple (loc, nhits)
    
    ## If the positions match between 2 but multiple hits
    pisc_comm_but_multihits={} ## Dictionary with key rec id and value tuple (loc, nhits)
    pisc_comm_hits={}
    dedup = {} ### Dictionary with deduplication
    
    ## If either barcode does not exist in piscem or the position does not map in chromap
    chrom_unique=[] ## List of tuples (loc, barcode)
    nlines = 1000000
    counter=0
    counter=find_diff(fchromap, fpiscem, pisc_diff_hits, pisc_comm_hits, chrom_unique, dedup, nlines=nlines, counter=counter)
    print("counter", counter)
    print("dedup")
    print(len(dedup))
    print("pisc diff")
    print(len(pisc_diff_hits))
    print("chrom unique")
    print(len(chrom_unique))
    sum1=0
    sum4=0
    sum2=0
    sum3=0
    for k in pisc_comm_hits.keys():
        vals = pisc_comm_hits[k]
        sum1+=vals[1]
        if vals[1] >= 2:
            pisc_comm_but_multihits[k]=vals
            sum4+=vals[1]
    for v2 in pisc_diff_hits.values():
        sum2+=v2[1]
    for v2 in dedup.values():
        sum3+=v2
    print("piscem all hits, miss hits, dedup, rem")
    print(sum1, sum2, sum3, sum4)
    with open(os.path.join(pick_path, pick_prefix+"_chrom_uniq.pickle"), "wb") as f:
        pickle.dump(chrom_unique, f)
    with open(os.path.join(pick_path, pick_prefix + "_pisc_uniq.pickle"), "wb") as f:
        pickle.dump(pisc_diff_hits, f)
    with open(os.path.join(pick_path, pick_prefix + "_pisc_comm_hits.pickle"), "wb") as f:
        pickle.dump(pisc_comm_hits, f)
    with open(os.path.join(pick_path, pick_prefix + "_pisc_mult_hits.pickle"), "wb") as f:
        pickle.dump(pisc_comm_but_multihits, f)
    with open(os.path.join(pick_path, pick_prefix + "_pisc_dedup.pickle"), "wb") as f:
        pickle.dump(dedup, f)

if __name__=="__main__":
    main()
