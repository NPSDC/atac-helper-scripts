import sys
import os
from itertools import islice
import pickle

def parse_files(barcode_file, read1_file, read2_file, size=1000):
    if not size %4 == 0:
        sys.exit("Size should be a multiple of 4")
    barcodes = {}
    read1_seq = []
    read2_seq = []
    with open(barcode_file, "r") as f1, open(read1_file, "r") as f2, open(read2_file, "r") as f3:
        while True:
            barcode_chunk = list(islice(f1, size))
            read1_chunk = list(islice(f2, size))
            read2_chunk = list(islice(f3, size))

            if not barcode_chunk:
                break
            if len(barcode_chunk) != len(read1_chunk) or len(read2_chunk) != len(read1_chunk):
                sys.exit("Chunks read incorrectly")
            ind = 0
            print(ind)
            while ind < len(barcode_chunk):
                header = barcode_chunk[ind].split(" ")[0][1:]
                barcodes[header] = barcode_chunk[ind+1][:-1]
                read1_seq.append(read1_chunk[ind+1][:-1])
                read2_seq.append(read2_chunk[ind+1][:-1])
                ind += 4
    return (barcodes, read1_seq, read2_seq)
        

def main():
    barcode_file = sys.argv[1]
    seq_file = sys.argv[2]
    seq2_file = sys.argv[3]
    pick_path = sys.argv[4]
    pick_prefix = sys.argv[5]

    seq_inf = parse_files(barcode_file, seq_file, seq2_file)
    with open(os.path.join(pick_path, pick_prefix) +".pi", "wb") as f:
        pickle.dump(seq_inf, f)

if __name__ == "__main__":
    main()