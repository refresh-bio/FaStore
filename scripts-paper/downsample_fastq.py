#!/usr/bin/python
import sys

quality_offset = 33
default_threshold = 20

def process_fq_illumina8(in_file, out_file):

    qua_trans_table =  [0, 0, 6, 6, 6,  6, 6, 6, 6, 6,  # 0x
                        15,15,15,15,15, 15,15,15,15,15, # 1x
                        22,22,22,22,22, 27,27,27,27,27, # 2x
                        33,33,33,33,33, 37,37,37,37,37, # 3x
                        40,40,40,40,40, 40,40,40,40,40, # 4x
                        40,40,40,40,40, 40,40,40,40,40, # 5x
                        40,40,40,40 ]

    with open(in_file, "r") as fq_in, open(out_file, "w") as fq_out:
        idx = 0
        for line in fq_in:

            if idx != 3:
                fq_out.write(line)
                idx += 1
                continue

            qua = line.strip()

            # bin q-scores
            new_q = []
            for q in qua:
                q = ord(q) - quality_offset
                qq = chr(qua_trans_table[q] + quality_offset)   # this can be done in separate fcn
                new_q.append(qq)

            fq_out.write(''.join(new_q) + '\n')
            idx = 0


def process_fq_threshold(in_file, out_file, threshold):
    min_value = 6
    max_value = 40

    with open(in_file, "r") as fq_in, open(out_file, "w") as fq_out:
        idx = 0
        for line in fq_in:

            # skip comments
            if idx != 3:
                fq_out.write(line)
                idx += 1
                continue

            qua = line.strip()
            
            # threshold q-scores
            new_q = []
            for q in qua:
                q = ord(q) - quality_offset
                qq = max_value if q >= threshold else min_value     # this can be done in separate fcn
                qq = chr(qq + quality_offset)
                new_q.append(qq)

            fq_out.write(''.join(new_q) + '\n')
            idx = 0



if __name__ == "__main__":
    if len(sys.argv) < 4 or (sys.argv[3] != 'B' and sys.argv[3] != 'T'):
        print >> sys.stderr, "usage: %s <in_fastq> <out_fastq> <B|T> [th]" % sys.argv[0]
        print >> sys.stderr, "where:  B - Illumina 8-level bining"
        print >> sys.stderr, "        T - binary thresholding using 'th' value (th default: %d)" % default_threshold
        exit()

    print >> sys.stderr, "processing..."
    if sys.argv[3] == 'B':
        process_fq_illumina8(sys.argv[1], sys.argv[2])
    else:
        threshold = default_threshold
        if len(sys.argv) == 5:
            threshold = int(sys.argv[4])
        process_fq_threshold(sys.argv[1], sys.argv[2], threshold)

    print >> sys.stderr, "done!"
