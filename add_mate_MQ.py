#!/usr/bin/env python3
"""
Add mate's MQ tags to SAM file.
"""

import sys


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    #parser.add_argument("-bam", required=True)
    args = parser.parse_args()

    current_qname = None
    read_cache = dict()

    def get_flags(line):
        toks = line.split()
        return int(toks[1])

    def key(line):
        toks = line.split()
        return (toks[0], toks[2], toks[3])

    def mate_key(line):
        toks = line.split()
        flags = get_flags(line)
        if flags & 0x8:
            # mate is unmapped
            return (toks[0], "*", "*")
        else:
            # mate is mapped!
            if toks[6] == "=": toks[6] = toks[2]
            return (toks[0], toks[6], toks[7])

    for line in sys.stdin:
        line = line.strip()

        if line.startswith("@"):
            print(line)
            continue

        qkey = key(line)
        qname = qkey[0]

        if current_qname is None: current_qname = qname

        if current_qname != qname:
            # time to print this read
            #print(">>>>>>", qname, current_qname)
            #print(list(read_cache.keys()))
            #print("-------------------")

            for _qkey, _line in read_cache.items():
                mkey = mate_key(_line)

                first_read = read_cache[_qkey]
                second_read = read_cache[mkey]

                toks = second_read.split()
                second_mq = int(toks[4])

                print(first_read + "\t" + f"MQ:i:{second_mq}")
            read_cache = dict()

        flags = get_flags(line)
        if flags & 0x4:
            # read is unmapped.
            # Note: BWA prints coords on unmapped read records 🤯
            qkey = (qkey[0], "*", "*")
        read_cache[qkey] = line
        current_qname = qkey[0]
