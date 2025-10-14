# imports
import csv
import numpy as np
import pandas as pd
import sys
import math
import collections


def get_region_dict(file):
    """
    retrieve the 500 bp bins from text file. bins must be sorted.
    """

    regions_dict = collections.defaultdict(list)  

    with open(file, "r") as input_file:
        regions_file = csv.reader(input_file, delimiter="\t")

        for line in regions_file:
            chrom, start, end = line[0], int(line[1]), int(line[2])

            regions_dict[chrom].append({
                "start": start,
                "end": end,
                "meth_frac": np.zeros(5)
            })

    return regions_dict


def get_methylation_counts(file, regions_dict):
    """
    add together the methylation values for all CpGs in the selected region
    """

    # file of CpGs
    with open(file, "r") as input_file:
        cpg_file = csv.reader(input_file, delimiter="\t")
        linecounter = 0
        # get methylation read counts for each position
        for line in cpg_file:
            linecounter += 1
            chrom, start, end = line[0], int(line[1]), int(line[2])
            meth = np.array(line[3:], dtype=np.float64)
            binstart = 500*math.floor(start/500)
            binend = 500*math.ceil(end/500)
            startrange = range(binstart, binend, 500)
            endrange = range(binstart+499, binend+1, 500)
            # if (chrom, binstart, binend) in regions_dict:
            #     regions_dict[(chrom, binstart)] += meth
            for region in regions_dict[chrom]:
                if region["start"] < binstart:
                    continue
                elif region["end"] > binend:
                    break
                else:
                    region["meth_frac"] += meth
            if linecounter % 5000 == 0:
                print(f"Processed {linecounter} lines.")


    return regions_dict


def write_bed_file(output_file, regions_dict):
    """
    write bed file of summed counts for all tissues
    """
    with open(output_file, "w") as output:
        bed_file = csv.writer(output, delimiter="\t",  lineterminator="\n")

        for chrom in regions_dict:
            for region in regions_dict[chrom]:
                values = region["meth_frac"]
                bed_file.writerow(
                    [chrom] + [region["start"]] + [region["end"]] + list(values)
                )


if __name__ == "__main__":

    regions_file = sys.argv[1]
    tissue_cpg_file = sys.argv[2] # the input file
    output_file_name = sys.argv[3]
    print("Reading regions file...")
    regions = get_region_dict(
        regions_file
    )  # get dictionary of regions to sum within
    print("Summing methylation...")
    get_methylation_counts(
        tissue_cpg_file, regions
    )  # get methylation read counts of Cpgs within region
    print("Writing output...")
    write_bed_file(output_file_name, regions)  # write output
