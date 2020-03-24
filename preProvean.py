#Function to read the .vcf file produced from Mutect2 and format it for use in PROVEAN
#
#File Header Reference: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TUMOR   NORMAL
#
runs = [<list of Sample IDs>]
#
for r in runs:
    infile = <path to run r>"/variants_PASS.vcf"
    outfile = <path to run r>"/prePROVEAN.txt"
    fi = open(infile, 'r')
    fo = open(outfile, 'w')
    print("Reading from " + infile)
    for i in fi:
        if i[0] != "#":
            i = i.split()
            chrom = i[0]
            chrom = chrom.lstrip("chr")
            pos = i[1]
            ID = i[2]
            ref = i[3]
            alt = i[4]
            qual = i[5]
            filt = i[6]
            info = i[7]
            fmt = i[8]
            tumor = i[9]
            normal = i[10]
            if filt == "PASS":
                fo.write(','.join([chrom, pos, ref, alt]) + '\n')
    fi.close()
    fo.close()
