#Function to read the .vcf file produced from Mutect2 and isolate only those variants with a "PASS" in their
#filter field.
#
#File Header Reference: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TUMOR   NORMAL
#
infile = "variants.vcf"
outfile = "variants_PASS.vcf"
fi = open(infile, 'r')
fo = open(outfile, 'w')
print("Reading from " + infile)
for i in fi:
    if i[0] == "#":
        fo.write(i)
    else:
        i = i.split()
        chrom = i[0]
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
            fo.write('\t'.join([chrom, pos, ID, ref, alt, qual, filt, info, fmt, tumor, normal]) + '\n')
fi.close()
fo.close()
