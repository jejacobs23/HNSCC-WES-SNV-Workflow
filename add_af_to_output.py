#This will take the provean results and add in the af from the .vcf file to the output.
#
#VCF File Header Reference: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TUMOR   NORMAL
#
runs = [<list of Sample IDs>]

for u in runs:
    infile1 = "variants_PASS.vcf"
    infile2 = "PROVEAN_RESULTS.tsv"
    outfile = "PROVEAN_RESULTS_with_AF.tsv"
    D = {}
    fi1 = open(infile1, 'r')
    fi2 = open(infile2, 'r')
    fo = open(outfile, 'w')
    print("Reading from " + infile1)
    print("Reading from " + infile2)
    for i in fi1: #read the .vcf file and create a dict. with the variant ID as the key and the af as the value
        if i[0] != "#":
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
            format = fmt.split(":")
            idx = format.index("AF")
            tumor_info = tumor.split(":")
            af = float(tumor_info[idx])
            unique = str(chrom) + ";" + str(pos) + ";" + ref + ":" + alt
            D[unique] = str(af)
    next(fi2) #This will skip the first line of the infile.  Only use if there is a header.
    fo.write('\t'.join(["INPUT", "PROTEIN_ID", "LENGTH", "STRAND", "CODON_CHANGE", "POS", "RESIDUE_REF", "RESIDUE_ALT", "AF_TUMOR", "TYPE", "SCORE", "PREDICTION_(cutoff=-2.5)", "SCORE", "PREDICTION_(cutoff=0.05)", "dbSNP_ID", "GENE_NAME", "DESCRIPTION"]) + '\n')
    for j in fi2: #add the af to the PROVEAN outpuf file
        j = j.split('\t')
        if len(j) != 21:
            print("invalid column number from input file: " + j[1])
        else:
            INPUT = j[1]
            PROTEIN_ID = j[2]
            LENGTH = j[3]
            STRAND = j[4]
            CODON_CHANGE = j[5]
            POS = j[6]
            RESIDUE_REF = j[7]
            RESIDUE_ALT = j[8]
            TYPE = j[9]
            pSCORE = j[10]
            pPREDICTION = j[11]
            sSCORE = j[14]
            sPREDICTION = j[15]
            dbSNP_ID = j[18]
            GENE_NAME = j[19]
            DESCRIPTION = j[20]
            sub = INPUT.split(",")
            UNIQUE = "chr" + str(sub[0]) + ";" + str(sub[1]) + ";" + sub[2] + ":" + sub[3]
            if UNIQUE in D:
                AF = D[UNIQUE]
                fo.write('\t'.join([INPUT, PROTEIN_ID, LENGTH, STRAND, CODON_CHANGE, POS, RESIDUE_REF, RESIDUE_ALT, AF, TYPE, pSCORE, pPREDICTION, sSCORE, sPREDICTION, dbSNP_ID, GENE_NAME, DESCRIPTION]))
            else:
                print("dictionary key-value problem :" + UNIQUE)
    fi1.close()
    fi2.close()
    fo.close()
