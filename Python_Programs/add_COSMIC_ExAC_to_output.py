#This will take the provean results and add in an COSMIC annotations if they exist.
#
#VCF File Header Reference: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TUMOR   NORMAL
#
runs = [<list of Sample IDs>]

import gzip
infile2 = "CosmicCodingMuts.vcf"
fi2 = open(infile2, 'r')
D = {}
print("Reading from " + infile2)
for i in fi2: #read the .vcf file and create a dict. with the variant ID as the key and the COSMIC info as the value
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
            info = info.split(";")
            for j in range(0,len(info)):
                if info[j][0:3] == "CDS":
                    C1 = info[j]
                if info[j][0:2] == "AA":
                    C2 = info[j]
            COSMIC = C1 + ";" + C2
            unique = str(chrom) + ";" + str(pos) + ";" + ref + ":" + alt
            D[unique] = COSMIC
fi2.close()

infile3 = "ExAC.r1.sites.vep.vcf.gz" #ExAC database
fi3 = gzip.open(infile3, 'r')
print("Reading from " + infile3)
E = {}
for e in fi3:
    if e[0] != "#":
        e = e.split()
        echrom = e[0]
        epos = e[1]
        eID = e[2]
        eref = e[3]
        ealt = e[4]
        equal = e[5]
        efilt = e[6]
        einfo = e[7]
        einfo = einfo.split(";")
        for h in range(0,len(einfo)):
            if einfo[h][0:2] == "AF":
                af = einfo[h].split("=")
                af = str(af[1])
        eunique = "chr" + str(echrom) + ";" + str(epos) + ";" + eref + ":" + ealt
        E[eunique] = af
fi3.close()

for u in runs:
    infile1 = <path to run u>"/PROVEAN_RESULTS_with_AF.tsv"
    outfile = <path to run u>"/PROVEAN_RESULTS_with_AF_COSMIC_ExAC.tsv"
    fi1 = open(infile1, 'r')
    fo = open(outfile, 'w')
    print("Reading from " + infile1)
    next(fi1) #This will skip the first line of the infile.  Only use if there is a header.
    fo.write('\t'.join(["INPUT", "PROTEIN_ID", "LENGTH", "STRAND", "CODON_CHANGE", "POS", "RESIDUE_REF", "RESIDUE_ALT", "AF_TUMOR", "TYPE", "SCORE", "PREDICTION_(cutoff=-2.5)", "SCORE", "PREDICTION_(cutoff=0.05)", "dbSNP_ID", "GENE_NAME$
    for j in fi1: #add the COSMIC mutations and ExAC allele frequency to the PROVEAN outpuf file
        j = j.split('\t')
        if len(j) != 17:
            print("invalid column number from input file: " + j[1])
        else:
            COSMIC = ""
            ExAC = ""
            INPUT = j[0]
            PROTEIN_ID = j[1]
            LENGTH = j[2]
            STRAND = j[3]
            CODON_CHANGE = j[4]
            POS = j[5]
            RESIDUE_REF = j[6]
            RESIDUE_ALT = j[7]
            AF_TUMOR = j[8]
            TYPE = j[9]
            pSCORE = j[10]
            pPREDICTION = j[11]
            sSCORE = j[12]
            sPREDICTION = j[13]
            dbSNP_ID = j[14]
            GENE_NAME = j[15]
            DESCRIPTION = j[16].rstrip()
            sub = INPUT.split(",")
            UNIQUE = "chr" + str(sub[0]) + ";" + str(sub[1]) + ";" + sub[2] + ":" + sub[3]
            if UNIQUE in D:
                COSMIC = D[UNIQUE]
            if UNIQUE in E:
                ExAC = E[UNIQUE]
            fo.write('\t'.join([INPUT, PROTEIN_ID, LENGTH, STRAND, CODON_CHANGE, POS, RESIDUE_REF, RESIDUE_ALT, AF_TUMOR, TYPE, pSCORE, pPREDICTION, sSCORE, sPREDICTION, dbSNP_ID, GENE_NAME, DESCRIPTION, COSMIC, ExAC]) + '\n')
    fi1.close()
    fo.close()
