import os, re, sys
import argparse
from argparse import RawTextHelpFormatter
import pandas as pd


def getArgs():
    parser = argparse.ArgumentParser(
        description="Convert TCGA RNA count value to TPM",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument("-i", "--infile", help="TCGA HiSeqV2 count file")
    parser.add_argument("-o", "--outdir", help="outdir")
    parser.add_argument(
        "-id",
        "--id",
        default="Ensembl_ID",
        help="id [gene_name or Ensembl_ID],default=Ensembl_ID",
    )
    parser.add_argument(
        "-gtf", "--gtffile", help="gtf file,if have geneLengthfile,dont need this file"
    )
    parser.add_argument(
        "-GC",
        "--geneColumn",
        default=4,
        type=int,
        help="genename column in gtf,In general, it is the third column in the emsemble GTF file and the fourth column in the genecode file,default=4",
    )
    parser.add_argument(
        "-g", "--geneLengthfile", help="Gene longest transcript length file"
    )
    parser.add_argument(
        "-outid",
        "--outid",
        default="Ensembl_ID",
        help="outfile id type[gene_name or Ensembl_ID],default=Ensembl_ID",
    )
    argv = parser.parse_args()
    return argv


def Outfile(outdir, name):
    outfile = os.path.join(outdir, name)
    f = open(outfile, "w")
    return (outfile, f)


def Gtf(gtffile, outdir, id, geneColumn):
    dic = {}
    name_id = {}
    with open(gtffile, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            lst = line.rstrip().split("\t")
            if lst[2] == "exon":
                transcript = re.search(r'["](.*?)["]', lst[8].split(";")[1]).group(1)
                genename = re.search(
                    r'["](.*?)["]', lst[8].split(";")[geneColumn]
                ).group(1)
                ensemble = re.search(r'["](.*?)["]', lst[8].split(";")[0]).group(1)
                if id == "gene_name":
                    gene = genename
                    name_id[genename] = ensemble
                elif id == "Ensembl_ID":
                    gene = ensemble
                    name_id[ensemble] = genename
                start = lst[3]
                end = lst[4]
                length = int(end) - int(start) + 1
                if gene not in dic:
                    dic[gene] = {}
                if transcript not in dic[gene]:
                    dic[gene][transcript] = 0
                dic[gene][transcript] += length
    outfile1, f1 = Outfile(outdir, id + ".gene_transcript_exonlength.xls")
    outfile2, f2 = Outfile(outdir, id + ".gene_longest_exonlength.xls")
    longestdic = {}
    for i in dic:
        longest = 0
        for t in dic[i]:
            l = dic[i][t]
            f1.write("\t".join([i, name_id[i], t, str(l)]) + "\n")
            if l > longest:
                longest = l
        longestdic[i] = longest
        f2.write("\t".join([i, name_id[i], str(longest)]) + "\n")
    f1.close()
    f2.close()
    return (longestdic, name_id)


def Longest(geneLengthfile):
    dic = {}
    name_id = {}
    with open(geneLengthfile, "r") as f:
        for line in f:
            lst = line.rstrip().split("\t")
            dic[lst[1]] = int(lst[2])
            name_id[lst[0]] = lst[1]
    return (dic, name_id)


def Count_TPM(infile, outdir, gtffile, geneLengthfile, id, outid, geneColumn):
    infile = os.path.abspath(infile)
    outdir = os.path.abspath(outdir)
    outfile = os.path.join(outdir, infile.rsplit("/", 1)[1] + "TPM.xls")
    if geneLengthfile is None:
        dic, name_id = Gtf(gtffile, outdir, id, geneColumn)
    else:
        dic, name_id = Longest(geneLengthfile)
    dic_adjust = {}
    dicsum = {}
    with open(infile, "r") as f:
        for line in f:
            lst = line.rstrip().split("\t")
            if (
                line.startswith("sample")
                or line.startswith("Ensembl_ID")
                or line.startswith("Gene")
            ):
                sample = lst[1:]
                num = len(sample)
                continue
            gene = lst[0]
            count = lst[1:]

            if id=="gene_name":
                pass
            else:
                if gene not in name_id:
                    continue
                else:
                    gene=name_id[gene]

            if gene not in dic:
                print(gene)
                continue
            exonlength = dic[gene]
            #print(gene)
            for i in range(0, num):
                c = count[i]
                s = sample[i]
                ad = int(c) / exonlength
                if s not in dic_adjust:
                    dic_adjust[s] = {}
                    dicsum[s] = 0
                if gene not in dic_adjust[s]:
                    dic_adjust[s][gene] = ad
                    dicsum[s] += ad
    dicTPM = {}
    for s in dic_adjust:
        dicTPM[s] = {}
        for g in dic_adjust[s]:
            # if (
            #     id == "Ensembl_ID"
            #     and outid == "gene_name"
            #     or id == "gene_name"
            #     and outid == "Ensembl_ID"
            # ):
            #     gene = name_id[g]
            dicTPM[s][g] = dic_adjust[s][g] * 1000000 / dicsum[s]
    df = pd.DataFrame.from_dict(dicTPM)
    df.to_csv(outfile, sep="\t")
    # print(dic_adjust)


if __name__ == "__main__":
    argv = getArgs()
    Count_TPM(
        argv.infile,
        argv.outdir,
        argv.gtffile,
        argv.geneLengthfile,
        argv.id,
        argv.outid,
        argv.geneColumn,
    )
