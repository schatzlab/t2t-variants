import argparse
import sys
import cyvcf2
import pandas as pd

default_columns = [
    "SVTYPe"
]


def vcf_record2dict(vcf_record: cyvcf2.Variant, columns=None):
    result = {}
    result["CHR1"] = vcf_record.CHROM
    result["CHR2"] = vcf_record.INFO["CHR2"]
    result["START"] = vcf_record.POS
    result["END"] = vcf_record.INFO["END"]
    result["SIZE"] = vcf_record.INFO["SVLEN"]
    result["SVTYPE"] = vcf_record.INFO["SVTYPE"]
    result["OR_SVTYPE"] = vcf_record.INFO.get("OLDTYPE", result["SVTYPE"])
    result["RE"] = vcf_record.INFO["RE"]
    result["SPECIFIC"] = vcf_record.INFO["IS_SPECIFIC"]
    result["PRECISE"] = vcf_record.INFO.get("PRECISE", False)
    result["STRANDS"] = vcf_record.INFO["STRANDS"]
    result["RE"] = vcf_record.INFO["RE"]
    result["SUPP_VEC"] = vcf_record.INFO.get("SUPP_VEC", "")
    gts = vcf_record.genotypes[0]
    result["GT"] = ("|" if gts[2] else "/").join(str(x) for x in gts[:2])
    return result


def get_vcf_df(vcf, sample=None, ref=None):
    entries = []
    for record in vcf:
        entries.append(vcf_record2dict(record))
    result = pd.DataFrame.from_records(entries)
    if sample is not None:
        result["sample"] = sample
    if ref is not None:
        result["ref"] = ref
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("VCF")
    parser.add_argument("--sample", default=None)
    parser.add_argument("--ref", default=None)
    parser.add_argument("--o-index", action="store_true", dest="output_index")
    parser.add_argument("-o", "--output", default=sys.stdout)
    args = parser.parse_args()
    vcf = cyvcf2.VCF(args.VCF)
    df = get_vcf_df(vcf, args.sample, args.ref)
    df.to_csv(args.output, index=args.output_index)

