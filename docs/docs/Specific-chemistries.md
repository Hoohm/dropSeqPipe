# Specific Chemistries

## SureCell/ddseq

This code has been provided by ([@TomKellyGenetics)](https://github.com/TomKellyGenetics)) in issue [#42](https://github.com/Hoohm/dropSeqPipe/issues/42)

You can run this on your files to get the structure that is expected by dropSeqPipe

```
Read1s=("Sample_S1_L001_R1_001.fastq" "Sample_S1_L002_R1_001.fastq")
Read2s=("Sample_S1_L001_R2_001.fastq" "Sample_S1_L002_R2_001.fastq")

    #remove adapter from SureCell (and correct phase blocks)
        for File in "${Read1s[@]}"; do
            #remove phase blocks and linkers
            sed -E '
                /.*(.{6})TAGCCATCGCATTGC(.{6})TACCTCTGAGCTGAA(.{6})ACG(.{8})GAC/ {
                s/.*(.{6})TAGCCATCGCATTGC(.{6})TACCTCTGAGCTGAA(.{6})ACG(.{8})GAC.*/\1\2\3\4/g
                n
                n
                s/.*(.{6}).{15}(.{6}).{15}(.{6}).{3}(.{8}).{3}/\1\2\3\4/g
                }' $File > .temp
            mv $.temp $File
        done
```