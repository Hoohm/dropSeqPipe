# FAQ

## 1. I get `error='Cannot allocate memory' (errno=12)`, what should I do. [Fixed]

This has been fixed by using a wrapper exposing the TMPDIR to the pipeline.

First, be sure that your TMPDIR from the first configuration yaml has at least 100Go.
If you still have problems, you should edit the following files in the Drop-seq_tools-1.12:

* TagBamWithReadSequenceExtended
* FilterBAM
* TrimStartingSequence
* PolyATrimmer
* TagReadWithGeneExon
* DetectBeadSynthesisErrors
* SingleCellRnaSeqMetricsCollector
* BAMTagHistogram

In each of those files, the last line should be something like:
`java -Xmx${xmx} -Djava.io.tmpdir=/path/to/temp/folder/ -jar $jar_deploy_dir/dropseq.jar $progname $*`

You can also use this simple bash script to do it:  
Replace `/path/to/temp/folder/` with your temp path and don't forget to use escapes for /
```
for f in  BAMTagHistogram SingleCellRnaSeqMetricsCollector DetectBeadSynthesisErrors TagReadWithGeneExon PolyATrimmer TrimStartingSequence FilterBAM TagBamWithReadSequenceExtended
do
 sed -i 's/java -Xmx${xmx}/java -Xmx${xmx} -Djava.io.tmpdir=/path/to/temp/folder/ /g' $f
done
```