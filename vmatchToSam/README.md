## Utility to convert vmatch output to sam format

[vmatch](http://www.vmatch.de/) a versatile software tool for eﬃciently solving large scale sequence matching tasks.
`vmatchToSam` converts the vmatch output to [sam format](http://samtools.github.io/hts-specs/SAMv1.pdf) so that
various tools for sam analysis can be used on vmatch. 

## Mics

```
cd ~/svn_git/Java-cafe/trunk/vmatchToSam/test_data/
vmatch -showdesc 0 -s -p -d -complete -q reads.fa ref.fa
```
