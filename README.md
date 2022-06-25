# SABER: Sequence Alignment using Block Edits and Rearrangements
_SABER_ is a pairwise sequence alignment algorithm under block edit distance models. It is currently capable of detecting block moves, reversals, and deletions together with the single character edits (insertion, deletion, substitution).    
_SABER_ can detect block (i.e., substring) rearrangement events by penalizing the same score for block operations and character operations and approximately finds the alignment that minimizes block edit distance between the two sequences.   
 
# Build
_SABER_ currently uses SeqAn and Edlib libraries to run.      
Make sure to download SeqAn natively by following the guide [here](https://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html#native-package-management) before using _SABER_.    

To download and build _SABER_, run the following commands:     

```
$ git clone --recurse-submodules https://github.com/BilkentCompGen/saber.git
$ cd saber/src    
$ make   
```
This will create a _saber_ executable.

# Run

After creating the _saber_ executable, run _SABER_ on the target and source sequences with:   
```
$ ./saber -s source.fa -t target.fa [-optional arguments]
```
SABER only accepts fasta and fastq files as inputs for source and target sequences.
Some useful optional arguments are as follows:    

***-h or --help:*** Display help menu      

***-r or --runtime*** : Display the runtime of the program.     
***-o "filename"***   : Specify the output path (default: stdout)     
***-i _integer_***    : Specify the number of iterations in the algorithm (default: 3)      
***-l _integer_***    : Specify the minimum block length (default: 8)       
***-m _integer_***    : Specify the maximum block length (default: 15)       
***-e _float_***      : Specify the error rate (default: 0.3)        

For more detailed information, run
```
$ ./saber --help
```

