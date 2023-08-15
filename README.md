# SABER: Sequence Alignment using Block Edits and Rearrangements
_SABER_ is a pairwise sequence alignment algorithm under block edit distance models. It is currently capable of detecting block moves, reversals, and deletions together with the single character edits (insertion, deletion, substitution).    
_SABER_ can detect block (i.e., substring) rearrangement events by penalizing the same score for block operations and character operations and approximately finds the alignment that minimizes block edit distance between the two sequences.   
 
# Dependencies
* [SeqAn](https://www.seqan.de/): C++ library for sequence analysis. Please install SeqAn in your system by following the guide [here](https://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html#native-package-management).
* [Edlib](https://martinsos.github.io/edlib/): Lightweight library for sequence alignment. Edlib package is included as a submodule in _SABER_.
 
# Build

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
SABER only accepts FASTA and FASTQ files as inputs for source and target sequences.
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

# Rearrangement Simulator

This is the source code for the testing of _SABER_ over different intensity rates. After creating the _rearrangement\_sim_ executable using _make sim_ command, run the tests by:   

```
$ ./rearrangement_sim sequence.fa no-samples m-min m-max l-min l-max move-remove-rate max-iterations error-rate step-interval
```   

This testing code only accepts FASTA files. The following code is:   

***no-samples***        : Number of samples generated and tested for each intensity   
***m-min and m-max***   : Size range of each generated sample   
***l-min and l-max***   : Size range of the blocks in the block operations   
***move-remove-rate***  : Ratio of block move operations to remove operations in the simulation   
***max-iterations***    : Maximum number of iterations to test _saber_ with.   
***error-rate***        : Error rate for character edits.   
***step-interval***     : Step interval for intensity testing (the intensity starts from 10, increases by _step-interval_ each step)   

