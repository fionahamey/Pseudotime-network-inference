# Pseudotime-network-inference
These scripts are provided for inferring Boolean regulatory rules from single-cell expression profiles. The algorithm has been designed for use 
with qRT-PCR data but could also be applied to other types of single-cell data measuring gene or protein expression. 
The user should provide an initial list of possible activators and repressors of each gene, a binary gene expression matrix,
and a pseudotime ordering of cells along a differentiation trajectory.

## Network inference scripts
Two python scripts are provided for applying the network inference algorithm. Python can be downloaded from https://www.python.org/.
These scripts use the Microsoft Z3 theorem prover to identify
high scoring Boolean rules. Z3 can be downloaded from https://github.com/Z3Prover/z3 and installed following the instructions provided at
this link. 

The script encodingFunctions.py provides Z3 encoding for the form and scoring of Boolean functions. Users can run the script
booleanRules.py to infer Boolean rules for a chosen gene. Users should change the line
```python
sys.path.append('/Users/fiona/z3/build')
```
to the path of the Z3 build on their own machine in both booleanRules.py and encodingFunctions.py. Additionally the line 
```python
sys.path.append('Users/fiona/Desktop/PhD/MRes/Rotation_2/Code/Python/')
```
in booleanRules.py should be changed to the path of the network inference scripts.


To run this rule inference on the command line use the following command.

```python booleanRules.py gene expressionFile step orderFile maxAct maxRep networkFile threshold thresholdStep```

* **gene** = Name of gene for which rules are to be inferred, e.g Gata1.
* **expressionFile** = Path to the expression file. This is a binary gene expression matrix with column names equal to genes and 
row names equal to cell names.
* **step** = An integer, e.g. 5. This is how far back in pseudotime to take as input for the boolean function evaluated at current time t. 
So for step=5 the input to the boolean function at time t is given by the gene expression states at time t=t-5.
* **orderFile** = Path to pseudotime ordering file. This should be a list of cell names indicating the pseudotime ordering of cells. 
These cells names should match those in expressionFile.
* **maxAct** = An integer from the set {0,1,2,3,4}. The number of activators allowed for each gene. Increasing this parameter can produce
more complex Boolean rules but will increase computation time.
* **maxRep** = An integer from the set {0,1,2,3,4}. The number of repressors allowed for each gene. Increasing this parameter can produce
more complex Boolean rules but will increase computation time.
* **networkFile** = Path to the network file. This indicates possible activators and repressors of each gene. Each line should be of the 
form "from.gene relation to.gene" where the relation is either 1 (indicating activation) or -1 (indicating repression).
* **threshold** = A decimal in the range [0,1] giving the starting threshold for agreement of rules. 
* **thresholdStep** = A decimal in the range [0,1] giving the increment in which the threshold will be lowered if no rules are found.

The rules found will be written to the file gene_boolean_rules_step.txt along with their score.

## Example script
Running this line of code would identify regulatory rules for *Bptf* in the MEP network. 

```python booleanRules.py Bptf binary_expression_MEP.txt 5 MEP_trajectory_order.txt 2 2 startingNetworkParCor.txt 0.95 0.05```

## Acknowledgements
I would like to thank Steven Woodhouse for help in writing the encoding of the Boolean functions - his scripts from the 
[Single Cell Network Synthesis Toolkit](https://github.com/swoodhouse/SCNS-Toolkit) were an extremely useful starting point.

## License
Scripts for the pseudotime network inference algorithm are released under the Apache 2.0 license; see LICENSE for more details.
