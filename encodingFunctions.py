#!/usr/bin/env python

"""encodingFunctions.py: Provides encoding functions for finding boolean rules governing gene expression with booleanRules.py script."""

__author__ = "Fiona Hamey"

import pandas as pd
import sys
sys.path.append('/Users/fiona/z3/build')
import z3

expressionFile = sys.argv[2]
networkFile = sys.argv[7]

# Read in expression matrix
expression = pd.read_table(expressionFile,  sep = "\t",
                                         index_col = 0, header = 0)

# Read in network 
network = pd.read_table(networkFile, sep = "\t", header = 0)
  
# How many genes are there
allGenes = list(expression.columns)
numGenes = len(allGenes)

# Define numbers corresponding to gates, genes and nothing variable
AND = 0
OR = 1
NOTHING = numGenes + 2

# Circuits are encoded as bitvectors
def makeCircuitVar(name):
    return z3.BitVec(name, 32)

    
def makeEnforcedVar(name):
    return z3.BitVec(name, 32)


# For each element in bitvector which gene/gate values is it allowed to take
def variableDomains(var, booleanAllowed, possibleInputs):
    
    if booleanAllowed == True:
        allowedValues = [0, 1] + possibleInputs
    else:
        allowedValues = possibleInputs
    
    return z3.Or([var == allowedValues[i] for i in range(len(allowedValues))])
    
    
# If a node = nothing then it is not allowed a parent as a gate
def parentsOfNothingArentGates(a, r):
    
    def f(c1,c2,p):
        return z3.Implies(z3.Or((c1 == NOTHING), (c2 == NOTHING)), z3.And(p != AND, p != OR))
    
    
    aParents = z3.And(z3.Implies(z3.Or(a[1] == NOTHING, a[2] == NOTHING), z3.And(a[0] != AND, a[0] != OR, a[0] != NOTHING)), \
                f(a[3], a[4], a[1]), \
                f(a[5], a[6], a[2]))
                
    rParents = z3.And(f(r[1], r[2], r[0]), \
                f(r[3], r[4], r[1]), \
                f(r[5], r[6], r[2]))   
    
    return z3.And(aParents, rParents)
    

# If a node is a gene then it must have a gate as a parent    
def parentsOfRestAreGates(a, r):
    
    def f(c1,c2,p):
        return z3.Implies(z3.Or((c1 != NOTHING), (c2 != NOTHING)), z3.Or(p == AND, p == OR))
    
    aParents = z3.And(f(a[1], a[2], a[0]), \
                f(a[3], a[4], a[1]), \
                f(a[5], a[6], a[2]))
                
    rParents = z3.And(f(r[1], r[2], r[0]), \
                f(r[3], r[4], r[1]), \
                f(r[5], r[6], r[2]))
                
    return z3.And(aParents, rParents)
    

# Can't have a gene more than once in the relation
def variablesDoNotAppearMoreThanOnce(symVars):
    
    def isVar(v):
        return z3.And(v != NOTHING, v != AND, v != OR)
    
    def notEqual(v, vars):
        return z3.And([v != i for i in vars if not z3.eq(v,i)])
        
    def doesNotAppearMoreThanOnce(v):
        return z3.Implies(isVar(v), notEqual(v, symVars))
        
    return z3.And([ doesNotAppearMoreThanOnce(j) for j in symVars])
    

# Don't waste time doing things multiple times
def enforceSiblingLexigraphicalOrdering(v1, v2):
    return (v1 <= v2)
    
    
def enforceLexigraphicalOrderingBetweenBranches(p1, p2, c1, c2):
    return z3.Implies(p1 == p2, c1 <= c2)
    
    
def enforceLexigraphicalOrderingNaryGate(vars):
    return z3.Implies(vars[0] == vars[1], vars[2] <= vars[3])


# Store the activator and repressor variables in a list    
activatorVars = ["a" + str(i) for i in range(7)]
repressorVars = ["r" + str(i) for i in range(7)]
circuitVars = activatorVars + repressorVars

# Depending on maximum number of inputs may want fewer nodes
def fixMaxInputs(v, max):
    if max == 0:
        return makeCircuitVar(v+"0") == NOTHING
    elif max == 1:
        return makeCircuitVar(v+"2") == NOTHING
    elif max == 2:
        return makeCircuitVar(v+"4") == NOTHING
    elif max == 3:
        return makeCircuitVar(v+"6") == NOTHING
    else:
        return True
        
def fixMaxActivators(max):
    return fixMaxInputs("a", max)

def fixMaxRepressors(max):
    return fixMaxInputs("r", max)


# This encodes the allowed update functions for a gene
def encodeUpdateFunction(gene, genes, maxActivators, maxRepressors, possAct, possRep):
    
    # Check all inputs are of right form
    assert (gene in genes and maxActivators > 0 and maxActivators <= 4 and maxRepressors >= 0 and maxRepressors <= 4), \
        "Incorrect arguments to encodeUpdateFunction"
        
    a = [makeCircuitVar("a%i" %i) for i in range(7)]
    r = [makeCircuitVar("r%i" %i) for i in range(7)]
    
    circuitEncoding = z3.And(variableDomains(a[0], True, possAct), \
                            variableDomains(r[0], True, possRep + [NOTHING]), \
                            
                            variableDomains(a[1], True, possAct + [NOTHING]), variableDomains(a[2], True, possAct + [NOTHING]), \
                            variableDomains(r[1], True, possRep + [NOTHING]), variableDomains(r[2], True, possRep + [NOTHING]), \
                            
                            variableDomains(a[3], False, possAct+ [NOTHING]), variableDomains(a[4], False, possAct + [NOTHING]), \
                            variableDomains(a[5], False, possAct + [NOTHING]), variableDomains(a[6], False, possAct + [NOTHING]), \
                            variableDomains(r[3], False, possRep + [NOTHING]), variableDomains(r[4], False, possRep + [NOTHING]), \
                            variableDomains(r[5], False, possRep + [NOTHING]), variableDomains(r[6], False, possRep + [NOTHING]), \
                            
                            parentsOfNothingArentGates(a, r), \
                            parentsOfRestAreGates(a, r), \
                            variablesDoNotAppearMoreThanOnce(a + r), \
                            
                            z3.And([enforceSiblingLexigraphicalOrdering(a[i], a[i+1]) for i in [1,3,5]]), \
                            z3.And([enforceSiblingLexigraphicalOrdering(r[i], r[i+1]) for i in [1,3,5]]), \
                            
                            enforceLexigraphicalOrderingBetweenBranches(a[1], a[2], a[3], a[5]), \
                            enforceLexigraphicalOrderingBetweenBranches(r[1], r[2], r[3], r[5]), \
                            
                            enforceLexigraphicalOrderingNaryGate(a), \
                            enforceLexigraphicalOrderingNaryGate(r), \
                            
                            fixMaxActivators(maxActivators), \
                            fixMaxRepressors(maxRepressors))
                            
    return (circuitEncoding, a, r)
    
    
# Given inputs evaluates function
def evaluateUpdateFunction(aVars, rVars, geneValues, counter):
    i = counter
    
    intermediateValueVariablesA = [ z3.Bool("va%i_%i" % (j, i)) for j in range(7)]
    intermediateValueVariablesR = [ z3.Bool("vr%i_%i" % (j, i)) for j in range(7)]
    
    def andConstraints(symVars, variables, pi, c1i, c2i):
        return z3.Implies(symVars[pi] == z3.BitVecVal(AND,32), variables[pi] == z3.And(variables[c1i], variables[c2i]))
    
    def orConstraints(symVars, variables, pi, c1i, c2i):
        return z3.Implies(symVars[pi] == z3.BitVecVal(OR,32), variables[pi] == z3.Or(variables[c1i], variables[c2i]))
        
        
    def variableConstraints(symVars, intermediateVars):
        
        def f(symVar, i):
            return z3.And([z3.Implies(symVar == v, intermediateVars[i] == z3.BoolVal(geneValues[v-2])) for v in range(2, NOTHING)])
            
        return z3.And([f(symVars[i], i) for i in range(7)])
        
    
    circuitVal = z3.Bool("circuit_%i" % i)

    def circuitValue():
        noRepressors = rVars[0] == NOTHING
        
        return z3.If(noRepressors, intermediateValueVariablesA[0], \
        z3.And(intermediateValueVariablesA[0],  z3.Not(intermediateValueVariablesR[0])))
          
    
    return (z3.And([variableConstraints(aVars, intermediateValueVariablesA), \
                    variableConstraints(rVars, intermediateValueVariablesR), \
                    andConstraints(aVars, intermediateValueVariablesA, 0, 1, 2), \
                    andConstraints(aVars, intermediateValueVariablesA, 1, 3, 4), \
                    andConstraints(aVars, intermediateValueVariablesA, 2, 5, 6), \
                    andConstraints(rVars, intermediateValueVariablesR, 0, 1, 2), \
                    andConstraints(rVars, intermediateValueVariablesR, 1, 3, 4), \
                    andConstraints(rVars, intermediateValueVariablesR, 2, 5, 6), \
                    orConstraints(aVars, intermediateValueVariablesA, 0, 1, 2), \
                    orConstraints(aVars, intermediateValueVariablesA, 1, 3, 4), \
                    orConstraints(aVars, intermediateValueVariablesA, 2, 5, 6), \
                    orConstraints(rVars, intermediateValueVariablesR, 0, 1, 2), \
                    orConstraints(rVars, intermediateValueVariablesR, 1, 3, 4), \
                    orConstraints(rVars, intermediateValueVariablesR, 2, 5, 6), \
                    circuitVal == circuitValue()]), circuitVal)
                    
    
def circuitEvaluatesTo(gene, aVars, rVars, input, output, counter):
    outValue = output[gene-2]
    inValues = input
    
    evaluationEncoding, circuitVal = evaluateUpdateFunction(aVars, rVars, inValues, counter)
    return (evaluationEncoding, circuitVal == z3.BoolVal(outValue))
    
