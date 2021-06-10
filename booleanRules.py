#!/usr/bin/env python
"""booleanRules.py: Infers network from pseudotime and binary expression.
To run type python booleanRules.py gene expressionFile step orderFile maxAct maxRep networkFile threshold thresholdStep
gene: name of gene e.g. Gata1. This should match gene names in expression and network files.
expressionFile: path to expression file. This is a binary gene expression matrix with colnames equal to genes and rownames equal to cell names.
step: e.g. 5. This is how far back in pseudotime to take as input for the boolean function evaluated at current time t. So for step=5 the input to the boolean function at time t is given by the gene expression states at time t=t-5.
orderFile: path to pseudotime order file. This should be a list of cell names indicating the pseudotime ordering of cells. These cells names should be the same form as those in expressionFile.
maxAct: e.g. 4. The number of activators allowed for each gene
maxRep: e.g. 2. The number of repressors allowed for each gene. Note, if these two parameters are high finding rules will be very slow.
networkFile: path to the network file. This indicates possible activators and repressorrs of each gene. Each line should be of the form "from.gene relation to.gene" where relation is either 1 (activation) or -1 (repression).
threshold: Starting threshold for agreement of rules. 
thresholdStep: Increments in which the threshold will be lowered if no rules are found.
"""

__author__ = "Fiona Hamey"


import pandas as pd
import numpy as np
import math
import sys
sys.path.append('/Users/fiona/z3/build')
import z3
sys.path.append('Users/fiona/Desktop/PhD/MRes/Rotation_2/Code/Python/')


gene = sys.argv[1]
expressionFile = sys.argv[2]
step = int(sys.argv[3])
orderFile = sys.argv[4]
maxAct = int(sys.argv[5])
maxRep = int(sys.argv[6])
networkFile = sys.argv[7]
threshold = float(sys.argv[8])
thresholdStep = float(sys.argv[9])


# Read in expression matrix
expression = pd.read_table(expressionFile,  sep = "\t", 
                                 index_col = 0, header = 0)

# Read in path cells names as list
pathNames = list(pd.read_table(orderFile, sep = "\n").iloc[:,0])

from encodingFunctions import *

step = step
inputOutput = [(pathNames[i-step-1], pathNames[i-step], pathNames[i-step+1], pathNames[i]) for i in range(step+1,len(pathNames))]
inputOutputSelection = inputOutput

# Read in network 
network = pd.read_table(networkFile, sep = "\t", header = 0)

# Define node class
class node:
    def __init__(self, pos_in, neg_in):
        self.p = pos_in
        self.n = neg_in
        
        
# Node list for network
nodeList = {}
nodeNames = list(set(network['to.gene']))

for n in nodeNames:
    allRelations = network[network['to.gene'] == n]
    posRelations = allRelations[allRelations['relation'] == 1]
    posGenes = list(posRelations['from.gene'])
    if 'FoxH1' in posGenes:
        posGenes.remove('FoxH1')
    if 'FoxO4' in posGenes:
        posGenes.remove('FoxO4')
    negRelations = allRelations[allRelations['relation'] == -1]
    negGenes = list(negRelations['from.gene'])
    if 'FoxH1' in negGenes:
        negGenes.remove('FoxH1')
    if 'FoxO4' in negGenes:
        negGenes.remove('FoxO4')
    nodeList[n] = node(posGenes, negGenes)
    
# Need to add some penalties in, remember she allowed self-activation 
# Penalties
penAll = 0
penSelf = 0.005    
    
# Dictionary to look up node names number equivalent
allNames = list(expression.columns)
nodeLookUp = {allNames[j]:j+2 for j in range(len(allNames))}

# Add constraints to the solver
def constraintsBitVec(ctor, model, d):
    x = z3.BitVecVal(int(str(model[d])), 32)
    return ctor(str(d)) != x

def addConstraintsCircuitVar(solver, model, ds):
    constraints = z3.Or([constraintsBitVec(makeCircuitVar, model, d) for d in ds])
    return constraints
 
def modeCalc(lst):
    return max(set(lst), key=lst.count)
     
    
# Function which enforces minimum agreeing counts
def totalAgreement(inputOutputPair, gene, aVars, rVars, expressionValues, counter):
    
    inputName0 = inputOutputPair[0]
    inputName1 = inputOutputPair[1]
    inputName2 = inputOutputPair[2]
    outputName = inputOutputPair[3]
    input0 = expressionValues.loc[inputName0,:]
    input1 = expressionValues.loc[inputName1,:]
    input2 = expressionValues.loc[inputName2,:]
    input = [modeCalc([input0[i], input1[i], input2[i]]) for i in range(len(input0))]
    output = expressionValues.loc[outputName,:]
    #counter += 1
    score = makeEnforcedVar("counter_%i" %counter)
    (encoding, match) = circuitEvaluatesTo(gene, aVars, rVars, input, output, counter)

    return (z3.And(encoding, z3.If(match, score == z3.BitVecVal(1,32), score == z3.BitVecVal(0,32))), score)

    
def agreesAtLeastNTimes(inputOutputList, expressionValues, gene, aVars, rVars, n):  


    both = [totalAgreement(inputOutputList[p], gene, aVars, rVars, expressionValues, p) for p in range(len(inputOutputList))]
    encodings = [both[i][0] for i in range(len(both))]
    scoreValues = [both[i][1] for i in range(len(both))]
    
    return z3.And(z3.And(encodings), z3.Sum(scoreValues) >= n)



# Genes seems a bit pointless?
def findFunctions(solver, gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors, inputOutput, expression, agreementThreshold):
    
    expressionBool = expression == 1
    
    possAct = nodeList[gene].p
    possAct = [nodeLookUp[i] for i in possAct]
    possRep = nodeList[gene].n
    possRep = [nodeLookUp[i] for i in possRep]
    
    gene = nodeLookUp[gene]
    genes = [nodeLookUp[genes[i]] for i in range(len(genes))]
    
    circuitEncoding, aVars, rVars = encodeUpdateFunction(gene, genes, maxActivators, maxRepressors, possAct, possRep)
    
    circuitVars = aVars + rVars
    
    # Choose number of agreement threshold
    agreementThreshold = int(agreementThreshold*len(inputOutput))
    
    solver.reset()
    allConstraints = z3.And(circuitEncoding, agreesAtLeastNTimes(inputOutput, expressionBool, gene, aVars, rVars, \
    agreementThreshold))
    solver.add(allConstraints)
    
    constraints = True
    
    possibleRules = []
    while str(solver.check()) == 'sat':
        m = solver.model()
        modelVariables = [m[i] for i in range(len(m))]
        circuitDecls = filter(lambda x: str(x) in [str(s) for s in circuitVars], modelVariables)
        enforceDecls = filter(lambda x: str(x)[0:7] == "counter", modelVariables)
        totalScore = sum([int(str( m[d])) for d in enforceDecls])
        
        newConstraints = addConstraintsCircuitVar(solver, m, circuitDecls)
        constraints = z3.And(constraints, newConstraints)
        
        solver.reset()
        solver.add(z3.And(allConstraints, constraints))
        
        possibleRules.append((m, totalScore))
        print('found rule')

        if len(possibleRules) >= 100:
            newThreshold = max([s[1] for s in possibleRules])
            displayThreshold = float(newThreshold)/len(inputOutput)
            print ('Finding too many rules. So now increase threshold to %f' %displayThreshold)
            
            solver.reset()
            allConstraints = z3.And(circuitEncoding, agreesAtLeastNTimes(inputOutput, expressionBool, gene, aVars, rVars, \
                                    newThreshold))
            solver.add(allConstraints)
                        
            constraints = True
                                
            possibleRules = []
            
            while str(solver.check()) == 'sat':
                m = solver.model()
                modelVariables = [m[i] for i in range(len(m))]
                circuitDecls = filter(lambda x: str(x) in [str(s) for s in circuitVars], modelVariables)
                enforceDecls = filter(lambda x: str(x)[0:7] == "counter", modelVariables)
                totalScore = sum([int(str( m[d])) for d in enforceDecls])
   
                newConstraints = addConstraintsCircuitVar(solver, m, circuitDecls)
                constraints = z3.And(constraints, newConstraints)
    
                solver.reset()
                solver.add(z3.And(allConstraints, constraints))

                possibleRules.append((m, totalScore))
                print("found rule")

        
    return possibleRules, z3.And(allConstraints)
    
def findBestRuleForGene(gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors, inputOutput, expression):
    
    s = z3.Solver()
    
    for k in range(20):
        s.reset()
        agreement = threshold - thresholdStep*k
        print('Lowered threshold to %f' %agreement)
        rules, allConstraints = findFunctions(s, gene, genes, nodeLookUp, nodeList, maxActivators, maxRepressors, inputOutput, expression, agreement)
        convertedRules = [ruleConvert(r) for r in rules]
        geneInRules = []
        conRules = [x[0] for x in convertedRules]
        for r in conRules:
            geneInRules.append(s[1] for s in r)
        if len(rules) > 0 and not all(g in r for r in geneInRules):
            break
    
    print('Found %d rules for gene satisfying threshold %f' %(len(rules),agreement))
    return rules, agreement, allConstraints
    

def ruleConvert(model):
    
    score = model[1]/len(inputOutput)
    modelRules = model[0]

    def int2gene(i):
        i = int(str(i))
        if i == AND:
            return 'and'
        elif i == OR:
            return 'or'
        else:
            return allNames[i-2]
        
    namedVariables = ["a%i" %i for i in range(7)] + ["r%i" %i for i in range(7)]
    modelVariables = [modelRules[i] for i in range(len(modelRules))]
    usefulVariables = filter(lambda x: str(x) in namedVariables, modelVariables)
    modelConstraints = [(str(v), int2gene(modelRules[v])) for v in usefulVariables if str(modelRules[v]) != str(NOTHING)]
    
    return modelConstraints, score
    
    
def evaluateRule(rule, input0, input1, input2):

    # Remember to input boolean expression
    rules = {'a%d'%i : NOTHING for i in range(7)}
    rules.update({'r%d'%i : NOTHING for i in range(7)})
    
    convertBack = {'and':0, 'or':1}
    convertBack.update({gene:nodeLookUp[gene] for gene in allNames})

    for r in rule:
        rules[r[0]] = convertBack[r[1]]
    
    def getValue(v):
        return modeCalc([input0[v-2], input1[v-2], input2[v-2]])
    
    # Find intermediate variables
    v = {'va%d'%i: getValue(rules['a%d'%i]) for i in range(7) if rules['a%d'%i] in range(2,NOTHING)}
    v.update({'vr%d'%i: getValue(rules['r%d'%i]) for i in range(7) if rules['r%d'%i] in range(2,NOTHING)})
    
    # Evaluate activators
    if rules['a1'] == 0:
        inter_a1 = v['va3'] and v['va4']
    elif rules['a1'] == 1:
        inter_a1 = v['va3'] or v['va4']
    elif rules['a1'] in range(2,NOTHING):
        inter_a1 = v['va1']
    else:
        inter_a1 = True
        
    if rules['a2'] == 0:
        inter_a2 = v['va5'] and v['va6']
    elif rules['a2'] == 1:
        inter_a2 = v['va5'] or v['va6']
    elif rules['a2'] in range(2,NOTHING):
        inter_a2 = v['va2']
    else:
        inter_a2 = True
        
    if rules['a0'] == 0:
        inter_a0 = inter_a1 and inter_a2
    elif rules['a0'] == 1:
        inter_a0 = inter_a1 or inter_a2
    else:
        inter_a0 = v['va0']
    
    # Evaluate repressors    
    if rules['r0'] != NOTHING:
        if rules['r1'] == 0:
            inter_r1 = v['vr3'] and v['vr4']
        elif rules['r1'] == 1:
            inter_r1 = v['vr3'] or v['vr4']
        elif rules['r1'] in range(2,NOTHING):
            inter_r1 = v['vr1']
        else:
            inter_r1 = True
    
        if rules['r2'] == 0:
            inter_r2 = v['vr5'] and v['vr6']
        elif rules['r2'] == 1:
            inter_r2 = v['vr5'] or v['vr6']
        elif rules['r2'] in range(2,NOTHING):
            inter_r2 = v['vr2']
        else:
            inter_r2 = True
    
        if rules['r0'] == 0:
            inter_r0 = inter_r1 and inter_r2
        elif rules['r0'] == 1:
            inter_r0 = inter_r1 or inter_r2
        else:
            inter_r0 = v['vr0']
        
        return inter_a0 and not inter_r0
    else:
        return inter_a0
           
           
g = sys.argv[1]
print('Trying to find rules for %s' %g)
expressBool = expression == 1

geneRules, agreement, solver = findBestRuleForGene(g, allNames, nodeLookUp, nodeList, maxAct, maxRep, inputOutputSelection, expression)
rulesForGene = [ruleConvert(r) for r in geneRules]
print('Converted to readable format')

def scoreForRule(rule):
    raw = 0
    
    for io in range(len(inputOutput)):
        inputName0 = inputOutput[io][0]
        inputName1 = inputOutput[io][1]
        inputName2 = inputOutput[io][2]
        outputName = inputOutput[io][3]
        input0 = expressBool.loc[inputName0,:]
        input1 = expressBool.loc[inputName1,:]
        input2 = expressBool.loc[inputName2,:]
        output = expressBool.loc[outputName,g]
        predictedOutput = evaluateRule(rule, input0, input1, input2)
        
        if predictedOutput == output:
            raw += 1
            
    score = penalty = len(rule)*penAll
    if any(g in x for x in rule):
        penalty = penalty + (penSelf - penAll)
    score = raw-penalty
    
    # Might want to return more things?
    return (score, rule)
    
scoreRules = [(scoreForRule(r[0]),'z3 score %f' %r[1]) for r in rulesForGene]
print('Found agreement level for each rule')
bestRule = max(scoreRules, key=lambda x: x[0])
maxValue = bestRule[0][0]
allBest = []
for rule in scoreRules:
    if rule[0][0] == maxValue:
        allBest.append(rule)
    
print('The best rules are %s' %str(allBest))

scoreRules = sorted(scoreRules, key=lambda x: x[0][0], reverse = True)

print('Writing the rules to a file')
f=open('%s_boolean_rules_%d.txt' % (g, step),'w')
f.write('Agreement level = %f \n' %agreement)
f.write('The best rules were:\n')
for item in allBest:
    f.write("%s\n" %str(item))
f.write('Other rules found were:\n')
for item in scoreRules:
    f.write("%s\n" %str(item))
f.close()

print('Found rules for %s' %g)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
