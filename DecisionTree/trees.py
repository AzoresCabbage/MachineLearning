from math import log
import operator

def createDataSet():
    dataSet = [[1, 1, 'yes'],
            [1, 1, 'yes'],
            [1, 0, 'no'],
            [0, 1, 'no'],
            [0, 1, 'no']]
    labels = ['no surfacing','flippers']
    return dataSet, labels

#in py3.4 we must use 'wb' mode of file handle
def storeTree(inputTree,filename):
    import pickle
    fw = open(filename,'wb')
    pickle.dump(inputTree,fw)
    fw.close()
    
def grabTree(filename):
    import pickle
    fr = open(filename, "rb")
    return pickle.load(fr)    

def calEntropy(dataSet):
    numEntries = len(dataSet)
    labelCount = {}
    #count the occur times of every key
    for vec in dataSet:
        curLabel = vec[-1]
        if curLabel not in labelCount.keys():
            labelCount[curLabel] = 0
        labelCount[curLabel] += 1
    entropy = 0
    for key in labelCount:
        prob = float(labelCount[key])/numEntries
        entropy -= prob * log(prob, 2)
    return entropy

#return a lists of all list that list[pos] == val and this pos should be delete
def splitDataSet(dataSet, pos, val):
    ret = []
    for vec in dataSet:
        if vec[pos] == val:
            spvec = vec[:pos]
            spvec.extend(vec[pos+1:])
            ret.append(spvec)
    return ret
    
def chooseBestFeatureToSplit(dataSet):
    #assume the last column of dataSet is className
    featureNum = len(dataSet[0]) - 1
    baseEntropy = calEntropy(dataSet)
    idx = -1
    infoGain = 0
    for i in range(featureNum):
        entropy = 0
        for val in set(eg[i] for eg in dataSet):
            subSet = splitDataSet(dataSet, i, val)
            prob = float(len(subSet)) / len(dataSet)
            entropy += prob * calEntropy(subSet)
        entropy = baseEntropy - entropy
        if infoGain < entropy:
            infoGain = entropy
            idx = i
    return idx        

#return the most frequent key in classList
def majorityCnt(classList):
    classCount = {}
    for vote in  classList:
        if vote not in classCount.keys():
            classCount[vote] = 0;
        classCount[vote] += 1
    sortedList = sorted(classCount.item(), key = operator.iteritem.itemgetter(1), reverse = True)
    return sortedList[0][0]
    
def createTree(dataSet, labels):
    classList = [eg[-1] for eg in dataSet]
    #only have one kind of classtype,then this is a leaf node
    if classList.count(classList[0]) == len(classList):
        return classList[0]
    #only have one row in dataSet,then we can't split it and this is a leaf node
    if len(dataSet) == 1:
        return majorityCnt(classList)
    bestFeatIdx = chooseBestFeatureToSplit(dataSet)
    bestFeatLabel = labels[bestFeatIdx]
    myTree = {bestFeatLabel:{}}
    del labels[bestFeatIdx]
    for val in set(eg[bestFeatIdx] for eg in dataSet):
        sublabel = labels[:]
        myTree[bestFeatLabel][val] = createTree(splitDataSet(dataSet, bestFeatIdx, val), sublabel)
    return myTree
    
def classify(myTree, labels, testVec):
    label = [i for i in myTree][0]
    rootDic = myTree[label]
    idx = labels.index(label)
    for son in rootDic.keys():
        if son == testVec[idx]:
            if type(rootDic[son]).__name__ == 'dict':
                classifyName = classify(rootDic[son], labels, testVec)
            else:
                classifyName = rootDic[son]
    return classifyName
    
if __name__ == '__main__':
    fr = open('lenses.txt','r')
    lenses = [i.strip().split('\t') for i in fr.readlines()]
    lenseLabels = ['age','prescript','astigmatic','tearRate']
    lenseTree = createTree(lenses,lenseLabels)
    print(lenseTree)
