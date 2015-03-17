from numpy import *
import time
import operator
import os

def createDataSet():
    group = array([[1.0,1.1],[1.0,1.0],[0,0],[0,0.1]])
    labels = ['A','A','B','B']
    return group,labels

def classify(inX, dataSet, labels, k):
    #the number of row in dataSet
    rowNum = dataSet.shape[0]
    
    #tile func to repeat inX rowNum row to be a matrix
    #then can get diff between inX and every row in dataSet via matrix minus
    diffMat = tile(inX, (rowNum,1)) - dataSet
    sqDiff = diffMat**2
    
    #get sum by row
    sqDist = sqDiff.sum(axis = 1)
    dist = sqDist**0.5
    
    #argsort sort the list and return the idx in list after sort them
    sortedDistIdx = dist.argsort()
    classCount = {}
    for i in range(k):
        #statictic the label
        vote = labels[sortedDistIdx[i]]
        classCount[vote] = classCount.get(vote,0)+1
    #sort the dict and get the label that appears most times
    sortedClassCount = sorted(classCount.items(), key=operator.itemgetter(1), reverse = True)
    return sortedClassCount[0][0]

#convert handwriting matrix(img) to vector
def img2vector(filename):
    row = 32
    col = 32
    imgVector = zeros((1, row*col))
    with open(filename) as f:
        for i in range(row):
            line = f.readline()
            for j in range(col):
                imgVector[0, j+col*i] = int(line[j])
    return imgVector
            
def loadData():
    print('get training set...')
    mainDir = r'E:\research\MachineLearning\Machine Learning in Action'
    SampleDir = mainDir + r'\digits\trainingDigits'
    SampleList = os.listdir(SampleDir)
    numSample = len(SampleList)
    sample_x = zeros((numSample, 1024))
    sample_y = []
    for i in range(numSample):
        filename = SampleList[i]
        sample_x[i, ] = img2vector(SampleDir+'\%s' % filename)
        label = filename.split('_')[0]
        sample_y.append(label)
    
    print('get test set...')
    testDir = mainDir + r'\digits\testDigits'
    testList = os.listdir(testDir)
    numTest = len(testList)
    test_x = zeros((numTest, 1024))
    test_y = []
    for i in range(numTest):
        filename = testList[i]
        test_x[i, ] = img2vector(testDir+'\%s' % filename)
        label = filename.split('_')[0]
        test_y.append(label)
    return sample_x, sample_y, test_x, test_y

def testHandWritingNum():
    print('1.Loading data...')
    sample_x, sample_y, test_x, test_y = loadData()
    
    #print(test_x)
    
    print('2.Training...')
    pass #nothing need to do in KNN in Training
    
    print('3.Testing...')
    numTest = test_x.shape[0]
    
    matched = 0
    for i in range(numTest):
        predict = classify(test_x[i], sample_x, sample_y, 3)
        if predict == test_y[i]:
            matched += 1
    accuracy = float(matched) / numTest
    print('4.the accuracy is %.2f%%' % (accuracy*100))
    
if __name__ == '__main__':
    cur = time.time()
    testHandWritingNum()
    print(str(time.time()-cur)+'s')
