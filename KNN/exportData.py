from numpy import *
import operator
import os

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

def export():
    sample_x, sample_y, test_x, test_y = loadData()
    numSample = sample_x.shape[0]
    print(numSample)
    print(sample_x.shape[1])
    
    with open("sample.txt","w") as f:
        for i in range(numSample):
            for j in range(32*32):
                f.write(str(int(sample_x[i][j])))
            f.write(sample_y[i]+'\n')
    numTest = test_x.shape[0]
    with open("test.txt","w") as f:
        for i in range(numTest):
            for j in range(32*32):
                f.write(str(int(test_x[i][j])))
            f.write(str(test_y[i])+'\n')
    
    
if __name__ == '__main__':
    export()
