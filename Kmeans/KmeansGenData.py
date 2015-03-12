'''
Name        : KmeansGenData.py
Author      : wyj

Date        : 2015.03.12

Version     : 1.0
Description : Generate data to test Kmeans algorithm. Generate random point in dim dimention 
rectangular [0,bound]...[bound,0] and [-bound,0]...[0,-bound]. tot points in each rectangular.
'''
import random

if __name__ == '__main__':
    bound = 2
    tot = 50
    dim = 2
    
    with open("input.txt","w") as f:
        f.write(str(tot*2)+' '+str(dim)+'\n')
        for i in range(0,tot):
            for d in range(0,dim):
                f.write(str(random.uniform(0,bound))+' ')
            f.write('\n')
            for d in range(0,dim):
                f.write(str(random.uniform(-bound,0))+' ')
            f.write('\n')

