import random
import math

if __name__ == '__main__':
    n = 1000
    m = 2
    r1 = 1
    r2 = 2
    sqr1 = r1**2
    sqr2 = r2**2
    ch = [1,-1]
    
    with open(r"in.txt","w") as f:
        f.write(str(n*2)+' '+str(m)+'\n')
        for i in range(n):
            x = random.uniform(-r1,r1)
            tmp = math.sqrt(sqr1-x**2)
            y = random.uniform(-tmp,tmp)
            f.write(str(x)+' '+str(y)+'\n')
        for i in range(n):
            x = random.uniform(-r2,r2)
            tmp = math.sqrt(sqr2-x**2)
            if x > 1 or x < -1:
                y = random.uniform(0,tmp)
            else:
                yy = math.sqrt(sqr1-x**2)
                y = random.uniform(yy,tmp)
            signY = random.choice(ch)
            f.write(str(x)+' '+str(signY*y)+'\n')
            #f.write(str(x)+' '+str(y)+'\n')
            
    with open(r'out.txt','w') as f:
        f.write(str(n*2)+'\n')
        for i in range(n):
            f.write('1\n')
        for i in range(n):
            f.write('-1\n')
