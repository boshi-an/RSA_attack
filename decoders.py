from cmath import sqrt
from glob import escape
import os
from sqlite3 import DataError
from zmq import NULL
import time
from gmpy2 import iroot, mpz, divm, powmod, gcd, t_div, div, t_mod
from gmpy2 import mpz_random, t_mod
import gmpy2
import primefac # a efficiently implemented library for primalilty test and decomposition
from random import randint as rdint
from sage.all_cmdline import *   # import sage library

def get_d(p, q, e) :

    phi_n = (p - 1) * (q - 1)
    return divm(1, e, phi_n)

def get_m(N, d, me) :

    return powmod(me, d, N)

def deal_ans(task_id, p, q, d, m) :

    with open(os.path.join("answer", str(task_id)+".txt"), "w") as f:
        f.write(str(p)+"\n")
        f.write(str(q)+"\n")
        f.write(str(d)+"\n")
        f.write(str(m)+"\n")
        s = mpz(m).digits(16)[-16:]
        t = ""
        for i in range(8):
            tmp = (chr(int(s[2*i:2*i+2], 16)))
            if tmp == ' ':
                t += '\ '
            else:
                t += tmp
        f.write("msg: "+t+"\n")

        t = ""
        p = mpz(m).digits(16)[:-16]
        f.write("pad: "+p+"\n")

def already_computed(task_id) :
    path = os.path.join("answer", str(task_id)+".txt")
    if os.path.exists(path) :
        with open(path, "r") as f:
            p = f.read()
            q = f.read()
            d = f.read()
            m = f.read()
            return p, q, d, m
    else :
        return None

class FactorizationDecoder() :
    '''
    This decoder tries to factorize N to obtain p and q
    '''

    def __init__(self, hints) :

        self.rs = gmpy2.random_state()  # initialize random number generator
        self.hints = hints              # magic numbers, trying to factorize N

    def decode(self, target):

        N, e, me, task_id = target

        last_ans = already_computed(task_id)
        if last_ans != None :
            print("task{} already solved!".format(task_id))
            return last_ans
        
        flag = False

        # try to factorize using magic numbers
        for x in self.hints :
            if N % x == 0:
                p = x
                q = N//x
                flag = True
        
        if not flag :
            # uses all possible decomposition methods for factorization
            p = primefac.multifactor(N)[0]
            q = N//p

        print("task{} NEWLY solved!".format(task_id))

        d = get_d(p, q, e)
        m = get_m(N, d, me)

        deal_ans(task_id, p, q, d, m)

        return p, q, d, m

class HastadDecoder():
    '''
    Hastad Broadcast Attack
    If the same message was transmitted by different N with the same e,
    one can recover the message with Chinese Remainder Theorem
    '''

    def __init__(self) -> None:
        pass

    def decode(self, targets) :
        
        flag = True
        for N, e, me, i in targets :
            last_ans = already_computed(i)
            if last_ans != None :
                print("task{} already solved!".format(i))
            else :
                flag = False
        if flag == True :
           return last_ans[-1]

        x = 0
        mod = 1

        for N, e, me, i in targets :
            x = x + divm(1, mod, N) * (me-x) * mod
            mod = mod * N
            x = t_mod(t_mod(x, mod) + mod, mod)
        
        m, exact = iroot(x, len(targets))

        if not exact :

            print("Failed!")
            return None

        for N, e, me, i in targets :
            print("task{} NEWLY solved!".format(i))
            deal_ans(i, None, None, None, m)
        
        return m

class LinearPaddingHastadDecoder() :
    '''
    Coppersimth theorem
    If the padding of the code is know, and the secrete is relativly small,
    one can restore the secrete by LLL algorithm in poly time.
    '''

    def __init__(self) -> None:
        pass

    def decode(self, targets, pad) :
        #e_exp = len(targets) # e need to be the same as length of targets
        nArr = []
        aArr = []
        bArr = []
        cArr = []
        l = len(targets)
        for (N, e, me, i), p in zip(targets, pad):
            cArr.append(Integer(me))
            nArr.append(Integer(N))
            aArr.append(Integer(1))
            bArr.append(Integer(p))
        TArray = [-Integer(1) ]*l
        for i in range(l):
            arrayToCRT = [Integer(0)]*l
            arrayToCRT[i] = 1
            TArray[i] = crt(arrayToCRT,nArr)
        P = PolynomialRing(Zmod(prod(nArr)), names=('x',)); (x,) = P._first_ngens(1)
        gArray = [-Integer(1) ]*l
        for i in range(l):
            gArray[i] = TArray[i]*(pow(aArr[i]*x + bArr[i],e) - cArr[i])
            # print(pow(aArr[i]*x + bArr[i],e) - cArr[i])
        g = sum(gArray)
        # g = g.monic()
        # Use Sage's inbuilt coppersmith method
        roots = g.small_roots(epsilon=1/6)
        if(len(roots) == Integer(0) ):
            print("No Solutions found")
            return -Integer(1)
        else :
            for (N, e, me, i), p in zip(targets, pad):
                last_ans = already_computed(i)
                if last_ans != None :
                    print("task{} already solved!".format(i))
                else :
                    deal_ans(i, None, None, None, p+roots[0])
                    print("task{} NEWLY solved!".format(i))

        return roots[0]

class SameModDecoder() :
    '''
    If a secrete x is transmitted twice, by <N,e1> and <N,e2>,
    while gcd(e1, e2)=1, one can efficiently restore x, using
    extented Euclidean algorithm. 
    '''

    def __init__(self) -> None:
        pass
    
    def exgcd(self, x, y) :

        if x == 0:
            return 0, 1
        else :
            a_,b_ = self.exgcd(y%x, x)
            return b_-a_*(y//x), a_

    def decode(self, target) :

        N, e1, me1, i1 = target[0]
        N, e2, me2, i2 = target[1]

        s1, s2 = self.exgcd(e1, e2)

        m = (pow(me1, s1, N) * pow(me2, s2, N)) % N

        for N, e, me, i in target :
            last_ans = already_computed(i)
            if last_ans != None :
                print("task{} already solved!".format(i))
            else :
                deal_ans(i, None, None, None, m)
                print("task{} NEWLY solved!".format(i))

        return m

class GcdDecoder() :
    '''
    If fortunately, two unkonw numbers has a common divisor,
    we can go and by lottery.
    '''
    def __init__(self) -> None:
        pass

    def decode(self, target) :

        N1, e1, me1, i1 = target[0]
        N2, e2, me2, i2 = target[1]

        g = gcd(N1, N2)

        if g == 1 :
            return
        
        for N, e, me, i in target :
            last_ans = already_computed(i)
            if last_ans != None :
                print("task{} already solved!".format(i))
            else :
                d = get_d(g, N//g, e)
                m = get_m(N, d, me)
                deal_ans(i, g, N//g, d, m)
                print("task{} NEWLY solved!".format(i))

        return 

if __name__ == "__main__" :

    samemod = SameModDecoder()