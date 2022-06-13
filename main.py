from ast import arg
from cmath import sqrt
from zmq import NULL
import time
from gmpy2 import iroot, mpz, divm, powmod, gcd, t_div, div
from gmpy2 import mpz_random, t_mod
import gmpy2
from multiprocessing import Process, Queue

from decoders import FactorizationDecoder, GcdDecoder, HastadDecoder, LinearPaddingHastadDecoder, SameModDecoder

def read_file(path) :
    with open(path) as f:
        seq = f.read()
    assert(len(seq) == 768)
    N = int(seq[:256], 16)
    e = int(seq[256:512], 16)
    me = int(seq[512:768], 16)
    return N, e, me

if __name__ == "__main__" :
    targets = []
    for i in range(21) :
        N, e, me = read_file("targets/data"+str(i))
        print(i, ":", N, e)
        targets.append((N, e, me, i))
    
    factorization_decoder = FactorizationDecoder([
        9686924917554805418937638872796017160525664579857640590160320300805115443578184985934338583303180178582009591634321755204008394655858254980766008932978699,
        10954856299233465126359914171500305822846165431085183673999109759449706417636519711881707731622506407722143163847672064459333431572992021257881551867597529
    ])
    hastad_decoder = HastadDecoder()
    linear_padding_decoder = LinearPaddingHastadDecoder()
    same_mod_decoder = SameModDecoder()
    gcd_decoder = GcdDecoder()

    factorization_decoder.decode(targets[1])
    factorization_decoder.decode(targets[3])
    factorization_decoder.decode(targets[8])
    factorization_decoder.decode(targets[10])
    factorization_decoder.decode(targets[12])

    same_mod_decoder.decode((targets[11], targets[14]))

    gcd_decoder.decode((targets[2], targets[19]))

    # 98764321261
    # 98764321277
    # 98764321291
    # 98764321301
    # 98764321307
    # 98764321361
    # padding: 123124123124123124123124123124123123
    def encode(m, pad, e, n) :
        return pow(m+pad, e, n)
    n1 = 98764321261*98764321277
    n2 = 98764321291*98764321301
    n3 = 98764321307*98764321361
    e = 3
    m = 142
    p1 = 1000000000
    p2 = 2222222222
    p3 = 9876543210

    # linear_padding_decoder.decode([
    #     (n1, e, encode(m, p1, e, n1), 111),
    #     (n2, e, encode(m, p2, e, n2), 111),
    #     (n3, e, encode(m, p3, e, n3), 111)
    # ], [p1, p2, p3])

    def get_pad(i) :
        pad_a = '0x9876543210abcdef'
        tmp = list('00000000')
        id = hex(i)[2:]
        for t in range(len(id)) :
            tmp[-len(id)+t] = id[t]
        pad_b = ''.join(tmp)
        pad_c = '0'*(128-16-8)
        pad = pad_a+pad_b+pad_c
        return int(pad, 16)

    # for id in [0, 6, 17] :
    #     N, e, me, i = targets[id]
    #     targets[i] = N, e, encode(0xffffffffffffffff, get_pad(id), e, N), i
    # linear_padding_decoder.decode([targets[0], targets[6], targets[17]], [get_pad(0), get_pad(6), get_pad(17)])
    for i in range(21) :
        linear_padding_decoder.decode([targets[0]], [get_pad(i)])
        linear_padding_decoder.decode([targets[6]], [get_pad(i)])
        linear_padding_decoder.decode([targets[17]], [get_pad(i)])
        # 0x9876543210abcdef0000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        # 0x9876543210abcdef0000000X00000000000000000000000000000000000000000000000000000000000000000000000000000000000000002054686174206973

    hastad_decoder.decode([targets[4], targets[5], targets[7], targets[13], targets[18]])

    