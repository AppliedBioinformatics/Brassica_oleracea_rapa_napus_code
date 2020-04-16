
with open('bla') as fh:
    for line in fh:
        ll = line.split()
        # print(this_pair, shared_zero, shared_one, different, a_zeros, b_zeros)
        # evm.model.chrA03.260 evm.model.chrA09.7893 0 70 2 1 1
        first, second = ll[0], ll[1]
        shared_zero, shared_one, different, a_zeros, b_zeros = map(int, ll[2:])
        if a_zeros > 20 and b_zeros > 20 and different > 20:
            print(ll)
        
