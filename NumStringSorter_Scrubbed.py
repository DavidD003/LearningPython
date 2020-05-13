def NumStringSorter(val):
    if type(val) != str:
        if val >= 0:
            val = val * (10 ** 6)
            ky = 'ab'
            for i in reversed(range(1, 15)):
                ky = ky + str(int(val) // (10 ** i))
                val = val-(int(val) // (10 ** i))*10**i
        else:
            ky = 'a'
            val = abs((val-(val % 10**-6))**-1)*(10**20)
            for i in reversed(range(1, 15)):
                ky = ky + str(int(val) // (10 ** i))
                val = val-(int(val) // (10 ** i))*10**i
        return ky
    else:
        return 'b' + val
