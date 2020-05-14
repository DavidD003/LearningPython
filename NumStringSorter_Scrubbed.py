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
            val = (abs((val-(val % 10**-6))**-1))*(10**20)
            for i in reversed(range(1, 30)):
                ky = ky + str(int(val) // (10 ** i))
                val = val-(int(val) // (10 ** i))*10**i
        return ky
    else:
        return 'b' + val


ls = [5, 4, 8, 'a', 'z', -1, 'ball', -15, 12.4, '!str', 'str',
      'aardvark', 'blast', 'bless', 7.8, 7.5, -.05, -5, -10.5, -10]

# lssorted = sorted(ls)

lssorted = sorted(ls, key=NumStringSorter)

print(lssorted)
