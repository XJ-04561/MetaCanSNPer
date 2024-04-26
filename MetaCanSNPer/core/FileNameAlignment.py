
class Table:
    def __init__(self, height, width):
        self._list = [[0 for _ in range(width)] for _ in range(height)]
    
    def __getitem__(self, key):
        return self._list[key[0]][key[1]]
    
    def __setitem__(self, key, value):
        self._list[key[0]][key[1]] = value

    def __str__(self):
        ret = []
        for l in self._list:
            ret.append(" ".join(["{:^5.1f}".format(e) for e in l]))
        return "\n".join(ret)

def checkDiag(table, y, x, value=0):
    n=0
    while y>0 and x>0:
        if table[y-1,x-1] == value:
            n+=1
        else:
            break
        y-=1
        x-=1
    return n

def align2(string1, string2, gapOpenCost=2, gapExtensionCost=1, matchReward=3, missMatchCost=2, contigBonus=0.1, gapSymbol="-", missMatchSymbol="X", trim=False, compressed=False):
    size = max(len(string1), len(string2))
    size1, size2 = len(string1), len(string2)
    scores = Table(size1+1, size2+1)
    outComes = Table(size1+1, size2+1)

    # Initiate starting conditions
    for i in range(size1+1):
        scores[i,0], outComes[i,0] = 0, -1
    for i in range(size2+1):
        scores[0,i], outComes[0,i] = 0, -1

    # Fill in Scores
    for y in range(size1):
        for x in range(size2):
            yt, xt = y+1, x+1
            gapString2  = scores[yt-1,xt]-(gapOpenCost if outComes[yt-1,xt] == 0 else gapExtensionCost)
            gapString1  = scores[yt,xt-1]-(gapOpenCost if outComes[yt,xt-1] == 0 else gapExtensionCost)
            noGap       = scores[yt-1,xt-1] + (matchReward if string1[y:y+1] == string2[x:x+1] else -missMatchCost)+contigBonus*checkDiag(outComes, yt, xt)
            scores[yt,xt], outComes[yt,xt] = max((gapString2, 1), (gapString1, 2), (noGap, 0), key=lambda x: x[0])
    
    m = (size1, size2)
    for y, x in list(zip(range(size1+1), [size2]*(size1+1))) + list(zip([size1]*(size2+1), range(size2+1))):
        if scores[y,x] > scores[m]:
            m = (y,x)
    
    y, x = m
    path = []
    while y > 0 and x > 0:
        if outComes[y,x] == 0:
            path.append(string1[y-1] if string1[y-1]==string2[x-1] else missMatchSymbol)
            y-=1
            x-=1
        elif outComes[y,x] == 1:
            path.append(gapSymbol)
            y-=1
        elif outComes[y,x] == 2:
            path.append(gapSymbol)
            x-=1
        else:
            raise ValueError("Not a recognized outCome symbol for align2()")

    if compressed:
        return gapSymbol.join([c for c in "".join(path[::-1]).strip(gapSymbol).split(gapSymbol) if c != ""])
    elif trim:
        return "".join(path[::-1]).strip(gapSymbol)
    else:
        return "".join(path[::-1])

def align(*strings, **kwargs):
    """def align(*strings, gapOpenCost=2, gapExtensionCost=1, matchReward=20, missMatchCost=1, gapSymbol="-", missMatchSymbol="#", trim=False)"""
    if len(strings) == 0:
        return -1
    elif len(strings) == 1:
        return strings[0]
    
    if len(strings) % 2:
        # Odd number of strings
        strings = (*strings, strings[0])
    
    nextStrings = []
    while len(strings) > 1:
        for i in range(len(strings)//2):
            nextStrings.append(align2(strings[2*i], strings[2*i+1], **kwargs))
        strings = nextStrings
        nextStrings = []
    return strings[0]


if __name__ == "__main__":
    print("Testing alignment scenarios.")
    aligned1 = align("FSC448_R1.fq", "FSC448_R2.fq")

    print("'{}' & '{}' = '{}'".format("FSC448_R1.fq", "FSC448_R2.fq", aligned1))

    aligned2 = align("FSC448.fq", "FSC448(1).fq")

    print("'{}' & '{}' = '{}'".format("FSC448.fq", "FSC448(1).fq", aligned2))
    
    aligned3 = align("FSC448.fq", "FSC658-[FSC448_R1].fq", "FSC567-[FSC448_R2].fq", compressed=True)

    print("'{}' & '{}' & '{}' = '{}'".format("FSC448.fq", "FSC658-[FSC448_R1].fq", "FSC567-[FSC448_R2].fq", aligned3))