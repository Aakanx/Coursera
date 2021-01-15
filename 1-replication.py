def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for i in range(n-k+1):      
            if Text[i:i+k] == Pattern:
                freq[Pattern] += 1
    return freq

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

def Reverse(Pattern):
    rev = ''
    for char in Pattern: 
        rev = char + rev # can also do... for char in range(len(Pattern)): rev = Pattern[char] + rev
    return rev

def Complement(Pattern):
    com = ''
    for char in Pattern:
        if char == 'A':
            com += 'T'
        elif char == 'T':
            com += 'A'
        elif char == 'G':
            com += 'C'
        elif char == 'C':
            com += 'G'
    return com

def ReverseComplement(Pattern):   
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

def PatternMatching(Pattern, Genome):
    positions = []
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2]) # don't need n//2 - 1 because b in [:b] is excluded
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol: # if want window of 4 when genome size is 8, adding 4 to 0 gives 5 indices, so must subtract 1
            array[i] = array[i]+1
    return array

def SkewArray(Genome):
    Skew = [0]
    for i in range(len(Genome)):
        if Genome[i] == 'G': # don't forget Genome in Genome[i], can't do i by itself
            Skew.append(Skew[i]+1)
        elif Genome[i] == 'C':
            Skew.append(Skew[i]-1)
        else:
            Skew.append(Skew[i])
    return Skew

def MinimumSkew(Genome):
    positions = []
    Skew = SkewArray(Genome)
    for i in range(len(Skew)): # need range(len())
        if Skew[i] == min(Skew):
            positions.append(i)
    return positions

def HammingDistance(p, q): # number of mismatches
    count = 0
    for i in range(len(p)): # can also use... for a, b in zip(p, q): if a != b:
        if p[i] != q[i]:
            count += 1
    return count

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count += 1
    return count
