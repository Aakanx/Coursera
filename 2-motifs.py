def Count(Motifs):  # takes list of strings, returns dictionary
    count = {}  # initializing the count matrix
    k = len(Motifs[0])  # length of the first kmer (but all are same length)
    for symbol in "ACGT":
        count[symbol] = []  # count matrix now has keys A, C, T, and G all with values of empty list
        for j in range(k):
            count[symbol].append(
                0)  # count matrix now has keys A, C, G, and T all with values of a list of zeroes of length equal to the length of a kmer
    t = len(Motifs)  # length of Motifs, a list of kmers (strings)
    for i in range(t):  # for each kmer in Motifs
        for j in range(k):  # for each element of the kmer
            symbol = Motifs[i][j]  # assigns a nucleotide (ACGT) in Motifs to the key (symbol)
            # count[symbol] corresponds to the key of the dictionary count
            # count[symbol][j] corresponds to the position in the list assigned to the key
            count[symbol][j] += 1  # adds 1 to the position in the list assigned to the key
    return count

def Profile(Motifs):  # takes list of strings, returns dictionary
    profile = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            profile[symbol][j] += 1 / t
    return profile


def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0  # resets after nested for loop when at next j
        frequentSymbol = ""  # resets after nested for loop when at next j
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(Motifs):
    score = 0
    consensus = Consensus(Motifs)  # string, just making variable name shorter here, not necessary
    t = len(Motifs)  # i
    k = len(Motifs[0])  # j
    for j in range(k):
        for i in range(t):
            if consensus[j] != Motifs[i][j]:
                score += 1  # can use score += HammingDistance(Motifs[i], Consensus(Motifs)) as well
    return score


def Entropy(Motifs):  # like score, but takes into account variation within column, treating like probability distribution
    from math import log
    profile = Profile(Motifs)
    entropy = 0
    k = len(Motifs[0])  # j, number of columns
    for symbol in "ACGT":
        for j in range(k):
            if profile[symbol][j] != 0:
                entropy += profile[symbol][j] * log(profile[symbol][j], 2)
    entropy *= -1
    return entropy


def Pr(Text, Profile):  # probability of text, given profile (given motifs)
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p


def ProfileMostProbableKmer(Text, k, Profile):  # most probable kmer in string Text given Profile
    probability = -1  # -1 to allow strings with probability = 0
    MostProbableKmer = ""
    for i in range(len(Text) - k + 1):
        if Pr(Text[i:i + k], Profile) > probability:  # each kmer in Text
            probability = Pr(Text[i:i + k], Profile)
            MostProbableKmer = Text[i:i + k]
    return MostProbableKmer


def GreedyMotifSearch(Dna, k, t):  # check score per string instead of per set of strings (which would be exhaustive)
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])  # sliding window through first string (0th string) in Dna
        for j in range(1, t):  # going through rest of strings in Dna, generating set of Motifs to be scored
            P = Profile(Motifs[0:j])  # until, excluding j
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))  # next string based on profile of previous string(s)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    countp = {}
    for symbol in "ACGT":
        countp[symbol] = []
        for j in range(k):
            countp[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            countp[symbol][j] += 1
    return countp


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    countp = CountWithPseudocounts(Motifs)
    profilep = {}
    for symbol in "ACGT":
        profilep[symbol] = []  # set up
        for j in range(k):
            profilep[symbol].append(0)  # set up
            profilep[symbol][j] = countp[symbol][j] / (t + 1 * 4)
    return profilep


def Consensus(Motifs):  # with pseudocounts!
    k = len(Motifs[0])
    countp = CountWithPseudocounts(Motifs)
    consensusp = ""
    for j in range(k):
        m = 0  # resets after nested for loop when at next j
        frequentSymbol = ""  # resets after nested for loop when at next j
        for symbol in "ACGT":
            if countp[symbol][j] > m:
                m = countp[symbol][j]
                frequentSymbol = symbol
        consensusp += frequentSymbol
    return consensusp


def Score(Motifs):  # with pseudocounts!
    scorep = 0
    consensusp = Consensus(Motifs)  # string, just making variable name shorter here, not necessary
    t = len(Motifs)  # i
    k = len(Motifs[0])  # j
    for j in range(k):
        for i in range(t):
            if consensusp[j] != Motifs[i][j]:
                scorep += 1  # can use score += HammingDistance(Motifs[i], Consensus(Motifs)) as well
    return scorep


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])  # sliding window through first string
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])  # until and excluding j
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))  # next string based on profile
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def Motifs(Profile, Dna):
    k = len(Profile["A"])
    t = len(Dna)
    Motifs = []
    for i in range(t):  # string by string in Dna
        Motifs.append(ProfileMostProbableKmer(Dna[i], k, Profile))
    return Motifs


def RandomMotifs(Dna, k, t):
    import random
    RandomMotifs = []
    for i in range(t):
        r = random.randint(0, len(Dna[0]) - k)
        RandomMotifs.append(Dna[i][r:r + k])
    return RandomMotifs


def RandomizedMotifSearch(Dna, k, t):  # works because Profile captures bias in Dna towards Consensus
    BestMotifs = RandomMotifs(Dna, k, t)
    while True:
        profilep = ProfileWithPseudocounts(BestMotifs)
        if Score(Motifs(profilep, Dna)) < Score(BestMotifs):
            BestMotifs = Motifs(profilep, Dna)
        else:
            return BestMotifs


def RandomizedMotifSearch(Dna, k, t, N):
    NBestMotifs = RandomMotifs(Dna, k, t)
    for i in range(N):
        BestMotifs = RandomMotifs(Dna, k, t)
        while True:  # until score of Motifs stops being better than score of BestMotifs
            profilep = ProfileWithPseudocounts(BestMotifs)
            if Score(Motifs(profilep, Dna)) < Score(BestMotifs):
                BestMotifs = Motifs(profilep, Dna)
            # if didn't generate Motifs of lower score, continue to elif
            elif Score(BestMotifs) < Score(NBestMotifs):
                NBestMotifs = BestMotifs
                break
            else:
                break
    return NBestMotifs


def Normalize(Probabilities):
    normalized = {}
    for kmer in Probabilities:
        normalized[kmer] = Probabilities[kmer] / sum(Probabilities.values())
    return normalized


def WeightedDie(Probabilities):
    import random
    p = random.uniform(0, 1)
    count = 0
    for k, v in Probabilities.items():
        if count <= p < count + v:
            return k
        count += v


# p = random.uniform(0, 1)
# for kmer in Probabilities:
#    p -= Probabilities[kmer]
#    if p <= 0:
#         return kmer


def ProfileGeneratedString(Text, Profile, k):
    n = len(Text)
    Probabilities = {}
    for i in range(0, n - k + 1):
        Probabilities[Text[i:i + k]] = Pr(Text[i:i + k], Profile)
    Probabilities = Normalize(Probabilities)
    return WeightedDie(Probabilities)


def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(N):
        i = random.randint(1, t - 1)
        LessMotifs = []
        for a in range(t):
            if a != i:
                LessMotifs.append(Motifs[a])
        Motifs[i] = ProfileGeneratedString(Dna[i], ProfileWithPseudocounts(LessMotifs), k)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
