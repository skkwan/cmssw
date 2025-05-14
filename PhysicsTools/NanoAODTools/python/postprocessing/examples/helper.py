import ROOT
from itertools import combinations, product

# helper function to find leading pair
def findLeadingPair(collection):
    """
    Find the object pair in the collection with opposite-sign and the highest scalar pT sum.
    """
    max_sum_pt = 0
    sum_pt = 0
    idx_best_pair = 0

    combsAsList = list(combinations(collection, 2))

    for i in range(len(combsAsList)):
        pair = combsAsList[i]
        # reject same-electric-charge pairs
        if ((pair[0].charge * pair[1].charge) > 0):
            continue 
        # calculate pT
        sum_pt = (pair[0].pt + pair[1].pt)
        if (sum_pt > max_sum_pt):
            max_sum_pt = sum_pt 
            idx_best_pair = i 
    return combsAsList[idx_best_pair]

def findLeadingPairFromTwoCollections(coll1, coll2):
    """
    Find the object pair combination from two collections (e.g. with different lepton flavour)
    with opposite signs and the highest scalar pT sum.
    """
    max_sum_pt = 0
    sum_pt = 0
    idx_best_pair = 0
    combsAsList = list(product(coll1, coll2))
    for i in range(len(combsAsList)):
        pair = combsAsList[i]
        # reject same-electric-charge pairs
        if ((pair[0].charge * pair[1].charge) > 0):
            continue 
        # calculate pT
        sum_pt = (pair[0].pt + pair[1].pt)
        if (sum_pt > max_sum_pt):
            max_sum_pt = sum_pt 
            idx_best_pair = i 
    return combsAsList[idx_best_pair]

def invariantMass(objTuple):
    """
    Calculate the total invariant mass of two objects in a tuple
    """
    total_p4 = ROOT.TLorentzVector()
    total_p4 += objTuple[0].p4()
    total_p4 += objTuple[1].p4()
    return total_p4.M()