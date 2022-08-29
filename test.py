from itertools import chain, combinations, permutations
def _powerset(iterable, allow_empty_set = True): 
    s = set(iterable)
    if allow_empty_set:
        return set(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))
    else: 
        return [set(subset) for subset in chain.from_iterable(combinations(s, r) for r in range(1,len(s)+1))]

""" 
for element in _powerset(set([1,2,3,4]), allow_empty_set=False):
    print(element)
    print(set([1]).intersection(element))

print(_powerset(set([1,2,3,4]), allow_empty_set=False))
print(set.intersection(*(_powerset([1,2,3,4], allow_empty_set=False))))


print(set([1])) """
clustering_system = [set([1,2,4]), set([1,3]),set([1]), set([3])]
for triple in permutations(clustering_system, 3):
    print(triple)


print(set([1]).__eq__(set([1])))


print(set() in {frozenset()})