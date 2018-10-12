import numpy as np
from rdkit.Chem import AllChem
from hashlib import sha1
import struct


prime = (1 << 61) - 1
max_hash = (1 << 32) - 1

print(prime)
print(max_hash)

def gen(n_permutations = 2048, seed = 42):
    """All minhashes created using this instance will use the arguments
    supplied to the constructor.

    Keyword arguments:
        n_permutations {int} -- The number of permutations used for hashing (default: {128})
        seed {int} -- The value used to seed numpy.random (default: {42})
    """  

    permutations_a = np.zeros([n_permutations], dtype=np.uint32)
    permutations_b = np.zeros([n_permutations], dtype=np.uint32)

    #print(permutations_a)
    #print(permutations_b)

    permutations_a_64 = np.zeros([n_permutations], dtype=np.uint64)
    permutations_b_64 = np.zeros([n_permutations], dtype=np.uint64)

    rand = np.random.RandomState(seed)

    # This is done in a loop as there shouldn't be any duplicate random numbers.
    # Also, numpy.random.choice seems to be implemented badly, as it throws
    # a memory error when supplied with a large n.
    for i in range(n_permutations):
      a = rand.randint(1, max_hash, dtype=np.uint32)
      #print(a)
      b = rand.randint(0, max_hash, dtype=np.uint32)

      while a in permutations_a:
        a = rand.randint(1, max_hash, dtype=np.uint32)
        #print("whilea:")
        #print(a)
      
      while b in permutations_b:
        b = rand.randint(0, max_hash, dtype=np.uint32)
        #print("whileb:")
        #print(b)
      
      permutations_a[i] = a
      permutations_b[i] = b


    # Reshape into column vectors
    permutations_a = permutations_a.reshape((n_permutations, 1))
    permutations_b = permutations_b.reshape((n_permutations, 1))
    return permutations_a, permutations_b


def shingling_from_mol(in_mol, radius=3, rings=True, kekulize=True):
    shingling = []
    if rings:
      for ring in AllChem.GetSymmSSSR(in_mol):
        bonds = set()
        ring = list(ring)
        for i in ring:
          for j in ring:
            if i != j:
              bond = in_mol.GetBondBetweenAtoms(i, j)
              print(bond)
              if bond != None:
                bonds.add(bond.GetIdx())
        shingling.append(AllChem.MolToSmiles(AllChem.PathToSubmol(in_mol, list(bonds)), canonical=True).encode('utf-8'))

    for index, _ in enumerate(in_mol.GetAtoms()):
      for i in range(1, radius + 1):
        p = AllChem.FindAtomEnvironmentOfRadiusN(in_mol, i, index)
        amap = {}
        submol = AllChem.PathToSubmol(in_mol, p, atomMap=amap)
        if index not in amap:
          continue
        smiles = AllChem.MolToSmiles(submol, rootedAtAtom=amap[index], canonical=True)
        if smiles is not '':
          shingling.append(smiles.encode('utf-8'))

    # Set ensures that the same shingle is not hashed multiple times
    # (which would not change the hash, since there would be no new minima)
    return list(set(shingling))


mol = AllChem.MolFromSmiles("c1ccccc1");

res = shingling_from_mol(mol, 3);
permutations_a,permutations_b =gen()
hash_values = np.zeros([2048, 1], dtype=np.uint32)
hash_values.fill(max_hash)

for t in res:
    t_h = struct.unpack('<I', sha1(t).digest()[:4])[0]
    print(t_h)
    hashes = np.remainder(np.remainder(permutations_a * t_h + permutations_b, prime), max_hash)
    hash_values = np.minimum(hash_values, hashes)

print(hash_values.reshape((1, 2048))[0])