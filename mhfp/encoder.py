import struct
import numpy as np

from hashlib import sha1
from rdkit.Chem import AllChem

class MHFPEncoder:
  """A class for encoding SMILES and RDKit molecule instances as MHFP fingerprints.
  """
  prime = (1 << 61) - 1
  max_hash = (1 << 32) - 1

  def __init__(self, n_permutations = 2048, seed = 42):
    """All minhashes created using this instance will use the arguments
    supplied to the constructor.

    Keyword arguments:
        n_permutations {int} -- The number of permutations used for hashing (default: {128})
        seed {int} -- The value used to seed numpy.random (default: {42})
    """  
    self.n_permutations = n_permutations
    self.seed = seed

    self.permutations_a = np.zeros([n_permutations], dtype=np.uint32)
    self.permutations_b = np.zeros([n_permutations], dtype=np.uint32)

    self.permutations_a_64 = np.zeros([n_permutations], dtype=np.uint64)
    self.permutations_b_64 = np.zeros([n_permutations], dtype=np.uint64)
    
    rand = np.random.RandomState(self.seed)

    # This is done in a loop as there shouldn't be any duplicate random numbers.
    # Also, numpy.random.choice seems to be implemented badly, as it throws
    # a memory error when supplied with a large n.
    for i in range(n_permutations):
      a = rand.randint(1, MHFPEncoder.max_hash, dtype=np.uint32)
      b = rand.randint(0, MHFPEncoder.max_hash, dtype=np.uint32)

      while a in self.permutations_a:
        a = rand.randint(1, MHFPEncoder.max_hash, dtype=np.uint32)
      
      while b in self.permutations_b:
        b = rand.randint(0, MHFPEncoder.max_hash, dtype=np.uint32)
      
      self.permutations_a[i] = a
      self.permutations_b[i] = b

    # Reshape into column vectors
    self.permutations_a = self.permutations_a.reshape((self.n_permutations, 1))
    self.permutations_b = self.permutations_b.reshape((self.n_permutations, 1))

  def encode(self, in_smiles, radius = 3, rings = True, kekulize = True, sanitize = True):
    """Creates an MHFP array from a SMILES string.
    
    Arguments:
      in_smiles {string} -- A valid SMILES string
      radius {int} -- The MHFP radius (a radius of 3 corresponds to MHFP6)  (default: {3})
      rings {boolean} -- Whether or not to include rings in the shingling (default: {True})
      kekulize {boolean} -- Whether or not to kekulize the extracted SMILES (default: {True})
      sanitize {boolean} -- Whether or not to sanitize the SMILES when parsing it using RDKit  (default: {True})
    
    Returns:
      numpy.ndarray -- An array containing the MHFP hash values.
    """

    return self.from_molecular_shingling(self.shingling_from_smiles(in_smiles, radius=radius, rings=rings, kekulize=kekulize, sanitize=sanitize))

  def encode_mol(self, in_mol, radius = 3, rings = True, kekulize = True):
    """ Creates an MHFP array from an RDKit molecule.
    
    Arguments:
      in_mol {rdkit.Chem.rdchem.Mol} -- A RDKit molecule instance
      radius {int} -- The MHFP radius (a radius of 3 corresponds to MHFP6)  (default: {3})
      rings {boolean} -- Whether or not to include rings in the shingling (default: {True})
      kekulize {boolean} -- Whether or not to kekulize the the extracted SMILES (default: {True})
    
    Returns:
      numpy.ndarray -- An array containing the MHFP hash values.
    """

    return self.from_molecular_shingling(self.shingling_from_mol(in_mol, radius=radius, rings=rings, kekulize=kekulize))

  def from_molecular_shingling(self, tokens):
    """Creates the hash set for a string array and returns it without changing the hash values of
    this instance.
    
    Arguments:
      a {numpy.ndarray} -- A string array.
    
    Returns:
      numpy.ndarray -- An array containing the hash values.
    """

    hash_values = np.zeros([self.n_permutations, 1], dtype=np.uint32)
    hash_values.fill(MHFPEncoder.max_hash)

    for t in tokens:
      t_h = struct.unpack('<I', sha1(t).digest()[:4])[0]
      hashes = np.remainder(np.remainder(self.permutations_a * t_h + self.permutations_b, MHFPEncoder.prime), self.max_hash)
      hash_values = np.minimum(hash_values, hashes)

    return hash_values.reshape((1, self.n_permutations))[0]

  def from_sparse_array(self, array):
    """Creates the hash set for a sparse binary array and returns it without changing the hash
    values of this instance. This is useful when a molecular shingling is already hashed.
    
    Arguments:
      s {numpy.ndarray} -- A sparse binary array.
    
    Returns:
      numpy.ndarray -- An array containing the hash values.
    """

    hash_values = np.zeros([self.n_permutations, 1], dtype=np.uint32)
    hash_values.fill(MHFPEncoder.max_hash)

    for i in array:
      hashes = np.remainder(np.remainder(self.permutations_a * i + self.permutations_b, MHFPEncoder.prime), self.max_hash)
      hash_values = np.minimum(hash_values, hashes)

    return hash_values.reshape((1, self.n_permutations))[0]

  @staticmethod
  def merge(a, b):
    """Merges (union) the two MHFP vectors.
    
    Arguments:
      a {numpy.ndarray} -- An array containing hash values.
      b {numpy.ndarray} -- An array containing hash values.
    Returns:
      numpy.ndarray -- An array containing the merged hash values.
    """
    return np.minimum(a, b)

  @staticmethod
  def merge_all(hash_values):
    """Merges (union) a list of hash_values.
    
    Arguments:
      hash_values {numpy.ndarray} -- An array of lists or arrays containing hash values.
    Returns:
      numpy.ndarray -- An array containing the merged hash values.
    """
    return np.amin(hash_values)

  @staticmethod
  def distance(a, b):
    """Estimates the Jaccard distance of two binary arrays based on their hashes.
    
    Arguments:
      a {numpy.ndarray} -- An array containing hash values.
      b {numpy.ndarray} -- An array containing hash values.
    
    Returns:
      float -- The estimated Jaccard distance.
    """
    return 1.0 - np.sum(a == b) / len(a)


  @staticmethod
  def shingling_from_mol(in_mol, radius=3, rings=True, kekulize=True):
    """Creates a molecular shingling from a RDKit molecule (rdkit.Chem.rdchem.Mol).
    
    Arguments:
      in_mol {rdkit.Chem.rdchem.Mol} -- A RDKit molecule instance
      radius {int} -- The MHFP radius (a radius of 3 corresponds to MHFP6)  (default: {3})
      rings {boolean} -- Whether or not to include rings in the shingling (default: {True})
      kekulize {boolean} -- Whether or not to kekulize the the extracted SMILES (default: {True})
    
    Returns:
      list -- The molecular shingling.
    """

    shingling = []
    
    if rings:
      for ring in AllChem.GetSymmSSSR(in_mol):
        bonds = set()
        ring = list(ring)
        for i in ring:
          for j in ring:
            if i != j:
              bond = in_mol.GetBondBetweenAtoms(i, j)
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

  @staticmethod
  def shingling_from_smiles(in_smiles, radius=3, rings=True, kekulize=True, sanitize=False):
    """Creates a molecular shingling from a SMILES string.
    
    Arguments:
      in_smiles {string} -- A valid SMILES string
      radius {int} -- The MHFP radius (a radius of 3 corresponds to MHFP6)  (default: {3})
      rings {boolean} -- Whether or not to include rings in the shingling (default: {True})
      kekulize {boolean} -- Whether or not to kekulize the extracted SMILES (default: {True})
      sanitize {boolean} -- Whether or not to sanitize the SMILES when parsing it using RDKit  (default: {False})
    
    Returns:
      list -- The molecular shingling.
    """

    return MHFPEncoder.shingling_from_mol(AllChem.MolFromSmiles(in_smiles, sanitize=sanitize), rings=rings, radius=radius, kekulize=True)