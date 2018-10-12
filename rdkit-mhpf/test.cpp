#include <iostream>
#include <algorithm>    // std::unique_copy, std::sort, std::distance
#include <random>
#include <GraphMol/MolOps.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include "sha1.hpp"

using namespace RDKit;
struct compare
{
    int key;
    compare(int const &i): key(i) { }

    bool operator()(int const &i)
    {
        return (i == key);
    }
};

struct distrib{
	std::vector<uint64_t> da;
	std::vector<uint64_t> db;
	long int max_hash;
	long long int prime;
};

// generator of distribution for Hashing
distrib generatordistrib(int seed, int distribsize){

	distrib dab;
	long int max_hash = std::llround(std::pow(2,32)) - 1; // fix
    long long int prime = std::llround(std::pow(2,61)) -1; // fix
 	std::vector<uint64_t> v1;
 	std::vector<uint64_t> v2;

	// generator 1:
    std::mt19937::result_type myseed = seed;
    std::mt19937 e1(myseed); //_64 to have 64 bits double not needed
    std::uniform_int_distribution<long long int> dist1(1, max_hash); // fix dist1
    while (v1.size()<distribsize) {
    	long long int a = dist1(e1);
        bool result = !std::any_of(v1.begin(),v1.end(),compare(a));
        if (result) {
        	v1.push_back(a);

    	}
    }

    // generator 2:
    std::mt19937 e2(myseed); //_64 to have 64 bits double not needed
    std::uniform_int_distribution<long long int> dist2(0, max_hash); // fix dist2
    while (v2.size()<distribsize) {
    	long long int b = dist2(e2);
        bool result = !std::any_of(v2.begin(),v2.end(),compare(b));
        if (result) {
        	v2.push_back(b);
    	}
    }

    dab.da=v1;
    dab.db=v2;
    dab.prime=prime;
    dab.max_hash=max_hash;
    return dab;


}


std::vector<std::string> shingling_from_mol(RDKit::ROMol *in_mol, int radius, bool rings, bool kekulize) {
    std::vector<std::string> shingling;
    if (rings){
    	PATH_TYPE bonds;
		RDKit::VECT_INT_VECT molrings;
  		RDKit::MolOps::symmetrizeSSSR(*in_mol, molrings);
  		for (const auto molring : molrings) {
  			for (const auto index1: molring) {
  				//std::cout << "index1 : "<< index1 << std::endl;
      			for (const auto index2: molring) {
      				//std::cout << "index2 : "<< index2 << std::endl;
      				if(index2!=index1) {
      					const Bond *bnd = in_mol->getBondBetweenAtoms(index1, index2);
      					if (bnd) {
      						bool result = !std::any_of(bonds.begin(),bonds.end(),compare(bnd->getIdx()));
      						if (result) {
      							bonds.push_back(bnd->getIdx());
      						}
                       	}
                    }
                }
    		}
  		}

  		// by default canonical MolToSmiles
  		RDKit::ROMol *submol = Subgraphs::pathToSubmol(*in_mol, bonds);
  		std::string subsmi = RDKit::MolToSmiles(*submol);
        shingling.push_back(subsmi);
	}



	PATH_TYPE p;
	std::string smiles;
    for (int index = 0; index < in_mol->getNumAtoms(); index++){
      for (int i = 1; i < radius + 1;  i++){
      	p = RDKit::findAtomEnvironmentOfRadiusN(*in_mol, i, index);
        INT_MAP_INT amap; // map
        ROMol *submol = Subgraphs::pathToSubmol(*in_mol, p, false, amap);
        // check if the atomindex match the current index value!
        if (amap.count(index) > 0){
	        smiles = RDKit::MolToSmiles(*submol, true, false, amap[index], true);
	        if (!smiles.empty()){
	           shingling.push_back(smiles);
	        }
        }
      }
	}
	// remove duplicates then!
	sort (shingling.begin(), shingling.end());
    shingling.erase(std::unique(shingling.begin(), shingling.end()), shingling.end());
    return shingling;
}


std::vector<uint32_t> from_molecular_shingling(std::vector<std::string> token, distrib dab, int distribsize){
	// set the hash_values
	std::vector<uint32_t> hash_values;

	for (int i=0;i<distribsize;i++){
		hash_values.push_back(dab.max_hash);
	}
	// generate hashes
	for(const auto res : token){
		SHA1 checksum;
    	checksum.update(res);
    	uint32_t t_h = checksum.Hash(); // not same as in python
        std::cout << "t_h : "<< t_h << std::endl;

    	for (int j=0;j<distribsize;j++){
			hash_values[j] = fmin( hash_values[j] , remainder( remainder( dab.da[j] * t_h + dab.db[j], dab.prime) , dab.max_hash) );
		}
	}
	return hash_values;
}



int main()
{
	int bitsize = 2048;
    distrib dab = generatordistrib(42, bitsize);
    //std::cout << dab.da.size()<< ", " << dab.db.size()<< std::endl ;
	std::cout << "prime : "<< dab.prime << std::endl;
	std::cout << "maxhash : "<< dab.max_hash << std::endl;

	RDKit::ROMol *mol1 = RDKit::SmilesToMol("c1ccccc1");

	std::vector<std::string> res = shingling_from_mol(mol1, 3, true, true);

    std::vector<uint32_t> hash_values = from_molecular_shingling(res, dab, bitsize);
	for (int i=0;i<bitsize;i++){
		std::cout <<  hash_values[i]  << ",";
	}

	std::cout << std::endl;

}

