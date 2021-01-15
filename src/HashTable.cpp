#include"HashTable.h"

//Initialize some optional variables and the hash table
bool clsHashTable::InitOptionData(Spaced_Qmer Q, hash_type expected_element)
{
	this->Q=Q;
	this->arrMyHash.reserve(expected_element);
	return true;
}

Spaced_Qmer clsHashTable::getQ() const {
	return this->Q;
}

//Extract l-mer from a string which belongs to sequence having id iSeqInd
void clsHashTable::ExtractFromString(id_seq_type iSeqInd, string& s_Str,
		RepeatBloomFilter::Ptr repeatQmer) {
	vector<bool> repeat;
	if(repeatQmer)
		repeatQmer->areRepeats(s_Str, repeat);

	vector<HashCorrect> vHash;
	GetHashes_naive(s_Str, this->Q, vHash, CharToInt);
//	GetHashes_speedup_unit(s_Str, this->Q, vHash, CharToInt);
//	GetHashes_speedup_previous(s_Str, this->Q, vHash, CharToInt);

	//Insert q-mer in hash and save ID read where is present
	#pragma omp parallel for ordered
	for(size_seq i = 0; i < vHash.size(); ++i)
	{
		if(vHash[i].second)
		{
			#pragma omp ordered
			if(repeatQmer && repeat[i])
			{
				MapHash::mapped_type& node = this->arrMyHash[vHash[i].first];
				node.push_back(iSeqInd);
			}
			else if(!repeatQmer)
			{
				MapHash::mapped_type& node = this->arrMyHash[vHash[i].first];
				node.push_back(iSeqInd);
			}
		}
	}
}

void GetMapSeqID(MapAdjSeqID& map, KmNode& kmnode) {
	map.reserve(kmnode.size());
	auto id_it = kmnode.begin();
	while(id_it != kmnode.end())
	{
		++map[*id_it];
		++id_it;
	}
}
