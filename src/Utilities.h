#ifndef	__UTILITIES_H__
#define __UTILITIES_H__ 

#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <forward_list>
#include <chrono>
#include <cstdint>
#include <memory>
#include <unordered_set>
#include <map>
#include <unordered_map>
//#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <boost/bimap/bimap.hpp> //Bidirectional map

using namespace std;

using namespace std;
using namespace std::chrono;

//Type definition for all project
typedef uint16_t size_seq; //max 2^16=65.536
typedef uint64_t size_seq_tot; //max 2^64=1,844674407×10¹⁹

typedef uint32_t id_seq_type;
typedef uint16_t id_specie_type;
typedef uint32_t id_grp_type;
typedef uint16_t id_cluster_type;

//typedef uint128_t hash_type; //max 2^128=4^64=3,402823669×10³⁸
typedef uint64_t hash_type; //max 2^64=1,844674407×10¹⁹
typedef pair<hash_type, bool> HashCorrect;
typedef uint16_t lMer_type;

//For sequences
enum SeedState {No_State_Seed, No_Seed, Seed, No_Other_Read};
typedef unordered_map<id_seq_type, size_seq> MapAdjSeqID;
enum TypeGraph {Paired = 0, Single, SingleUnion};

typedef pair<string, string> SequenceHeader;
typedef map<id_seq_type, SequenceHeader> MapIDFile_Header;

//////////////////////////////////////////////////////////////////////////////////////////////////////
//Necessari per compilazione codice Assessment.h e Assessment.cpp
typedef map<id_specie_type, id_seq_type> Map_Specie__NumRead;

typedef map<id_grp_type, id_seq_type> Map_Grp__Size;
typedef map<id_grp_type, Map_Specie__NumRead> Map_Grp__Map_Specie__IDSeq;
typedef map<id_grp_type, Map_Specie__NumRead::value_type> Map_Grp__Max_Pair_Specie__IDSeq;

typedef map<id_cluster_type, id_seq_type> Map_Cluster__Size;
typedef map<id_cluster_type, Map_Specie__NumRead> Map_Cluster__Map_Specie__IDSeq;
typedef map<id_cluster_type, Map_Specie__NumRead::value_type> Map_Cluster__Max_Pair_Specie__IDSeq;

typedef map<id_specie_type, unordered_set<id_seq_type>> Map_Specie__SetIdSeq;
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct Lmer
{
	typedef shared_ptr<Lmer> Ptr;
	string l_mer;
	size_seq count = 1; //1 = non ha in sé il complemento inverso, 2 = ha il complemento inverso
	double prob = 0;
};

struct LMerVectorCompress
{
	typedef shared_ptr<LMerVectorCompress> Ptr;
	typedef map<lMer_type, lMer_type> Map_HashLMer_IndexVector;
	LMerVectorCompress(size_seq L);
	const Lmer::Ptr& GetWithIndex(lMer_type index);
	const Lmer::Ptr& GetWithHash(lMer_type hash);
	lMer_type GetIndexWithHash(lMer_type hash);
	size_seq getL() const;
	const Map_HashLMer_IndexVector& getMapHash() const;
	const vector<Lmer::Ptr>& getLmer() const;

private:
	size_seq L;
	Map_HashLMer_IndexVector mapHash;
	vector<Lmer::Ptr> vLmer;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////

class Spaced_Qmer {
public:
	//For unit
	struct Position_NOneBefore {
		size_t pos_start = 0;
		size_t n_one_before = 0;//numero di 1 che precedono il primo 1 dello unit
	};
	struct NOne_Position_NOneBefore {
		size_t n_one = 0;
		vector<Position_NOneBefore> v_pos_bef;
	};
	typedef vector<NOne_Position_NOneBefore> V_Unit;

	//for shift
	typedef vector<size_t> Position;
	struct PreviusShiftMinChange {
		vector<Position> one_to_change; //pos_to_change[i] -> i = index dello shift, il contenuto indica l'indice rispetto al vettore pos_one di quali uno cambiare
		Position one_exit; //numero uno da shiftare
		Position pos_min; //shift precedente migliore da cui computare hash
	};

	Spaced_Qmer();
	Spaced_Qmer(string spaced_qmer);
	inline size_t GetWeight() const {
		return this->pos_one.size();
	}
	inline size_t GetQ() const {
		return this->spaced_q.length();
	}
	inline bool isOne(size_t index) const {
		return this->spaced_q[index] == '1';
	}
	inline const Position& GetPosOne() const {
		return this->pos_one;
	}
	inline const V_Unit& GetUnitV() const {
		return this->unit_map;
	}
	inline const PreviusShiftMinChange& GetShiftMinChange() const {
		return this->shift_min_change;
	}
	inline const string& toString() const {
		return this->spaced_q;
	}
	void reset(string spaced_qmer);
private:
	string spaced_q;
	Position pos_one;
	V_Unit unit_map;
	PreviusShiftMinChange shift_min_change;

	void SaveIndexOne();
	void GetUnitV(V_Unit& v_unit);
	void GetShiftMax(PreviusShiftMinChange& shift_max);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////

inline static hash_type CharToInt(char ch)
{
	if(ch == 'A')
		return 0;
	if(ch == 'C')
		return 1;
	if(ch == 'G')
		return 2;
	if(ch == 'T')
		return 3;
	return 4; //ERROR CODE
}

inline static hash_type CharToIntComplement(char ch)
{
	if(ch == 'A')
		return 3;
	if(ch == 'C')
		return 2;
	if(ch == 'G')
		return 1;
	if(ch == 'T')
		return 0;
	return 4; //ERROR CODE
}

//Hash per tutti 1 su spaced qmer
inline static void GetHash(const string& s_Str, size_t startQmer, size_t length,
					HashCorrect& hash, hash_type (*fConvertion)(char))
{
	hash.first = 0;
	hash.second = true;
//	#pragma omp parallel for ordered
	for(size_t i = startQmer; i < startQmer + length; ++i)
	{
		hash_type ch = (*fConvertion)(s_Str[i]);
//		#pragma omp ordered
		if(hash.second)
		{
			if(ch == 4) //Errore conversione
				hash.second = false;
			if(hash.second)
				hash.first |= ch << ((i - startQmer) * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
		}
	}
	if(!hash.second)
		hash.first =  0;
}

//Hash per spaced qmer con *
inline static void GetHash(const string& s_Str, size_t startQmer, const Spaced_Qmer& spaced_qmer,
					HashCorrect& hash, hash_type (*fConvertion)(char))
{
	hash.first = 0;
	hash.second = true;
	const Spaced_Qmer::Position& pos_one = spaced_qmer.GetPosOne();
	for(size_t j = 0; j < pos_one.size(); ++j)
	{
		hash_type ch = (*fConvertion)(s_Str[startQmer+pos_one[j]]);
		if(hash.second)
		{
			if(ch == 4) //Errore conversione
				hash.second = false;
			if(hash.second)
				hash.first |= ch << (j * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
		}
	}
	if(!hash.second)
		hash.first =  0;
}

//Hash veloce con spaced qmer tutti 1
inline static void GetHashes_speedup_previous(const string& s_Str, size_t length,
						vector<HashCorrect>& vHash, hash_type (*fConvertion)(char))
{
	if(s_Str.size() >= length)
	{
		size_t n_hashes = s_Str.size() - length + 1;
		vHash.resize(n_hashes, HashCorrect(0, true)); //Crea vettore
		vector<size_t> err(n_hashes, 0); //Vettore che conta errori su qmer

		//Computa vettore che mi indica quanti errori ci sono nei qmer, quindi poi posso decidere se calcolare o no.
//		#pragma omp parallel for ordered
		for(size_t pos=0; pos < s_Str.size(); ++pos)
		{
			bool newErr = (*fConvertion)(s_Str[pos]) == 4;
//			#pragma omp ordered
			if(pos < length)
			{
				if(newErr)
					++err[0];
			}
			else
			{
				size_t actual = pos-length+1;
				size_t prev = pos-length;
				bool exitErr = (*fConvertion)(s_Str[prev]) == 4;
				err[actual] = err[prev];
				if(exitErr)
					--err[actual];
				if(newErr)
					++err[actual];
			}
		}
		//Imposto se computare o no l'hash
		#pragma omp parallel for
		for(size_t pos=0; pos < vHash.size(); ++pos)
		{
			if(err[pos] > 0)
			{
				vHash[pos].first = 0;
				vHash[pos].second = false;
			}
		}
		if(!vHash.empty())
		{
			GetHash(s_Str, 0, length, vHash[0], fConvertion);//primo da computare a parte
			for(size_t pos=1; pos < vHash.size(); ++pos)
			{
				if(vHash[pos].second) //Se devo computare
				{
					if(!vHash[pos-1].second)
						GetHash(s_Str, pos, length, vHash[pos], fConvertion);
					else
					{
						vHash[pos].first = vHash[pos - 1].first;
						vHash[pos].first -= (*fConvertion)(s_Str[pos - 1]); //sottrai primo elemento che viene eliminato
						vHash[pos].first >>= 2; //dividi per 4, sposta 2 bit
						vHash[pos].first |= ((*fConvertion)(s_Str[pos + length - 1]) << ((length - 1) * 2));	//aggiungi ultimo elemento OR possibile perchè prima ho
																												//diviso per 4 e la posizione dove scrivo ha sicuramente 0
					}
				}
			}
		}
	}
}

inline static void GetHashes_naive(const string& s_Str, const Spaced_Qmer& spaced_qmer,
						vector<HashCorrect>& vHash, hash_type (*fConvertion)(char))
{
	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
	if(isAllOne)
		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
	else
	{
		if(s_Str.size() >= spaced_qmer.GetQ())
		{
			size_t n_hashes = s_Str.size() - spaced_qmer.GetQ() + 1;
			vHash.resize(n_hashes, HashCorrect(0, true)); //Crea vettore
			#pragma omp parallel for
			for(size_t pos=0; pos < vHash.size(); ++pos)
				GetHash(s_Str, pos, spaced_qmer, vHash[pos], fConvertion);
		}
	}
}


inline static void GetHashes_speedup_unit(const string& s_Str, const Spaced_Qmer& spaced_qmer, vector<HashCorrect>& vHash, hash_type (*fConvertion)(char)) {
	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
	if(isAllOne)
		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
	else
	{
		if(s_Str.size() >= spaced_qmer.GetQ())
		{
			const Spaced_Qmer::V_Unit& spaced_v = spaced_qmer.GetUnitV();

			//Get hash v for all unit present
			//TODO: si può migliorare, dato che ordinate si può prendere quelle più piccole per comporre quelle più grandi
			vector<vector<HashCorrect>> hash_v(spaced_v.size());
			#pragma omp parallel for
			for(size_t i = 0; i < spaced_v.size(); ++i)//parallel computation
				GetHashes_speedup_previous(s_Str, spaced_v[i].n_one, hash_v[i], fConvertion);

			//Combine different hash
			size_t n_hashes = s_Str.size() - spaced_qmer.GetQ() + 1;
			vHash.resize(n_hashes, HashCorrect(0, true)); //Crea vettore
			#pragma omp parallel for
			for(size_t i = 0; i < vHash.size(); ++i)
			{
				for(size_t j = 0; j < spaced_v.size(); ++j)
				{
					vector<HashCorrect>& hash_unit = hash_v[j];
					for(size_t h = 0; h < spaced_v[j].v_pos_bef.size(); ++h)
					{
						size_t pos_unit = spaced_v[j].v_pos_bef[h].pos_start;
						size_t shift =  spaced_v[j].v_pos_bef[h].n_one_before*2;
						HashCorrect& hash = hash_unit[i+pos_unit];
						if(vHash[i].second && hash.second)
							vHash[i].first |= (hash.first << shift);
						else
						{
							vHash[i].first = 0;
							vHash[i].second = false;
						}
					}
				}
			}
		}
	}
}

inline static void GetHashes_speedup_previous(const string& s_Str, const Spaced_Qmer& spaced_qmer, vector<HashCorrect>& vHash, hash_type (*fConvertion)(char)) {
	bool isAllOne = spaced_qmer.GetWeight() == spaced_qmer.GetQ();
	if(isAllOne)
		GetHashes_speedup_previous(s_Str, spaced_qmer.GetQ(), vHash, fConvertion);
	else
	{
		auto get_hash = [](const string& s_Str, const Spaced_Qmer& spaced_qmer, vector<HashCorrect>& vHash,
				vector<Spaced_Qmer::Position>& vErr,
				const Spaced_Qmer::PreviusShiftMinChange& shift, size_t index_prev,
				const Spaced_Qmer::Position& pos_one, size_t i, hash_type (*fConvertion)(char)){
			size_t shift_of_one = shift.one_exit[index_prev];

			const Spaced_Qmer::Position one_to_change = shift.one_to_change[index_prev];//Recupera uno da cambiare

			size_t pos_hash_get = i-index_prev;//la posizione dell'hash presa è la posizione attuale meno l'indice dello shift dove si fan meno cambiamenti
			if(vHash[pos_hash_get].second)//se quello da cui devo prendere è corretto computa
			{
				//copia hash
				vHash[i] = vHash[pos_hash_get]; //Copia hash
				vHash[i].first >>= 2*shift_of_one;//Shifta correttamente

				//Controlla se attualmente hash è corretto
				Spaced_Qmer::Position& err = vErr[pos_hash_get];
				for(size_t e = 0; e < err.size(); ++e)
					if(err[e]>=shift_of_one)//1 non esce
						vErr[i].push_back(err[e]-shift_of_one);
				vHash[i].second = vErr[i].empty();

				//reset one
				hash_type reset_one = 0;
				for(size_t j = 0; j < one_to_change.size(); ++j)
					reset_one |= (hash_type)3 << (one_to_change[j] * 2);
				vHash[i].first &= ~reset_one;

				//aggiorna rimanenti posizioni su hash
				for(size_t j = 0; j < one_to_change.size(); ++j)
				{
					size_t index_char = i+pos_one[one_to_change[j]];
					hash_type ch = (*fConvertion)(s_Str[index_char]);
					if(ch == 4) //Errore conversione
					{
						vHash[i].second = false;
						vErr[i].push_back(one_to_change[j]);
					}
					if(vHash[i].second)
						vHash[i].first |= ch << (one_to_change[j] * 2);//OR possibile perchè sommo potenze di 4, OR su posizioni diverse, non c'è riporto
				}
			}
			else //computa da zero
				GetHash(s_Str, i, spaced_qmer, vHash[i], fConvertion);
		};
		if(s_Str.size() >= spaced_qmer.GetQ())
		{
			const Spaced_Qmer::PreviusShiftMinChange& shift = spaced_qmer.GetShiftMinChange();
			const Spaced_Qmer::Position& pos_one = spaced_qmer.GetPosOne();

			//Compute hash
			size_t n_hashes = s_Str.size() - spaced_qmer.GetQ() + 1;
			vHash.resize(n_hashes, HashCorrect(0, true)); //Crea vettore
			vector<Spaced_Qmer::Position> vErr(n_hashes);
			if(!vHash.empty())
			{
				GetHash(s_Str, 0, spaced_qmer, vHash[0], fConvertion);//primo da computare a parte
				for(size_t i = 1; i < shift.pos_min.size() && i < vHash.size(); ++i)//Per tutte le posizioni che contemplano gli shift nel primo pezzo di sequenza
					get_hash(s_Str, spaced_qmer, vHash, vErr, shift, shift.pos_min[i], pos_one, i, fConvertion);
				for(size_t i = shift.pos_min.size(); i < vHash.size(); ++i)
					get_hash(s_Str, spaced_qmer, vHash, vErr, shift, shift.pos_min.back(), pos_one, i, fConvertion);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void GetKmer(hash_type index, size_seq K, string& Kmer);
vector<string> GetAllKmers(size_seq K);

void createDirAndSubDir(string path);

int parseLineForMemory(char* line);
int getVirtualMemoryUsed();
int getPeakVirtualMemoryUsed();
int getPhysicalMemoryUsed();

#endif
