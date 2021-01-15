#ifndef	__UTILITIES_CPP__
#define __UTILITIES_CPP__

#include"Utilities.h"

LMerVectorCompress::LMerVectorCompress(size_seq L) {
	this->L = L;
	vector<string> allLmers = GetAllKmers(L);
	for(lMer_type i = 0; i < allLmers.size(); i++)
	{
		Lmer::Ptr l_mer = make_shared<Lmer>();
		l_mer->l_mer = allLmers[i];
		this->vLmer.push_back(l_mer);
	}

	for(lMer_type i = 0; i < this->vLmer.size(); ++i)
	{
		//Sicuramente hashLmer assegnato a indice i
		HashCorrect hashLmer;
		GetHash(this->vLmer[i]->l_mer, 0, L, hashLmer, CharToInt);
		this->mapHash[hashLmer.first] = i;

		//Vediamo dove inserire il reverse, creo hash del reverse
		string revLm = "";
		reverse_copy(this->vLmer[i]->l_mer.begin(), this->vLmer[i]->l_mer.end(), revLm.begin());
		HashCorrect hashReverseLmer;
		GetHash(revLm, 0, L, hashReverseLmer, CharToIntComplement);

		for(lMer_type j = i + 1; j < this->vLmer.size(); ++j)
		{
			HashCorrect hashLmerNext;
			GetHash(this->vLmer[j]->l_mer, 0, L, hashLmerNext, CharToInt);
			if(hashLmerNext.first == hashReverseLmer.first)//if(revLm == this->vLmer[j]->l_mer) //Se il reverse è presente lo cancello e metto tutto su indice hashLmer
			{
				++this->vLmer[i]->count;
				this->vLmer.erase(this->vLmer.begin() + j);
				this->mapHash[hashReverseLmer.first] = i;
				break;
			}
		}
	}
}


const Lmer::Ptr& LMerVectorCompress::GetWithIndex(lMer_type index) {
	return this->vLmer[index];
}

const Lmer::Ptr& LMerVectorCompress::GetWithHash(lMer_type hash) {
	return this->vLmer[this->mapHash.at(hash)];
}

lMer_type LMerVectorCompress::GetIndexWithHash(lMer_type hash) {
	return this->mapHash.at(hash);
}

size_seq LMerVectorCompress::getL() const {
	return L;
}

const LMerVectorCompress::Map_HashLMer_IndexVector& LMerVectorCompress::getMapHash() const {
	return mapHash;
}

const vector<Lmer::Ptr>& LMerVectorCompress::getLmer() const {
	return vLmer;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

Spaced_Qmer::Spaced_Qmer()
{
	this->reset("");
}

Spaced_Qmer::Spaced_Qmer(string spaced_qmer) {
	this->reset(spaced_qmer);
}

void Spaced_Qmer::reset(string spaced_qmer) {
	this->spaced_q = spaced_qmer;
	this->SaveIndexOne();
	this->GetUnitV(this->unit_map);
	this->GetShiftMax(this->shift_min_change);
}

void Spaced_Qmer::SaveIndexOne() {
	this->pos_one.clear();this->pos_one.shrink_to_fit();
	for(size_t i = 0; i < this->spaced_q.length(); ++i)
		if(this->isOne(i))
			this->pos_one.push_back(i);
}

void Spaced_Qmer::GetUnitV(V_Unit& v_unit) {
	//reset
	v_unit.clear();v_unit.resize(this->pos_one.size()+1); //Dimensione massima è il weight+1
	//Calcola
	size_t init_unit = 0;
	size_t size_unit = 0;
	for(size_t i = 0; i < this->pos_one.size(); ++i)
	{
		if(this->pos_one[i] != this->pos_one[init_unit]+size_unit)
		{
			Position_NOneBefore new_pos_nonebef;
			v_unit[size_unit].n_one = size_unit; //inserisco su vettore, poi elimino quindi perdo posizione
			new_pos_nonebef.pos_start = this->pos_one[init_unit];
			new_pos_nonebef.n_one_before = init_unit;
			v_unit[size_unit].v_pos_bef.push_back(move(new_pos_nonebef));

			init_unit = i;
			size_unit = 0;
		}
		++size_unit;
	}
	if(size_unit > 0)
	{
		Position_NOneBefore new_pos_nonebef;
		v_unit[size_unit].n_one = size_unit;
		new_pos_nonebef.pos_start = this->pos_one[init_unit];
		new_pos_nonebef.n_one_before = init_unit;
		v_unit[size_unit].v_pos_bef.push_back(move(new_pos_nonebef));
	}
	//Remove zero
	auto it = v_unit.begin();
	while(it != v_unit.end())
	{
		if(it->n_one == 0)
			it = v_unit.erase(it);
		else
			++it;
	}
	v_unit.shrink_to_fit();
}

void Spaced_Qmer::GetShiftMax(PreviusShiftMinChange& shift_max) {
	//reset
	this->shift_min_change.pos_min.clear();
	this->shift_min_change.pos_min.resize(this->spaced_q.size(), 0);
	this->shift_min_change.pos_min.shrink_to_fit();

	this->shift_min_change.one_exit.clear();
	this->shift_min_change.one_exit.resize(this->spaced_q.size(), 0);
	this->shift_min_change.one_exit.shrink_to_fit();

	this->shift_min_change.one_to_change.clear();
	this->shift_min_change.one_to_change.resize(this->spaced_q.size());
	this->shift_min_change.one_to_change.shrink_to_fit();
	//calcola posizioni da salvare con i vari shift e poi salva
	size_t init = 0;
	for(size_t i = 1; i < this->spaced_q.size(); ++i)//per tutti gli shift possibili a parte il primo
	{
		//Cerco indice di partenza su vettore uno
		bool find = false;
		for(size_t j = init; j < this->pos_one.size(); ++j)
		{
			if(this->pos_one[j] >= i)
			{
				init = j;
				find = true;
				break;
			}
		}
		if(!find)
			init = this->pos_one.size();//Serve per saltare prossimo ciclo senza altri controlli
		for(size_t j = init; j < this->pos_one.size(); ++j) //per tutte le posizioni del secondo vettore
		{
			//confronta indici incolonnati tenendo conto dello shift, verifica se diversi, se si c'è un operazione da fare in quel punto
			if(this->pos_one[j-init] != this->pos_one[j]-i)
				this->shift_min_change.one_to_change[i].push_back(j-init);
		}
		for(size_t j = this->pos_one.size()-init; j < this->pos_one.size(); ++j)//rimanenti da inserire tutti
			this->shift_min_change.one_to_change[i].push_back(j);

		//Aggiorna minimo e pos_one exit, attenzione con i==1
		this->shift_min_change.one_exit[i] = init;
		size_t previous_shift_min = this->shift_min_change.pos_min[i-1];
		Position& previous_shift_change_min = this->shift_min_change.one_to_change[previous_shift_min];
		Position& current_shift_change = this->shift_min_change.one_to_change[i];
		size_t current_shift_min = (previous_shift_change_min.size() < current_shift_change.size()) ? previous_shift_min : i;
		this->shift_min_change.pos_min[i] = (i == 1) ? i : current_shift_min;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetKmer(hash_type index, size_seq K, string& Kmer) {
	if(index >= ((hash_type)1 << K*2))//(hash_type)boost::multiprecision::pow(boost::multiprecision::uint1024_t(4), K)) //Superato indice massimo possibile
		return;
	vector<string> nucleotide = {"A","C","G","T"};
	vector<int> indexNucleotide;
	for(size_seq i = 0; i < K; ++i)
	{
		indexNucleotide.push_back(index%4);
		index /= 4;
	}
	for(size_seq i = K; i > 0; --i)
		Kmer.append(nucleotide[indexNucleotide[i-1]]);
}

vector<string> GetAllKmers(size_seq K) {
	hash_type max = (hash_type)1 << K*2;//(hash_type)boost::multiprecision::pow(boost::multiprecision::uint1024_t(4), K);
	vector<string> allKmer;
	for(hash_type index = 0; index < max; ++index)
	{
		string LMer = "";
		GetKmer(index, K, LMer);
		allKmer.push_back(LMer);
	}
	return allKmer;
}

void createDirAndSubDir(string path) {
	string dir = "";
	string delimiter = "/";
	size_t pos = path.find(delimiter);
	while(pos != path.npos)
	{
		dir += path.substr(0, pos + delimiter.length());
		path.erase(0, pos + delimiter.length());
		mkdir(dir.c_str(), S_IRWXU);
		pos = path.find(delimiter);
	}
}

int parseLineForMemory(char* line){
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}


int getVirtualMemoryUsed(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPeakVirtualMemoryUsed() { //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmPeak:", 7) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPhysicalMemoryUsed(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLineForMemory(line);
            break;
        }
    }
    fclose(file);
    return result;
}

#endif
