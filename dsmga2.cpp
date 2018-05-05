/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <list>
#include <vector>
#include <algorithm>
#include <iterator>

#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"

#include <iomanip>
using namespace std;


DSMGA2::DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff) {


    previousFitnessMean = -INF;
    ell = n_ell;
    nCurrent = (n_nInitial/2)*2;  // has to be even
    nOrig = nCurrent;

    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    selectionPressure = 2;
    maxGen = n_maxGen;
    maxFe = n_maxFe;

    graph.init(ell);
    graph_size.init(ell);
    orig_graph.init(ell);
    orig_graph_size.init(ell);
     
    bestIndex = -1;
    masks = new list<int>[ell];
    selectionIndex.resize(nCurrent);
    orig_selectionIndex.resize(nCurrent);
    orderN.resize(nCurrent);
    orderELL = new int[ell];
    fastCounting = new FastCounting[ell];
    orig_fc = new FastCounting[ell];
    population.clear();
    usedIndex = new bool[ell];

    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);


    pHash.clear();
    pHashOrig.clear();

    for (int i=0; i<nCurrent; ++i) {
        Chromosome ch;
        population.push_back(ch);
        population[i].initR(ell);
    }

    if (GHC) {
        for (int i=0; i < nCurrent; i++)
            population[i].GHC();
    }

    for (int i=0; i<nCurrent; ++i) {
        Chromosome ch = population[i];
        orig_popu.push_back(ch);
        double f = population[i].getFitness();
        pHash[population[i].getKey()] = f;
        pHashOrig[population[i].getKey()] = f;
    }

}


DSMGA2::~DSMGA2 () {
    delete []masks;
    // delete []orig_masks;
    delete []orderELL;
    delete []fastCounting;
    delete []orig_fc;
    delete []usedIndex;
}



bool DSMGA2::isSteadyState () {

    if (stFitness.getNumber () <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean ()) {
        previousFitnessMean = stFitness.getMean () + EPSILON;
        return false;
    }

    return true;
}

bool DSMGA2::converged() {
    if (stFitness.getMax() == lastMax &&
        stFitness.getMean() == lastMean &&
        stFitness.getMin() == lastMin)
        convergeCount++;
    else
        convergeCount = 0;

    lastMax = stFitness.getMax();
    lastMean = stFitness.getMean();
    lastMin = stFitness.getMin();
    return (convergeCount > 300) ? true : false;
}

int DSMGA2::doIt (bool output) {
    generation = 0;
    while (!shouldTerminate ()) {
        oneRun (output);
        #ifdef DEBUG
        cin.get();
        #endif
    }
    return generation;
}


void DSMGA2::oneRun (bool output) {

    if (CACHE)
        Chromosome::cache.clear();

    mixing();


    double max = -INF;
    stFitness.reset ();

    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);

    }

    if (output)
        showStatistics ();

    ++generation;
}


bool DSMGA2::shouldTerminate () {
    bool  termination = false;

    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }

    if (population[0].getMaxFitness() <= stFitness.getMax() )
        termination = true;


    if (stFitness.getMax() - EPSILON <= stFitness.getMean() )
        termination = true;

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {

    printf ("Gen: %d, NFE: %d  Fitness:(Max/Mean/Min):%f/%f/%f \n ",
            generation, Chromosome::nfe, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
    printf ("best chromosome:");
    population[bestIndex].printOut();
    printf ("\n");


    fflush(NULL);
}



void DSMGA2::buildFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
    }

}

int DSMGA2::countOne(int x) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;


        val = fastCounting[x].gene[i];

        val ^= fastCounting[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}

//2016-03-09
// Almost identical to DSMGA2::findClique
// except check 00 or 01 before adding connection
void DSMGA2::findMask(Chromosome& ch, list<int>& result,int startNode){
    result.clear();

    
	DLLA rest(ell);
	genOrderELL();
	for( int i = 0; i < ell; i++){
		if(orderELL[i] == startNode)
			result.push_back(orderELL[i]);
		else
			rest.insert(orderELL[i]);
	}

	double *connection = new double[ell];

	for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
	    pair<double, double> p = graph(startNode, *iter);
		int i = ch.getVal(startNode);
		int j = ch.getVal(*iter);
		if(i == j)//p00 or p11
			connection[*iter] = p.first;
		else      //p01 or p10
			connection[*iter] = p.second;
	}
   
    while(!rest.isEmpty()){

	    double max = -INF;
		int index = -1;
		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
		    if(max < connection[*iter]){
			    max = connection[*iter];
				index = *iter;
			}
		}

		rest.erase(index);
		result.push_back(index);

		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
			pair<double, double> p = graph(index, *iter);
			int i = ch.getVal(index);
			int j = ch.getVal(*iter);
			if(i == j)//p00 or p11
				connection[*iter] += p.first;
			else      //p01 or p10
				connection[*iter] += p.second;
		}
	}

	delete []connection;
  
}

void DSMGA2::restrictedMixing(Chromosome& ch, int startNode) {
    
    if (startNode == -1)
        startNode = myRand.uniformInt(0, ell - 1);


    list<int> mask;
	findMask(ch, mask,startNode);
    size_t size = findSize(ch, mask);
   
    list<int> mask_size; 
    findMask_size(ch,mask_size,startNode,size);
    size_t size_original = findSize(ch,mask_size);

    if (size > size_original)
        size = size_original;
    while (mask.size() > size)
        mask.pop_back();
   
    bool taken = restrictedMixing(ch, mask);


    EQ = true;
    if (taken) {
        for (auto i: mask)
            usedIndex[i] = true;
    
        genOrderN();

        for (int i=0; i<nCurrent; ++i) {

            if (EQ)
                backMixingE(ch, mask, population[orderN[i]]);
            else
                backMixing(ch, mask, population[orderN[i]]);
        }

        // BMhistory.push_back(BMRecord(ch, mask, EQ, (double)sum/(double)nCurrent));
        BMhistory.push_back(BMRecord(ch, mask, EQ, 0));
    }

}
void DSMGA2::findMask_size(Chromosome& ch, list<int>& result,int startNode,int bound){
    result.clear();

    
	DLLA rest(ell);

	for( int i = 0; i < ell; i++){
		if(orderELL[i] == startNode)
			result.push_back(orderELL[i]);
		else
			rest.insert(orderELL[i]);
	}

	double *connection = new double[ell];

	for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
	    pair<double, double> p = graph_size(startNode, *iter);
		int i = ch.getVal(startNode);
		int j = ch.getVal(*iter);
		if(i == j)//p00 or p11
			connection[*iter] = p.first;
		else      //p01 or p10
			connection[*iter] = p.second;
	}
    bound--;
    while(!rest.isEmpty()&&bound>0){
        bound--;
	    double max = -INF;
		int index = -1;
		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
		    if(max < connection[*iter]){
			    max = connection[*iter];
				index = *iter;
			}
		}

		rest.erase(index);
		result.push_back(index);

		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
			pair<double, double> p = graph_size(index, *iter);
			int i = ch.getVal(index);
			int j = ch.getVal(*iter);
			if(i == j)//p00 or p11
				connection[*iter] += p.first;
			else      //p01 or p10
				connection[*iter] += p.second;
		}
	}

    delete []connection;
  
}

void DSMGA2::backMixing(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
        ++(des.layer);
          
        return;
    }

}

void DSMGA2::backMixingE(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        EQ = false;
        des = trial;
        ++(des.layer);

        return;
    }

    //2016-10-21
    if (trial.getFitness() >= des.getFitness() - EPSILON) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        des = trial;
        // ++(des.layer);

        return;
    }

}

bool DSMGA2::restrictedMixing(Chromosome& ch, list<int>& mask) {

    bool taken = false;
    size_t lastUB = 0;

    for (size_t ub = 1; ub <= mask.size(); ++ub) {

        size_t size = 1;
        // Chromosome trial(ell);
        Chromosome trial = ch;
	    
		//2016-03-03
	    vector<int> takenMask;

        for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
            
            //2016-03-03
			takenMask.push_back(*it);

            trial.flip(*it);

            ++size;
            if (size > ub) break;
        }

        //if (isInP(trial)) continue;
        //2016-10-21
        if (isInP(trial)) break;

        if (trial.getFitness() >= ch.getFitness() - EPSILON) {
            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();

            taken = true;
            ch = trial;
            ++(ch.layer);
        }

        if (taken) {
            lastUB = ub;
            break;
        }
    }

    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }

    return taken;

}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask, Chromosome& ch2) const {

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        if (ch.getVal(*it) == ch2.getVal(*it)) break;
        ++size;
    }
    return size;
}

void DSMGA2::mixing() {

    if (SELECTION)
        selection();

    buildFastCounting();
    buildGraph();
    buildGraph_sizecheck();

    genOrderELL();
    for (int i=0; i<ell; ++i)
        usedIndex[i] = false;

    for (int i=0; i<ell; ++i) {
        // if (usedIndex[orderELL[i]])
        //     continue;
        restrictedMixing(population[nCurrent - 1], orderELL[i]);
    }

    while (true) {
        if (SELECTION)
            selection();

        buildFastCounting();
        buildGraph();
        buildGraph_sizecheck();

        // Assume no chromosome added or deleted during RM
        int *old = new int[nCurrent];
        for (int i=0; i<nCurrent; ++i)
            old[i] = population[i].layer;

        genOrderN();
        for (int i=0; i<nCurrent; ++i) {
            restrictedMixing(population[orderN[i]]);
        }

        int count = 0;
        for (int i=0; i<nCurrent; ++i)
            if (old[i] < population[i].layer)
                ++count;

        delete []old;
        if (Chromosome::hit || 2 * count < nCurrent)
            break;
    }

    increaseOne();
}

inline bool DSMGA2::isInP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline bool DSMGA2::isInOrigP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHashOrig.find(ch.getKey());
    return (it != pHashOrig.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent-1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void DSMGA2::buildGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;

            double linkage = computeMI(p00,p01,p10,p11);
            
            //2016-04-08_computeMI_entropy
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00*log(p00/p_0/p0_);
            if (p11 > EPSILON)
                linkage00 += p11*log(p11/p_1/p1_);
            if (p01 > EPSILON)
                linkage01 += p01*log(p01/p0_/p_1);
            if (p10 > EPSILON)
                linkage01 += p10*log(p10/p1_/p_0);
           
            if(Chromosome::nfe < 0){
                pair<double, double> p(linkage, linkage);
                graph.write(i, j, p);
            }
            else{
                pair<double, double> p(linkage00, linkage01);
                graph.write(i, j, p);
            }
				
        }
    }


    delete []one;

}
void DSMGA2::buildGraph_sizecheck() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;

            double linkage = computeMI(p00,p01,p10,p11);
            
            //2016-04-08_computeMI_entropy
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00*log(p00/p_0/p0_);
            if (p11 > EPSILON)
                linkage00 += p11*log(p11/p_1/p1_);
            if (p01 > EPSILON)
                linkage01 += p01*log(p01/p0_/p_1);
            if (p10 > EPSILON)
                linkage01 += p10*log(p10/p1_/p_0);
        
	
            pair<double, double> p(linkage, linkage);
            graph_size.write(i, j, p);
			
        }
    }


    delete []one;

}


// from 1 to ell, pick by max edge
void DSMGA2::findClique(int startNode, list<int>& result) {
    
   }

    
    double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1-p0;
    double q1 = 1-q0;

    double join = 0.0;
    if (a00 > EPSILON)
        join += a00*log(a00);
    if (a01 > EPSILON)
        join += a01*log(a01);
    if (a10 > EPSILON)
        join += a10*log(a10);
    if (a11 > EPSILON)
        join += a11*log(a11);

    double p = 0.0;
    if (p0 > EPSILON)
        p += p0*log(p0);
    if (p1 > EPSILON)
        p += p1*log(p1);


    double q = 0.0;
    if (q0 > EPSILON)
        q += q0*log(q0);
    if (q1 > EPSILON)
        q += q1*log(q1);

    return -p-q+join;

}


void DSMGA2::selection () {
    tournamentSelection ();
}


// tournamentSelection without replacement
void DSMGA2::tournamentSelection () {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;
    }
}

void DSMGA2::increaseOne() {
    double max = -INF;
    int bestIndex = 0;

    genOrderN();
    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[orderN[i]].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = orderN[i];
        }
    }

    double *one = new double[ell];
    for (int i=0; i<ell; ++i) {
        one[i] = 0.0;
        for (int j=0; j<nCurrent; ++j)
            //if (orig_popu[j].getVal(i)==1)
            if (population[bestIndex].getVal(i)==1)
                one[i] += 1.0;
        //one[i] = nCurrent*0.5;
    }

    Chromosome ch;

    do {
        ch.initR(ell);
        do {
            for (int i=0; i<ell; ++i)
                ch.setVal(i, myRand.flip(one[i]/(double)nCurrent *0.33333 + 0.33333)? 0:1);
        } while (isInOrigP(ch) || isInP(ch));
        //ch.setVal(i, myRand.flip((one[i]+1.0)/(double)(nCurrent+2.0))? 0:1);
        ch.GHC();
    } while (isInOrigP(ch) || isInP(ch));

    delete []one;

    ++nCurrent;
    ++nOrig;

    population.push_back(ch);

    orig_popu.push_back(population[nCurrent-1]);

    pHash[population[nCurrent-1].getKey()] = population[nCurrent-1].getFitness();
    pHashOrig[population[nCurrent-1].getKey()] = population[nCurrent-1].getFitness();

    selectionIndex.emplace_back();
    orig_selectionIndex.emplace_back();
    orderN.push_back(nCurrent-1);

    for (int i = 0; i < ell; i++) {
        fastCounting[i].init(nCurrent);
        orig_fc[i].init(nOrig);
    }

    int size = BMhistory.size();
    int *rrr = new int[size];
    myRand.uniformArray(rrr, size, 0, size-1);
    for (int i=0; i<size; ++i) {
        //int r = i;
        int r = rrr[i];
        if (BMhistory[r].eq)
            backMixingE(BMhistory[r].pattern, BMhistory[r].mask, population[nCurrent-1]);
        else
            backMixing(BMhistory[r].pattern, BMhistory[r].mask, population[nCurrent-1]);
    }

    //population[nCurrent-1].GHC();

    delete []rrr;


}

