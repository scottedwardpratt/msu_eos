#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	double T0=0.10,T,delT=0.01;
	double rhoB,rhoQ,rhoS;
	Crandy *randy=new Crandy(-1234);
	string parfilename="parameters/parameters.txt";
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(parfilename);
	CresList *reslist=new CresList(parmap);
	Csampler *sampler=new Csampler(T0,0.093,parmap,reslist,randy);
	//Csampler *sampler=0;
	
	CinteractingHadronGas *IntHadronGas=new CinteractingHadronGas(parmap);
	
	rhoB=8*0.16/11.0;
	rhoQ=0.4*rhoB;
	rhoS=0.0;	
	
	for(T=0.10;T<0.171;T+=delT){
		
		IntHadronGas->CalcQuantities(T0,rhoB,rhoQ,rhoS,sampler);
		IntHadronGas->PrintQuantities();
	}
	
	delete sampler;
	delete reslist;
	
	return 0;
}
