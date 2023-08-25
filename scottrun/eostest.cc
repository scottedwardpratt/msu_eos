#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	double T0=0.10;
	double rhoB,rhoQ,rhoS;
	Crandy *randy=new Crandy(-1234);
	string parfilename="parameters/parameters.txt";
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(parfilename);
	CresList *reslist=new CresList(parmap);
	CinteractingHadronGas *IntHadronGas=new CinteractingHadronGas(parmap);
	Csampler *sampler=new Csampler(T0,0.093,parmap,reslist,randy);
	IntHadronGas->hgasinfo->sampler=sampler;
	
	//IntHadronGas->CalcQuantitiesVsEpsilon(epsilon,rhoB,rhoQ,rhoS);
	
	cout << IntHadronGas->hgasinfo->sampler << endl;
	
	double T,delT=0.01;
	for(T=0.10;T<0.171;T+=delT){
		rhoB=8*0.16/11.0;
		rhoQ=0.4*rhoB;
		rhoS=0.0;
		printf("T=%g\n",T);
		cout << IntHadronGas->hgasinfo->sampler << endl;
		IntHadronGas->CalcQuantities(T,rhoB,rhoQ,rhoS);
		IntHadronGas->PrintQuantities();
	}
	
	delete sampler;
	delete reslist;
	
	return 0;
}
