#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	double T0=0.05,T,delT=0.01,oldP=-1.0,P;
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
	
	for(T=T0;T<0.171;T+=delT){
		oldP=-1.0;
		printf("--------- T=%g -----------\n",T);
		for(rhoB=0.05;rhoB<1.0;rhoB+=0.02){
		
			IntHadronGas->CalcQuantities(T,rhoB,rhoQ,rhoS,sampler);
			//IntHadronGas->PrintQuantities();
			P=IntHadronGas->P;
			if(P<oldP){
				printf("-------------------- unstable -------------------------\n");
			}
			printf("%7.4f %7.4f %7.4f %7.4f %7.4f\n",T,rhoB,P,IntHadronGas->cs2,(P+IntHadronGas->epsilon)*T/(IntHadronGas->s));
			oldP=P;
		}
	}
	
	delete sampler;
	delete reslist;
	
	return 0;
}
