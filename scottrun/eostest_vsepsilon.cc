#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	double T0=0.190;
	double rhoB,rhoQ,rhoS,epsilon,etarget;
	epsilon=0.25;
	Crandy *randy=new Crandy(-1234);
	string parfilename="parameters/parameters.txt";
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(parfilename);
	CresList *reslist=new CresList(parmap);
	CinteractingHadronGas *IntHadronGas=new CinteractingHadronGas(parmap);
	Csampler *sampler=new Csampler(T0,0.093,parmap,reslist,randy);
	IntHadronGas->hgasinfo->sampler=sampler;
	
	//IntHadronGas->CalcQuantitiesVsEpsilon(epsilon,rhoB,rhoQ,rhoS);
	
	double deltau=0.002,tau0=1.0;
	for(double tau=1.0;tau<11.0+0.5*deltau;tau+=deltau){
		rhoB=8*0.16/tau;
		rhoQ=0.4*rhoB;
		rhoS=0.0;
		epsilon=3.0*pow(tau/tau0,-1.25);
		etarget=epsilon;
		IntHadronGas->CalcQuantitiesVsEpsilon(epsilon,rhoB,rhoQ,rhoS);
		
		printf("tau=%5.2f, T=%8.5f, e=%8.5f=?%8.5f, e_h=%8.5f, e_int=%8.5f, rhoB=%8.5f, cs2=%8.5f, P=%8.5f, s=%8.5f\n",
		tau,IntHadronGas->T,IntHadronGas->epsilon,etarget,
		IntHadronGas->hgasinfo->epsilon,IntHadronGas->hintinfo->epsilon,IntHadronGas->rhoB,IntHadronGas->cs2,IntHadronGas->P,IntHadronGas->s);
	}
	
	delete sampler;
	delete reslist;
	
	return 0;
}
