#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	MSU_EOS::USE_POLE_MASS=true;
	double T0=0.190;
	double rhoB,rhoQ,rhoS,epsilon;
	epsilon=0.25;
	rhoB=(8.0/11.0)*0.16;
	rhoQ=0.4*rhoB;
	rhoS=0.0;
	printf("rho=(%g,%g,%g)\n",rhoB,rhoQ,rhoS);
	Crandy *randy=new Crandy(-1234);
	string parfilename="parameters/parameters.txt";
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(parfilename);
	CresList *reslist=new CresList(parmap);
	CinteractingHadronGas *IntHadronGas=new CinteractingHadronGas(parmap);
	Csampler *sampler=new Csampler(T0,0.093,parmap,reslist,randy);
	IntHadronGas->hgasinfo->sampler=sampler;
	
	//IntHadronGas->CalcQuantitiesVsEpsilon(epsilon,rhoB,rhoQ,rhoS);
	
	double deltau=0.1,tau0=1.0;
	for(double tau=1.0;tau<11.0+0.5*deltau;tau+=deltau){
		rhoB=8*0.16/tau;
		rhoQ=0.4*rhoB;
		rhoS=0.0;
		epsilon=3.0*pow(tau/tau0,-1.25);
		printf("tau=%5g, epsilon=%g\n",tau,epsilon);
		IntHadronGas->CalcQuantitiesVsEpsilon(epsilon,rhoB,rhoQ,rhoS);
	}
	
	delete sampler;
	delete reslist;
	
	return 0;
}
