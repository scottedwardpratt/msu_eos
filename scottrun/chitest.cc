#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

using namespace std;
using namespace NMSUPratt;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	double T0=0.10,T,delT=0.001,oldP=-1.0,P;
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
	T=0.150;
	IntHadronGas->CalcQuantities(T,rhoB,rhoQ,rhoS,sampler);
	printf("T=%g, rhoB=%g:\n",T,rhoB);
	cout << IntHadronGas->chi << endl;
	
	printf("--------------------\n");
	rhoB=8*0.16;
	rhoQ=0.4*rhoB;
	rhoS=0.0;	
	T=0.181837;
	//printf("Enter T: ");
	//scanf("%lf",&T);
	IntHadronGas->CalcQuantities(T,rhoB,rhoQ,rhoS,sampler);
	printf("T=%g, epsilon=%g, rhoB=%g:\n",T,IntHadronGas->epsilon, rhoB);
	cout << IntHadronGas->chi << endl;
	
	delete sampler;
	delete reslist;
	
	return 0;
}
