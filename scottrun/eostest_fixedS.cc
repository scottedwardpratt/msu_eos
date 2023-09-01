#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	double T0=0.15,T,delT=0.0001;
	double rhoB0,rhoB,rhoQ,rhoS,soverrhoB,starget,oldsoverrhoB,dS;
	double P,epsilon,s,chiEE,chiEB,chiEQ,chiES,chiPP,cs2;
	Eigen::Matrix3d chi;
	Crandy *randy=new Crandy(-1234);
	string parfilename="parameters/parameters.txt";
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(parfilename);
	CresList *reslist=new CresList(parmap);
	Csampler *sampler=new Csampler(T0,0.093,parmap,reslist,randy);
	//Csampler *sampler=0;
	
	CinteractingHadronGas *IntHadronGas=new CinteractingHadronGas(parmap);
	
	rhoB0=8*0.16/11.0;
	T=T0;
	rhoS=0.0;
	starget=22.5311;
	soverrhoB=starget;
	
	printf("   T   rhoB  epsilon    P   s/rhoB    cs2    chiEE/s  chiPP/s chiEB/s chiEQ/s chiES/s chiBB/s  chiQQ/s  chiSS/s\n");
	
	for(rhoB=rhoB0;rhoB<1.0;rhoB+=0.03){
		rhoQ=0.4*rhoB;
		
		do{
			oldsoverrhoB=soverrhoB;
			IntHadronGas->CalcQuantities(T,rhoB,rhoQ,rhoS,sampler);
			soverrhoB=IntHadronGas->s/rhoB;
			T+=delT;
		}while(soverrhoB<starget);
		dS=soverrhoB-oldsoverrhoB;
		T=(T-delT)*(soverrhoB-starget)/dS + T*(starget-oldsoverrhoB)/dS;
		
		//printf("oldsoverrhoB=%g, starget=%g, soverrhoB=%g, ds=%g=?%g\n",oldsoverrhoB,starget,soverrhoB,soverrhoB-oldsoverrhoB,dS);
		
		IntHadronGas->CalcQuantities(T,rhoB,rhoQ,rhoS,sampler);
		soverrhoB=IntHadronGas->s/rhoB;
		
		//printf("oldsoverrhoB=%g, starget=%g, soverrhoB=%g, ds=%g=?%g\n",
		//oldsoverrhoB,starget,soverrhoB,soverrhoB-oldsoverrhoB,dS);
		
		P=IntHadronGas->P;
		epsilon=IntHadronGas->epsilon;
		T=IntHadronGas->T;
		s=IntHadronGas->s;
		chiEE=IntHadronGas->chiEE;
		chiPP=(P+epsilon)/T;
		chiEB=IntHadronGas->chiEQ[0];
		chiEQ=IntHadronGas->chiEQ[1];
		chiES=IntHadronGas->chiEQ[2];
		chi=IntHadronGas->chi;
		cs2=IntHadronGas->cs2;
		
		printf("%6.4f %6.4f %6.4f %6.4f %7.4f %7.4f  %7.4f %7.4f %7.4f %7.4f %7.4f   %7.4f %7.4f %7.4f\n",
		T,rhoB,epsilon,P,s/rhoB,cs2,chiEE/s,chiPP/s,chiEB/s,chiEQ/s,chiES/s,chi(0,0)/s,chi(1,1)/s,chi(2,2)/s);
		
		//IntHadronGas->PrintQuantities();
	}
	
	delete sampler;
	delete reslist;
	
	return 0;
}
