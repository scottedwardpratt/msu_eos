#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	MSU_EOS::USE_POLE_MASS=true;
	double T0=0.140;
	double rhoB,rhoQ,rhoS,epsilon;
	epsilon=0.25;
	rhoB=(8.0/11.0)*0.16;
	rhoQ=0.4*rhoB;
	rhoS=0.0;
	Crandy *randy=new Crandy(-1234);
	string parfilename="parameters/parameters.txt";
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(parfilename);
	CresList *reslist=new CresList(parmap);
	CinteractingHadronGas *IntHadronGas=new CinteractingHadronGas(parmap);
	Csampler *sampler=new Csampler(T0,0.093,parmap,reslist,randy);
	IntHadronGas->hgasinfo->sampler=sampler;
	
	//ChIntInfo_Scott *int_scott=new ChIntInfo_Scott(parmap); // note this does (must) not depend on T
	//int_scott->CalcQuantities(T0,rhoB,rhoQ,rhoS);
	
	IntHadronGas->CalcQuantitiesVsEpsilon(epsilon,rhoB,rhoQ,rhoS);
	

	
	
	/*
	IntHadronGas->hgasinfo->sampler=sampler;
	
	string filename;
	FILE *fptr;
	for(T=T0;T<0.161;T+=0.01){
	filename="results/T"+to_string(lrint(1000*T))+".txt";
	fptr=fopen(filename.c_str(),"w");
	printf("---------- T=%g ------------\n",T);
	printf("   rho.    cs2.      P.       f.      mu/T.   chiinv(0,0)\n");
	fprintf(fptr,"   rho.    cs2.      P.       f.      mu/T.   chiinv(0,0)\n");
	for(rhoB=0.0;rhoB<1.0;rhoB+=0.02){
	rhoS=0.0;
	rhoQ=0.4*rhoB;
	IntHadronGas->CalcQuantities(T,rhoB,rhoQ,rhoS);
	printf("%6.3f  %7.4f  %7.4f  %7.4f   %7.4f   %7.4f\n",
	rhoB,IntHadronGas->cs2,IntHadronGas->P,
	IntHadronGas->hintinfo->f,IntHadronGas->hintinfo->muB,IntHadronGas->hintinfo->chiinv(0,0));
	fprintf(fptr,"%6.3f  %7.4f  %7.4f  %7.4f   %7.4f   %7.4f\n",
	rhoB,IntHadronGas->cs2,IntHadronGas->P,
	IntHadronGas->hintinfo->f,IntHadronGas->hintinfo->muB,IntHadronGas->hintinfo->chiinv(0,0));
	//IntHadronGas->PrintQuantities();
	}
	fclose(fptr);
	}
		
	*/
		
	delete sampler;
	delete reslist;
	
	return 0;
}
