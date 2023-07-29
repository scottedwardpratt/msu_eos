#include "msu_eos/resonances.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/log.h"
#include "msu_sampler/sampler.h"
#include "msu_eos/eos.h"

CcanonicalHadronGasInfo::CcanonicalHadronGasInfo(){
	muB=muS=muQ=0.0;
	sampler=NULL;
	chi.resize(3,3);
	chiinv.resize(3,3);
	A.resize(4,4);
	chiEQ.resize(3);
}

CcanonicalHadronGasInfo::~CcanonicalHadronGasInfo(){
	chi.resize(0,0);
	chiinv.resize(0,0);
	A.resize(0,0);
	chiEQ.resize(0);
}

void CcanonicalHadronGasInfo::CalcQuantities(double Tset,double rhoBset,double rhoQset,double rhoSset){
	T=Tset;
	rhoB=rhoBset;
	rhoQ=rhoQset;
	rhoS=rhoSset;
	double rhoII,rhoIIcheck,rhoBcheck,rhoScheck,muII;
	if(sampler==NULL){
		CLog::Fatal("In CcanonicalHadronGasInfo::CalcHadronicQuantities, sampler pointer is still NULL\n");
	}
	rhoII=2*rhoQ-rhoB-rhoS;
	if(fabs(sampler->Tf-T)>0.0001 || !sampler->forMU0_calculated){
		//CLog::Info("Resetting sampler T in CcanonicalHadronGasInfo to "+to_string(T)+"\n");
		sampler->Tf=T;
		sampler->GetNHMu0();
	}
	sampler->GetMuNH(rhoB,rhoII,rhoS,muB,muII,muS,nhadrons);
	//sampler->GetEpsilonRhoChi(muB,muII,muS,epsilon,rhoBcheck,rhoIIcheck,rhoScheck,chi);
	sampler->GetEpsilonRhoDerivatives(muB,muII,muS,epsilon,rhoBcheck,rhoIIcheck,rhoScheck,A);
	double etest=0.0;
	sampler->CalcNHadronsEpsilonP(muB,muII,muS,nhadrons,etest,P);
	s=(etest+P)/T-muB*rhoB-muII*rhoII-muS*rhoS;
	f=etest-T*s;
	
	chi(0,0)=A(1,1);
	chi(0,1)=0.5*(A(1,1)+A(1,2)+A(1,3));
	chi(0,2)=A(1,3);
	
	chi(1,0)=chi(0,1);
	chi(1,1)=0.25*(A(1,1)+A(2,2)+A(3,3))+0.5*(A(1,2)+A(1,3)+A(2,3));
	chi(1,2)=0.5*(A(1,3)+A(2,3)+A(3,3));
	
	chi(2,0)=chi(0,2);
	chi(2,1)=chi(1,2);
	chi(2,2)=A(3,3);
	
	chiEE=A(0,0);
	chiEQ(0)=A(0,1);
	chiEQ(1)=0.5*(A(0,1)+A(0,2)+A(0,3));
	chiEQ(2)=A(0,3);
	
	muQ=2*muII;
	muB=-muII+muB;
	muS=-muII+muS;
	
	chiinv=chi.inverse();
	
}

CinteractingHadronGas::CinteractingHadronGas(CparameterMap *parmap_set){
	parmap=parmap_set;
	muB=muS=muQ=0.0;
	chi.resize(3,3);
	chiinv.resize(3,3);
	chiEQ.resize(3);
	dPdrho_T.resize(3);
	dedrho_T.resize(3);
	dPdrho_e.resize(3);
	hgasinfo=new CcanonicalHadronGasInfo();
	//hintinfo=new ChadronInteractionInfo();
	hintinfo=new ChIntInfo_Scott(parmap);
}

CinteractingHadronGas::~CinteractingHadronGas(){
	chi.resize(0,0);
	chiinv.resize(0,0);
	chiEQ.resize(0);
	//delete chi;
	//delete 	chiinv;
	//delete chiEQ;
	//delete hgasinfo;
	//delete htinfo;
}

void CinteractingHadronGas::CalcQuantities(double Tset,double rhoBset,double rhoQset,double rhoSset){
	double dedt_rho;
	int a,b;
	T=Tset;
	rhoB=rhoBset;
	rhoQ=rhoQset;
	rhoS=rhoSset;
	hgasinfo->CalcQuantities(T,rhoB,rhoQ,rhoS);
	hintinfo->CalcQuantities(T,rhoB,rhoQ,rhoS);
	chiinv=hgasinfo->chiinv+hintinfo->chiinv;
	chi=chiinv.inverse();
	muB=hgasinfo->muB+hintinfo->muB;
	muQ=hgasinfo->muQ+hintinfo->muQ;
	muS=hgasinfo->muS+hintinfo->muS;
	epsilon=hgasinfo->epsilon+hintinfo->epsilon;
	P=hgasinfo->P+hintinfo->P;
	s=hgasinfo->s+hintinfo->s;
	f=hgasinfo->f+hintinfo->f;
	Eigen::MatrixXd dedrho(3,3);
	
	dedrho=hgasinfo->chiinv*hgasinfo->chiEQ;
	dedrho=dedrho+hintinfo->dedrho;
	chiEQ=chi*dedrho;
	
	dedt_rho=hgasinfo->chiEE-hgasinfo->chiEQ.transpose()*hgasinfo->chiinv*hgasinfo->chiEQ;
	dedt_rho+=hintinfo->dedT;
	chiEE=dedt_rho+chiEQ.transpose()*chiinv*chiEQ;
	
	//chiEE=hgasinfo->chiEE+hintinfo->dedT*T*T;
	//chiEE-=T*T*hintinfo->dedrho.transpose()*chiinv*hintinfo->dmudT;
	
	// Calculate speed of sound
	Eigen::Matrix3d A(3,3),Ainv(3,3);
	Eigen::Vector3d rho(3),mu(3),chimu(3),M(3);

	rho(0)=rhoB;
	rho(1)=rhoQ;
	rho(2)=rhoS;
	mu(0)=muB;
	mu(1)=muQ;
	mu(2)=muS;
	A=s*chi;
	chimu=chi*mu;
	for(a=0;a<3;a++){
		for(b=0;b<3;b++){
			A(a,b)+=rho(a)*chimu(b);
			A(a,b)-=(1.0/T)*rho(a)*chiEQ(b);
		}
	}
	Ainv=A.inverse();
	double muzeta=chiEQ.transpose()*mu;
	M=Ainv*( (chiEE/(T*T)-muzeta/T)*rho -(s/T)*chiEQ);
	//for(a=0;a<3;a++){
	//	M(a)=0.0;
	//	for(b=0;b<3;b++){
		//	M(a)+=Ainv(a,b)*( rho(b)*((chiEE/(T*T))-muzeta/T) -(s/T)*chiEQ(b) );
	//	}
	//}
	double rhoM,chiEQM;
	rhoM=rho.transpose()*M;
	chiEQM=chiEQ.transpose()*M;
	cs2=((P+epsilon)/T+rhoM)/(chiEE/(T*T)+chiEQM/T);
	
	dPdT_rho=(P+epsilon)/T-(1.0/T)*rho.transpose()*chiinv*chiEQ;
	dedT_rho=chiEE/(T*T)-(1.0/(T*T))*chiEQ.transpose()*chiinv*chiEQ;
	dPde_rho=dPdT_rho/dedT_rho;
	dPdrho_T=chiinv*rho*T;
	dedrho_T=chiinv*chiEQ;
	dPdrho_e=dPdrho_T-(dPdT_rho/dedT_rho)*dedrho_T;
	
}

void CinteractingHadronGas::CalcQuantitiesVsEpsilon(double epsilontarget,double rhoBtarget,double rhoQtarget,double rhoStarget){
	char message[CLog::CHARLENGTH];
	//4D Newton's Method
	// Here rhoII refers to rho_u-rho_d = 2*I3 and mu[1]=muII/2
	double rhoII,rhoIItarget,smb,cmb,Tf,epsilonh,epsilon,muII_h,muB_h,muS_h;
	Eigen::MatrixXd A(4,4);
	Eigen::VectorXd dmu(4),drho(4);
	int ntries=0;
	
	rhoIItarget=2*rhoQtarget-rhoBtarget-rhoStarget;
	// change basis of chemical potentials
	muII_h=0.5*hgasinfo->muQ;
	muB_h=hgasinfo->muB+muII_h;
	muS_h=hgasinfo->muS+muII_h;
	
	Csampler *sampler=hgasinfo->sampler;
	Tf=sampler->Tf;
	sampler->GetNHMu0();
	
	do{
		ntries+=1;
		if(ntries>30){
			snprintf(message,CLog::CHARLENGTH,"FAILURE, ntries=%d\n",ntries);
			CLog::Fatal(message);
		}
		smb=sinh(muB_h);
		cmb=cosh(muB_h);
		
		hintinfo->CalcQuantities(Tf,rhoBtarget,rhoQtarget,rhoStarget);
		sampler->GetEpsilonRhoDerivatives(muB_h,muII_h,muS_h,epsilonh,rhoB,rhoII,rhoS,A);
		for(int i=0;i<4;i++){
			A(i,1)=A(i,1)/cmb;
			A(i,0)=A(i,0)/(Tf*Tf);
		}
		epsilon=epsilonh+hintinfo->epsilon;
		A(0,0)+=hintinfo->dedT;
		
		drho[0]=epsilontarget-epsilon;
		drho[1]=rhoBtarget-rhoB;
		drho[2]=rhoIItarget-rhoII;
		drho[3]=rhoStarget-rhoS;
		dmu=A.colPivHouseholderQr().solve(drho);
		Tf+=dmu[0];
		sampler->Tf=Tf;
		sampler->GetNHMu0();
		smb+=dmu[1];
		muB_h=asinh(smb);
		muII_h+=dmu[2];
		muS_h+=dmu[3];

	}while(fabs(drho[0])>1.0E-5 || fabs(drho[1])>1.0E-6 || fabs(drho[2])>1.0E-6 || fabs(drho[3])>1.0E-6);
	
	epsilon=epsilontarget;
	rhoB=rhoBtarget;
	rhoQ=rhoQtarget;
	rhoS=rhoStarget;
	T=Tf;
	// change back to BQS basis
	hgasinfo->muQ=2*muII_h;
	hgasinfo->muB=-muII_h+muB_h;
	hgasinfo->muS=-muII_h+muS_h;
	
	CalcQuantities(T,rhoB,rhoQ,rhoS);
	
	//snprintf(message,CLog::CHARLENGTH,"TEST: T=%g, epsilon=%g, rho=(%g,%g,%g)\n",T,epsilon,rhoB,rhoQ,rhoS);
	//CLog::Info(message);
}
	
void CinteractingHadronGas::PrintQuantities(){
	int a;
	char output[CLog::CHARLENGTH];
	snprintf(output,CLog::CHARLENGTH,"T=%g, rhoB=%g, rhoQ=%g, rhoS=%g\n",T,rhoB,rhoQ,rhoS);
	CLog::Info(output);
	snprintf(output,CLog::CHARLENGTH,"P=%g, epsilon=%g, f=%g, s=%g, c_s^2=%g\n",P,epsilon,f,s,cs2);
	CLog::Info(output);
	snprintf(output,CLog::CHARLENGTH,"chi=\n");
	CLog::Info(output);
	for(a=0;a<3;a++){
		snprintf(output,CLog::CHARLENGTH,"%10.3e %10.3e %10.3e\n",chi(a,0),chi(a,1),chi(a,2));
		CLog::Info(output);
	}
}
