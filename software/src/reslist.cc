#include "msu_eos/resonances.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/log.h"
using namespace NMSUPratt;

double CresList::MIN_DECAY_WIDTH=0.0001;
char *CresList::message=new char[CLog::CHARLENGTH];

CresList::CresList(){
}

CresList::~CresList(){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos=resmap.begin();
	while(rpos!=resmap.end()){
		resinfo=rpos->second;
		delete resinfo;
		resmap.erase(rpos);
		rpos=resmap.begin();
	}
	resmap.clear();
	massmap.clear();
}

CresInfo* CresList::GetResInfoPtr(int pid){
	CresInfoMap::iterator rpos;
	rpos=resmap.find(pid);
	if(rpos!=resmap.end())
		return rpos->second;
	else{
		snprintf(message,CLog::CHARLENGTH,"GetResInfoPtr() can't find match for PID=%d\n",pid);
		CLog::Fatal(message);
		return NULL;
	}
}

void CresList::CalcMinMasses(){
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		resinfo->CalcMinMass();
	}
}

CresList::CresList(CparameterMap* parmap_in){
	parmap=parmap_in;
	CresInfo::NSPECTRAL=parmap->getI("MSU_SAMPLER_NSPECTRAL",100);
	CresInfo::SFDIRNAME=parmap->getS("MSU_SAMPLER_SFDIRNAME","../progdata/resinfo/spectralfunctions");
	IGNORE_CHARM_BOTTOM=parmap->getB("MSU_SAMPLER_IGNORE_CHARM_BOTTOM",true);
	//RESONANCE_DECAYS=parmap->getB("RESONANCE_DECAYS",true);
	readmsu=parmap->getB("MSU_SAMPLER_READ_MSU",false);
	CresInfo::reslist=this;
	if(readmsu){
		ReadResInfo_MSU();
	}
	else{
		ReadResInfo();
	}
	//CalcSpectralFunctions();
	if(!readmsu){
		ReadSpectralFunctions();
	}
	CalcMinMasses();
}

void CresList::CalcSpectralFunctions(){
	CresMassMap::iterator rpos;
	CresInfo *resinfo;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->decay && !resinfo->SFcalculated){
			resinfo->CalcSpectralFunction();
		}
		else
			resinfo->SFcalculated=true;
	}
}

void CresList::PrintMassMaps(){
	CresMassMap::iterator rpos;
	map<double,double>::iterator it;
	CresInfo *resinfo;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->decay){
			it=resinfo->sfmassmap.begin();
			snprintf(message,CLog::CHARLENGTH," ----- SF massmap for pid=%d ----- \n",resinfo->pid);
			CLog::Info(message);
			while(it!=resinfo->sfmassmap.end()){
				snprintf(message,CLog::CHARLENGTH,"%g   %g\n",it->first,it->second);
				CLog::Info(message);
				it++;
			}
		}
	}
}
