#include "msu_eos/resonances.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/log.h"
using namespace std;
using namespace NMSUPratt;

void CresList::ReadResInfo(){
	//Cmerge *merge;
	int motherpid,pid;
	bool nophotons,ignore_resonance;
	double bmax,bsum;
	int ires,ichannelres,ichannel,ibody,nbodies,NResonances,LDecay=1;
	int degenread,baryonread,sread,cread,bread,Iread,qread,nchannelsread;
	double mread,wread;
	int ires1,ires2,iresflip,Nu,Nd,Ns,Nc,Nb;
	//int netq,netb,nets;
	string name, filename;
	CresInfo *resinfo=NULL,*aresinfo=NULL,*temp=NULL;
	CdecayInfo *decayinfo, *adecayinfo;
	CbranchInfo *bptr=NULL,*firstbptr=NULL;
	FILE *resinfofile;
	char cname[200],dummy[200];
	CdecayInfoMap decaymap;
	int dummy_int;
	CresInfoMap::iterator iter;
	CdecayInfoMap::iterator diter;
	Cmerge *merge;
	filename=parmap->getS("RESONANCES_INFO_FILE",string("../software/resinfo/pdg-SMASH.dat"));
	snprintf(message,CLog::CHARLENGTH,"will read resonance info from %s\n",filename.c_str());
	resinfofile=fopen(filename.c_str(),"r");
	CLog::Info("will read resonance info from "+filename+"\n");
	if (resinfofile==NULL) {
		snprintf(message,CLog::CHARLENGTH,"Can't open resinfofile\n");
		CLog::Fatal(message);
	}

	ires=0;
	NResonances=0;
	while(fscanf(resinfofile," %d",&pid)!=EOF && NResonances<20000){
		fscanf(resinfofile,"%s %lf %lf %d %d %d %d %d %d %d %d",
		cname,&mread,&wread,&degenread,&baryonread,&sread,&cread,&bread,&Iread,&qread,&nchannelsread);
		fgets(dummy,200,resinfofile);
		if(cread==0 && bread==0     // neglect any charmed or bottomed hadrons
			&& pid!=441 && pid!=443 && pid!=445 && pid!=100441 && pid!=100441
				&& pid!=551 && pid!=553 && pid!=555 && pid!=100553 && pid!=200553
		&& abs(baryonread)<2){
			
			resinfo=new CresInfo();
			NResonances+=1;
			resinfo->pid=pid;
			resinfo->mass=mread;
			resinfo->width=wread;
			resinfo->degen=degenread;
			resinfo->baryon=baryonread;
			resinfo->strange=sread;
			resinfo->charm=cread;
			resinfo->bottom=bread;
			resinfo->total_isospin=Iread;
			resinfo->charge=qread;
			resinfo->nchannels=nchannelsread;
			
			//fscanf(resinfofile,"%s %lf %lf %d %d %d %d %d %d %d %d", cname,&resinfo->mass,&resinfo->width,&resinfo->degen,&resinfo->baryon,&resinfo->strange,&resinfo->charm,&resinfo->bottom,&resinfo->total_isospin,&resinfo->charge,&resinfo->nchannels);
			resinfo->minmass=resinfo->mass;
			if(resinfo->width<MIN_DECAY_WIDTH){
				resinfo->decay=false;
				decayinfo=NULL;
			}
			else{
				resinfo->decay=true;
				decayinfo=new CdecayInfo();
			}
			cname[int(strlen(cname))-1]='\0';
			resinfo->name=cname;
			resinfo->q[0]=resinfo->baryon+resinfo->charge-resinfo->charm+resinfo->bottom;
			resinfo->q[1]=2.0*resinfo->baryon+resinfo->strange-resinfo->charge;
			resinfo->q[2]=-resinfo->strange;
			
			if(abs(resinfo->baryon)==1){
				Nu=abs(resinfo->q[0]);
				Nd=abs(resinfo->q[1]);
				Ns=abs(resinfo->q[2]);
				Nc=abs(resinfo->charm);
				Nb=abs(resinfo->bottom);
			}
			else{
				if(resinfo->q[0]==0 && resinfo->q[1]==0 && resinfo->strange==0 && resinfo->charm==0 && resinfo->bottom==0){
					if(pid==331 || pid==10331 || pid==100331 || pid==333 || 100333 || 337){
						Ns=2; Nu=Nd=Nc=Nb=0;
					}
					if(pid==443 || pid==10441 || pid==100441){
						Nc=2; Nu=Nd=Ns=Nb=0;
					}
					if(pid==551 || pid==553){
						Nb=2; Nu=Nd=Ns=Nc=0;
					}
				}
				else{
					Nu=abs(resinfo->q[0]);
					Nd=abs(resinfo->q[1]);
					Ns=abs(resinfo->q[2]);
					Nc=abs(resinfo->charm);
					Nb=abs(resinfo->bottom);
				}
			}
			resinfo->Nu=Nu;
			resinfo->Nd=Nd;
			resinfo->Ns=Ns;
			resinfo->Nc=Nc;
			resinfo->Nb=Nb;
			//
			//
			//decay reading
			//reads into map values: will access for decays when done creating resonances
			for (ichannel=0; ichannel<resinfo->nchannels; ichannel++){
				if(resinfo->decay){
					fscanf(resinfofile, "%d %d %lf %d %d %d %d %d %d", 
					&dummy_int,&decayinfo->Nparts[ichannel],&decayinfo->branchratio[ichannel],&decayinfo->products[ichannel][0],
					&decayinfo->products[ichannel][1],&decayinfo->products[ichannel][2],&decayinfo->products[ichannel][3],
					&decayinfo->products[ichannel][4],&decayinfo->d_L[ichannel]);
				}
				else{
					fscanf(resinfofile,"%d",&dummy_int);
					fgets(dummy,200,resinfofile);
				}
			}
			if(resinfo->decay){
				bsum=0.0;
				for(ichannel=0;ichannel<resinfo->nchannels;ichannel++)
					bsum+=decayinfo->branchratio[ichannel];
				for(ichannel=0;ichannel<resinfo->nchannels;ichannel++)
					decayinfo->branchratio[ichannel]=decayinfo->branchratio[ichannel]/bsum;
			}

			if(resinfo->pid!=22){ //copied from old pid
				resinfo->ires=ires;
				ires+=1;
			}
			if(resinfo->baryon!=0){
				resinfo->SetBtype();
			}
			resinfo->branchlist.clear();

			ignore_resonance=false;
			if(!IGNORE_CHARM_BOTTOM || (resinfo->Nc==0 && resinfo->Nb==0)){
				resmap.insert(CresInfoPair(resinfo->pid,resinfo));
				massmap.insert(CresMassPair(resinfo->mass,resinfo));
				if(resinfo->decay){
					if(!IGNORE_CHARM_BOTTOM || (resinfo->Nc==0 && resinfo->Nb==0)){
						decaymap.insert(CdecayInfoPair(resinfo->pid,decayinfo));
					}
				}
			}
			else{
				ignore_resonance=true;
				if(resinfo->decay){
					delete decayinfo;
				}
				delete resinfo;
			}
			
			//antiparticle creation

			if(!ignore_resonance){
				if(resinfo->baryon!=0){
					aresinfo=new CresInfo();
					NResonances+=1;

					aresinfo->pid=-resinfo->pid;
					aresinfo->mass=resinfo->mass;
					aresinfo->minmass=resinfo->minmass;
					aresinfo->width=resinfo->width;
					aresinfo->degen=resinfo->degen;
					aresinfo->baryon=-resinfo->baryon;
					aresinfo->strange=-resinfo->strange;
					aresinfo->charm=-resinfo->charm;
					aresinfo->bottom=-resinfo->bottom;
					aresinfo->charge=-resinfo->charge;
					aresinfo->decay=resinfo->decay;
					aresinfo->nchannels=resinfo->nchannels;
					aresinfo->q[0]=-resinfo->q[0];
					aresinfo->q[1]=-resinfo->q[1];
					aresinfo->q[2]=-resinfo->q[2];
					cname[int(strlen(cname))-1]='\0';
					string s(cname);
					aresinfo->name="Anti-"+s;
					aresinfo->ires=ires;
					ires+=1;
					resmap.insert(CresInfoPair(aresinfo->pid,aresinfo));
					massmap.insert(CresMassPair(aresinfo->mass,aresinfo));
					if(aresinfo->decay)
						adecayinfo=new CdecayInfo();
					else
						adecayinfo=NULL;

					aresinfo->branchlist.clear();
					if(aresinfo->decay){
						for(ichannel=0; ichannel<resinfo->nchannels; ichannel++) { //reads into map values: will access for decays when done creating resonances
							for(int i=0; i<decayinfo->Nparts[ichannel]; i++) {
								pid=decayinfo->products[ichannel][i];
								if(pid!=0){
									temp=GetResInfoPtr(pid);
									if(temp->baryon==0 && temp->charge==0 && temp->strange==0){
										adecayinfo->products[ichannel][i]=decayinfo->products[ichannel][i];
									}
									else{
										adecayinfo->products[ichannel][i]=-decayinfo->products[ichannel][i];
									}
								}
								else adecayinfo->products[ichannel][i]=0;
							}
							adecayinfo->Nparts[ichannel]=decayinfo->Nparts[ichannel];
							adecayinfo->branchratio[ichannel]=decayinfo->branchratio[ichannel];
							adecayinfo->d_L[ichannel]=decayinfo->d_L[ichannel];
						}
						decaymap.insert(CdecayInfoPair(aresinfo->pid,adecayinfo));
					}
					aresinfo->Nu=Nu;
					aresinfo->Nd=Nd;
					aresinfo->Ns=Ns;
					aresinfo->Nc=Nc;
					aresinfo->Nb=Nb;
				}
			}
		}
		else{
			for(ichannel=0;ichannel<nchannelsread;ichannel++){
				fgets(dummy,200,resinfofile);
			}
		}
	}
	fclose(resinfofile);
	CLog::Info("NResonances:"+to_string(NResonances)+"\n");
	//------------------------------------------
	MergeArray=new Cmerge **[NResonances];
	for(ires=0;ires<NResonances;ires++){
		MergeArray[ires]=new Cmerge *[NResonances];
		for(ichannelres=0;ichannelres<NResonances;ichannelres++){
			MergeArray[ires][ichannelres]=NULL;
		}
	}

	//now, use the stored decay information to create branchlists
	for(iter=resmap.begin();iter!=resmap.end();++iter){
		resinfo=iter->second;
		if(resinfo->decay){
			motherpid=iter->first;
			decayinfo=(decaymap.find(motherpid))->second; //decaymap[motherpid];
			bmax=0.0;
			for(ichannel=0; ichannel<resinfo->nchannels; ichannel++) {
				nbodies=decayinfo->Nparts[ichannel];
				bptr=new CbranchInfo();
				bptr->resinfo.clear();
				resinfo->branchlist.push_back(bptr);
				bptr->branching=decayinfo->branchratio[ichannel];
				//netq=-resinfo->charge;
				//netb=-resinfo->baryon;
				//nets=-resinfo->strange;

				nophotons=true;
				for(ibody=0; ibody<nbodies; ibody++) {
					pid=decayinfo->products[ichannel][ibody];
					if(pid==22)
						nophotons=false;
					bptr->resinfo.push_back(GetResInfoPtr(pid));
					//netq+=bptr->resinfo[ibody]->charge;
					//netb+=bptr->resinfo[ibody]->baryon;
					//nets+=bptr->resinfo[ibody]->strange;
				}
				bptr->L=decayinfo->d_L[ichannel];

				if(nophotons){
					ires1=bptr->resinfo[0]->ires;
					ires2=bptr->resinfo[1]->ires;
					if(ires1>ires2){
						iresflip=ires1; ires1=ires2; ires2=iresflip;
					}
					merge=MergeArray[ires1][ires2];
					if(merge==NULL){
						MergeArray[ires1][ires2]=new Cmerge(resinfo,bptr->branching, LDecay);
					}
					else{
						while(merge->next!=NULL){
							merge=merge->next;
						}
						merge->next=new Cmerge(resinfo,bptr->branching, LDecay);
					}
				}

				// switch places to make sure first branch has largest
				if(bptr->branching>bmax){
					bmax=bptr->branching;
					if(ichannel>0){
						firstbptr=resinfo->branchlist[0];
						resinfo->branchlist[0]=bptr;
						resinfo->branchlist[ichannel]=firstbptr;
					}
				}
			}  //out of channel loops
		}
	}
	for(iter=resmap.begin();iter!=resmap.end();++iter){
		resinfo=iter->second;
		motherpid=iter->first;
		decayinfo=decaymap[motherpid];
		delete decayinfo;
	}
	decaymap.clear();
}

void CresList::ReadSpectralFunctions(){
	CresMassMap::iterator rpos;
	CresInfo *resinfo;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->decay && !resinfo->SFcalculated){
			resinfo->ReadSpectralFunction();
			resinfo->NormalizeSF();
		}
		else
			resinfo->SFcalculated=true;
	}
}

void CresList::ReadResInfo_MSU(){
	Cmerge *merge;
	int mothercode,code,decay,NResonances;
	double mothermass,netm,bmax,spin;
	int ires,jres,ires1,ires2,iresflip,ichannel,nchannels,ibody,nbodies,LDecay;
	int netq,netb,nets;
	string name, filename;
	CresInfo *resinfoptr=NULL;
	CbranchInfo *bptr=NULL,*firstbptr=NULL;
	FILE *resinfofile;
	FILE * decayinfofile;
	char dummy[200],cname[200];
	filename=parmap->getS("RESONANCES_INFO_FILE",string("../resinfo/resonances_standardhadrons.dat"));
	resinfofile=fopen(filename.c_str(),"r");
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fscanf(resinfofile,"%d",&NResonances);
	fgets(dummy,200,resinfofile);
	MergeArray=new Cmerge **[NResonances];
	for(ires=0;ires<NResonances;ires++){
		MergeArray[ires]=new Cmerge *[NResonances];
		for(jres=0;jres<NResonances;jres++){
			MergeArray[ires][jres]=NULL;
		}
	}
	for(ires=0;ires<NResonances;ires++){
		resinfoptr=new CresInfo();
		fscanf(resinfofile,"%d %lf %d %d %d %lf %d %d %lf", &resinfoptr->pid,&resinfoptr->mass,&resinfoptr->charge,&resinfoptr->baryon, &resinfoptr->strange,&spin,&resinfoptr->G_Parity,&decay,&resinfoptr->width);
		resinfoptr->mass*=0.001;
		resinfoptr->minmass=resinfoptr->mass;
		resinfoptr->width*=0.001;
		resinfoptr->degen=lrint(2.0*spin+1);

		resinfoptr->q[0]=resinfoptr->baryon+resinfoptr->charge;
		resinfoptr->q[1]=2.0*resinfoptr->baryon+resinfoptr->strange-resinfoptr->charge;
		resinfoptr->q[2]=-resinfoptr->strange;
		
		fgets(cname,100,resinfofile);
		cname[int(strlen(cname))-1]='\0';
		resinfoptr->name=cname;
		resinfoptr->decay=bool(decay);
		resinfoptr->ires=ires;
		resinfoptr->branchlist.clear();
		resmap.insert(CresInfoPair(resinfoptr->pid,resinfoptr));
		massmap.insert(CresMassPair(resinfoptr->mass,resinfoptr));
	} 
	fclose(resinfofile);

	filename=parmap->getS("RESONANCES_DECAYS_FILE",string("../resinfo/decays_pdg_weak.dat"));
	decayinfofile=fopen(filename.c_str(),"r");
	while(fscanf(decayinfofile,"%d %lf",&mothercode,&mothermass) && !feof(decayinfofile)){
		mothermass*=0.001;
		fgets(dummy,200,decayinfofile);
		fscanf(decayinfofile,"%d %d",&mothercode,&nchannels);
		resinfoptr=GetResInfoPtr(mothercode);
		resinfoptr->minmass=1.0E10;
		bmax=0.0;
		for(ichannel=0;ichannel<nchannels;ichannel++){
			bptr=new CbranchInfo();
			bptr->resinfo.clear();
			resinfoptr->branchlist.push_back(bptr);
			fscanf(decayinfofile,"%d",&nbodies);
			netq=-resinfoptr->charge;
			netb=-resinfoptr->baryon;
			nets=-resinfoptr->strange;
			netm=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				fscanf(decayinfofile,"%d",&code);
				bptr->resinfo.push_back(GetResInfoPtr(code));
				netq+=bptr->resinfo[ibody]->charge;
				netb+=bptr->resinfo[ibody]->baryon;
				nets+=bptr->resinfo[ibody]->strange;
				netm+=bptr->resinfo[ibody]->mass;
			}
			//total charge and baryon number should be conserved, and shouldn't be larger than single strangeness
			if(netq!=0 || netb!=0 || abs(nets)>1){
				CLog::Info("Charge conservation failure while reading decay info,\nnetq="+to_string(netq)+", netb"+to_string(netb)+", nets="+to_string(nets)+"\n");
				CLog::Info("MOTHER (ichannel="+to_string(ichannel)+", nbodies="+to_string(nbodies)+":\n");
				resinfoptr->Print();
				CLog::Info("DAUGHTERS:\n");
				for(ibody=0;ibody<nbodies;ibody++)
					bptr->resinfo[ibody]->Print();
				if(netq!=0 || netb!=0)
					exit(1);
			}
			fscanf(decayinfofile,"%lf %d",&bptr->branching,&LDecay);
			
			//store two body decays only
			if(nbodies==2){
				ires1=bptr->resinfo[0]->ires;
				ires2=bptr->resinfo[1]->ires;
				if(ires1>ires2){
					iresflip=ires1; ires1=ires2; ires2=iresflip;
				}
				merge=MergeArray[ires1][ires2];
				if(merge==NULL){
					MergeArray[ires1][ires2]=new Cmerge(resinfoptr,bptr->branching, LDecay);
				}
				else{
					while(merge->next!=NULL){
						merge=merge->next;
					}
					merge->next=new Cmerge(resinfoptr,bptr->branching, LDecay);
				}
			}
			//if the total mass is smaller than the minimum required mass, replace it
			if(netm<resinfoptr->minmass){
				resinfoptr->minmass=netm;
				resinfoptr->bptr_minmass=bptr;
			}
			// switch places to make sure first branch has largest 
			if(bptr->branching>bmax){
				bmax=bptr->branching>bmax;
				if(ichannel>0){
					firstbptr=resinfoptr->branchlist[0];
					resinfoptr->branchlist[0]=bptr;
					resinfoptr->branchlist[ichannel]=firstbptr;
				}
			}
		}
	}
	
	fclose(decayinfofile);
}


