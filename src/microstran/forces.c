/*	FRAME3DD: Static and dynamic structural analysis of 2D & 3D frames and trusses
	Copyright (C) 2007-2008 John Pye

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define MSTRANP_BUILD

#include "forces.h"
#include "new.h"
#include "array.h"
#include "error.h"

#include <stdio.h>
#include <string.h>

modelforces *modelforces_create(){
	modelforces *mf;
	mf = NEW(modelforces);
	if(!mf)abort();
	strcpy(mf->filename,"");
	mf->cases = ARRAY_CREATE(caseforces,10);
	return mf;
}

void modelforces_destroy(modelforces *mf){
	unsigned i;
	for(i=0;i<ARRAY_NUM(mf->cases);++i){
		caseforces_destroy(array_get(&(mf->cases),i));
	}
	array_destroy(&(mf->cases));
	free(mf);
}

caseforces *caseforces_create(unsigned id, const char *name){
	caseforces *cf;
	cf = NEW(caseforces);
	if(!cf)abort();
	cf->members = ARRAY_CREATE(memberforce,10);
	cf->id = id;
	strncpy(cf->name,name,MAXCASENAME);
	cf->name[MAXCASENAME-1]='\0';
	return cf;
}

void caseforces_destroy(caseforces *cf){
	array_destroy(&(cf->members));
	free(cf);
}

membernodeforce membernodeforce_create(unsigned nodeid, vec3 F, vec3 M){
	membernodeforce mnf;
	mnf.node = nodeid;
	mnf.F = F;
	mnf.M = M;
	return mnf;
}

cbool caseforces_add_member(caseforces *cf, unsigned memberid, membernodeforce fromnode, membernodeforce tonode){
	memberforce mf;
	//fprintf(stderr,"ADDING MEMBER %d to case %d\n",memberid,cf->id);
	mf.member = memberid;
	mf.fromnode = fromnode;
	mf.tonode = tonode;
	array_append(&(cf->members),&mf);
	return 1;
}

cbool modelforces_add_case(modelforces *mf, caseforces *cf){
	//fprintf(stderr,"Model forces: contains %d cases\n",ARRAY_NUM(mf->cases));
	array_append(&(mf->cases),cf);
	return 1;
}

/* QUERY FUNCTIONS */

caseforces *modelforces_find_case(modelforces *mf, unsigned caseid){
	unsigned i;
	//fprintf(stderr,"Model displacement contains %d cases\n",ARRAY_NUM(md->cases));
	caseforces *cf;
	for(i=0; i<ARRAY_NUM(mf->cases); ++i){
		cf = array_get(&(mf->cases), i);
		if(cf->id == caseid){
			return cf;
		}
	}
	return NULL;
}	

memberforce *caseforces_find_member(caseforces *cf, const unsigned memberid){
	unsigned i;
	memberforce *mf;
	for(i=0; i<ARRAY_NUM(cf->members); ++i){
		mf = array_get(&(cf->members), i);
		if(mf->member == memberid){
			return mf;
		}
	}
	return NULL;
}

caseforces *modelforces_get_case(modelforces *mf, unsigned index){
	return array_get(&(mf->cases), index);
}

