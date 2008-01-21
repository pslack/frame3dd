#include "array.h"
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>

static void *array_grow(array *a,unsigned index);

array array_create(size_t eachsize, unsigned cap){
	array a;
	a.eachsize = eachsize;
	a.num = 0;
	a.cap = cap;
	a.data = malloc(eachsize * cap);
	//fprintf(stderr,"Create array with capacity %d for objects of size %d\n",cap, eachsize);
	return a;
}

void array_destroy(array *a){
	if(a->data)free(a->data);
	a->data = NULL;
	a->num = 0;
	a->cap = 0;
}	

void *array_set(array *a, unsigned index, void *val){
	void *pos;

	//fprintf(stderr,"Adding to array at position %d\n",index);

	if(index >= a->cap){
		if(!array_grow(a, index)){
			return NULL;
		}
	}
	if(index >= a->num){
		//fprintf(stderr,"Incremented num from %d to %d\n",a->num,index+1);
		a->num = index + 1;
	}

	//fprintf(stderr,"a.num = %d, a.data = %p\n",a->num, a->data);
	//fprintf(stderr,"a.eachsize = %d\n",a->eachsize);

	pos = a->data + (a->eachsize * index);

	//fprintf(stderr,"pos = %p\n",pos);

	if(!memcpy(pos, val, a->eachsize)){
		return NULL;
	}

	return pos;
}

void *array_get(array *a, unsigned index){
	if(index >= a->num){
		return NULL;
	}
	return a->data + (a->eachsize * index);
}

void *array_grow(array *a,unsigned index){
	unsigned newcap;
	void *newdata;

	if(index < a->cap){
		return a->data;
	}
	if(index < a->cap * 2){
		newcap = a->cap * 2;
	}else{
		newcap = index + 1;
	}
	
	newdata = malloc(a->eachsize * newcap);
	if(!newdata){
		return NULL;
	}
	if(!memcpy(newdata,a->data,a->eachsize * a->num)){
		return NULL;
	}
	free(a->data);
	a->data = newdata;
	return a->data;
}

void *array_append(array *a,void *val){
	return array_set(a, a->num, val);
}

array array_copy(const array oa){
	array a;
	a = array_create(oa.eachsize, oa.num);
	//fprintf(stderr,"eachsize = %d in array_copy\n",oa.eachsize);
	a.data = malloc(a.eachsize * oa.num);
	if(!a.data){
		fprintf(stderr,"unable to allocate memory in array_copy\n");
		return a;
	}
	if(!memcpy(a.data,oa.data,oa.eachsize*oa.num)){
		fprintf(stderr,"unable to memcpy in array_copy\n");
		return a;
	}
	a.num = oa.num;
	return a;
}



