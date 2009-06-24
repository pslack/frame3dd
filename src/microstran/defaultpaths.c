#define MSTRANP_BUILD
#include "defaultpaths.h"

#ifdef WIN32
#include <windows.h>
#endif

#define FRAME3DD_REGKEY_ROOT HKEY_LOCAL_MACHINE
#define FRAME3DD_REGKEY_MAIN "SOFTWARE\\FRAME3DD"

const char *get_default_data_path(){
#ifdef WIN32
# define MAXLEN 3000
	HKEY key;
	DWORD datatype, len = MAXLEN;
	long res;
	static char value[MAXLEN];
	res = RegOpenKeyEx(FRAME3DD_REGKEY_ROOT, FRAME3DD_REGKEY_MAIN, 0L, KEY_QUERY_VALUE, &key);
	if(res==ERROR_SUCCESS){
		res = RegQueryValueEx(key, "Install_Dir", NULL, &datatype, value, &len);
		if(res==ERROR_SUCCESS){
			RegCloseKey(key);
			return value;
		}
	}
#endif
	return FRAME3DD_DEFAULT_DATA_DIR;
}
