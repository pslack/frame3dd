"""	FRAME: Static and dynamic structural analysis of 2D & 3D frames and trusses
	Copyright (C) 1992-2007  Henri P. Gavin

	This program is free software; you may redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

version = '0.20080911'

import platform
deftools = ['default']
if platform.system()=="Windows":
	deftools = ['mingw']
	default_itoa=1
else:
	default_itoa=0

env = Environment(
	tools=deftools + ['disttar','substinfile','soqt','nsis']
	,toolpath=['scons']
)

opts = Options(['options.cache'])

opts.Add(
	"CC"
	,"C compiler"
	,"gcc"
)

opts.Add(BoolOption(
	"DEBUG"
	,"Whether to add debugger symbols to the binary"
	,1
))

opts.Add(BoolOption(
	"HAVE_ITOA"
	,"Do you standard C libraries include the function 'itoa'?"
	,default_itoa
))

opts.Add(BoolOption(
	"WITH_GCCVISIBILITY"
	,"Set true if you want to use the GCC 'visibility' feature."
	,True
))

opts.Add("INSTALL_PREFIX","Install location prefix (usually /usr or /usr/local)","/usr/local")
opts.Add("INSTALL_BIN","Install location for binaries","$INSTALL_PREFIX/bin")
opts.Add("INSTALL_DATA","Install location for general data files","$INSTALL_PREFIX/share");
opts.Add("INSTALL_FRAMEDATA","Install location FRAME's data files","$INSTALL_DATA/frame");
opts.Add("INSTALL_ROOT","Install root (for building RPMs etc)","")

# NSIS TARGET FILENAME

opts.Add(
	'WIN_INSTALLER_NAME'
	,'Windows Installer name'
	,'frame3dd-%s.exe' % version
)

opts.Update(env)
opts.Save('options.cache',env)
Help(opts.GenerateHelpText(env))

env['CCFLAGS']=['-O', '-Wall']
env['VERSION'] = version

if env.get('DEBUG'):
	print "DEBUGGING TURNED ON"
	env.Append(CCFLAGS=['-g'])

#====================
# CONFIGURATION TESTS

#----------------
# GCC

gcc_test_text = """
#ifndef __GNUC__
# error "Not using GCC"
#endif

int main(void){
	return __GNUC__;
}
"""

def CheckGcc(context):
	context.Message("Checking for GCC... ")
	is_ok = context.TryCompile(gcc_test_text,".c")
	context.Result(is_ok)
	return is_ok

#----------------
# GCC VISIBILITY feature

gccvisibility_test_text = """
#if __GNUC__ < 4
# error "Require GCC version 4 or newer"
#endif

__attribute__ ((visibility("default"))) int x;

int main(void){
	extern int x;
	x = 4;
}
"""

def CheckGccVisibility(context):
	context.Message("Checking for GCC 'visibility' capability... ")
	if not context.env.has_key('WITH_GCCVISIBILITY') or not env['WITH_GCCVISIBILITY']:
		context.Result("disabled")
		return 0
	is_ok = context.TryCompile(gccvisibility_test_text,".c")
	context.Result(is_ok)
	return is_ok

#--------
# CPPUNIT

if platform.system()=="Windows":
	cppunit_config_command = ["cppunit-config"]
else:	
	cppunit_config_command = ["cppunit-config"]

def CheckCppUnitConfig(context):
	res = 0
	context.Message("Checking for cppunit-config... ")
	if context.env.WhereIs(cppunit_config_command[0]):
		res = 1	
		context.env.ParseConfig(cppunit_config_command + ["--libs","--cflags"])

	context.Result(res)
	return res

def CheckCppUnit(context):
	envtemp = conf.env.Copy()
	conf.env.Parseconfig
	
		
	context.Result(1)
	return 1

conf = Configure(env
	, custom_tests = { 
		'CheckGcc' : CheckGcc
		, 'CheckGccVisibility' : CheckGccVisibility
		, 'CheckCppUnitConfig' : CheckCppUnitConfig
	} 
)

if conf.CheckGcc():
	conf.env['HAVE_GCC']=True;
	if env['WITH_GCCVISIBILITY'] and conf.CheckGccVisibility():
		conf.env['HAVE_GCCVISIBILITY']=True;
		conf.env.Append(CCFLAGS=['-fvisibility=hidden'])
		conf.env.Append(CPPDEFINES=['HAVE_GCCVISIBILITY'])
	conf.env.Append(CCFLAGS=['-Wall'])

if conf.CheckCppUnitConfig():
	conf.env['HAVE_CPPUNIT']=True;

#-------------
# build the program

env['installdirs'] = []
env['PROGS'] = []

env.BuildDir('build','src')
env.SConscript('build/SConscript',['env'])

#------------
# test suite

env.BuildDir('build/test','test')
env.SConscript('build/test/SConscript',['env'])

#------------
# install example files

examples = Split("""
	exA.frm  exC.frm  exE.frm  exG.frm  exI.frm
	exB.frm  exD.frm  exF.frm  exH.frm  
""")

datadir=Dir(env.subst("$INSTALL_ROOT$INSTALL_FRAMEDATA"))
env.Install(datadir,['examples/%s'%e for e in examples])
env['installdirs']+=[datadir]

env.Alias("install",env['installdirs'])

#------------
# create the RPM .spec file
#---------------------------------------------------
# CPPUNIT

if platform.system()=="Windows":
	cppunit_config_command = ["cppunit-config"]
else:	
	cppunit_config_command = ["cppunit-config"]

def CheckCppUnitConfig(context):
	res = 0
	context.Message("Checking for cppunit-config... ")
	if context.env.WhereIs(cppunit_config_command[0]):
		res = 1	
		context.env.ParseConfig(cppunit_config_command + ["--libs","--cflags"])

	context.Result(res)
	return res

def CheckCppUnit(context):
	envtemp = conf.env.Copy()
	conf.env.Parseconfig
	
		
	context.Result(1)
	return 1
env.Append(SUBST_DICT= {
	'@VERSION@':version
})

env.SubstInFile('frame3dd.spec.in')

#------------
# create distribution zip-file

env['DISTTAR_FORMAT']='bz2'
env.Append(
	DISTTAR_EXCLUDEEXTS=['.o','.os','.so','.a','.dll','.cc','.cache','.pyc'
		,'.cvsignore','.dblite','.log', '.gz', '.bz2', '.zip', '.patch', '.mm'
		,'.out','.tmp','.swp']
	,DISTTAR_EXCLUDEDIRS=['CVS','.svn','.sconf_temp', 'dist','build'
		,'development','buildings','gui','images','tests']
)

tar = env.DistTar("dist/frame3dd-"+version
        , [env.Dir('#')]
)

env.Depends(tar,"frame3dd.spec")

if platform.system()=="Windows":
	env.Append(NSISDEFINES={
		'OUTFILE':"#dist/"+env['WIN_INSTALLER_NAME']
	})
	
	installer = env.Installer('installer.nsi')
	Depends(installer,env['PROGS'])
	env.Alias('installer',installer)


#-------

env.Default(env['PROGS'])


# vim: set syntax=python:
