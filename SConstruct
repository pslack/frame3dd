#!/usr/bin/python SCons
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

version = '0.20100105'

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

vars = Variables(['options.cache'])

vars.Add(
	"CC"
	,"C compiler"
	,"gcc"
)
vars.Add(
	"CXX"
	,"C++ compiler"
	,"g++"
)

vars.Add(BoolVariable(
	"DEBUG"
	,"Whether to add debugger symbols to the binary"
	,1
))

vars.Add(BoolVariable(
	"HAVE_ITOA"
	,"Do you standard C libraries include the function 'itoa'?"
	,default_itoa
))

vars.Add(BoolVariable(
	"WITH_GCCVISIBILITY"
	,"Set true if you want to use the GCC 'visibility' feature."
	,True
))

vars.Add("INSTALL_PREFIX","Install location prefix (usually /usr or /usr/local)","/usr/local")
vars.Add("INSTALL_BIN","Install location for binaries (Linux)","$INSTALL_PREFIX/bin")
vars.Add("INSTALL_LIB","Install location for libraries (Linux)","$INSTALL_PREFIX/lib")
vars.Add("INSTALL_INCLUDE","Install location for header files (Linux)","$INSTALL_PREFIX/include")
vars.Add("INSTALL_DATA","Install location for general data files (Linux)","$INSTALL_PREFIX/share");
vars.Add("INSTALL_DOC","Install location for documentation files (Linux)","$INSTALL_DATA/doc/frame3dd");
vars.Add("INSTALL_FRAMEDATA","Install location FRAME's data files (Linux)","$INSTALL_DATA/frame3dd");
vars.Add("INSTALL_ROOT","Install root (for building RPMs etc)","")

# NSIS TARGET FILENAME

vars.Add(
	'WIN_INSTALLER_NAME'
	,'Windows Installer name'
	,'frame3dd-%s.exe' % version
)

vars.Update(env)
vars.Save('options.cache',env)
Help(vars.GenerateHelpText(env))

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
		env1 = context.env.Clone()
		env1.ParseConfig(cppunit_config_command + ["--libs","--cflags"])
		context.env["CPPUNIT_LIBS"] = env1.get("LIBS")
		context.env["CPPUNIT_LIBPATH"] = env1.get("LIBPATH")
		context.env["CPPUNIT_CPPPATH"] = env1.get("CPPPATH")
	context.Result(res)
	return res

#--------
# SOQT

soqt_test_code = """
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoCube.h>

int main(int argc, char **argv){
	QWidget * mainwin = SoQt::init(argv[0]);
	SoCube * cube = new SoCube;
   return 0;
}
"""
	
def CheckSoQt(context):
	context.Message("Checking for SoQt... ")
	if not context.env.get("SOQT_LIBS"):		
		context.Result(False)
		return False
	old_env = context.env.Clone()
	context.env.Append(
		CPPPATH = env.get('SOQT_CPPPATH')
		, LIBS = env.get('SOQT_LIBS')
		, LIBPATH = env.get('SOQT_LIBPATH')
		, CPPDEFINES = env.get('SOQT_CPPDEFINES')
	)
	res = context.TryLink(soqt_test_code,".cpp")
	context.Result(res)
	for i in ['LIBS','CPPPATH','LIBPATH','CPPDEFINES']:
		if old_env.get(i) is not None:
			context.env[i] = old_env[i]
		else:
			del context.env[i]
	return res

#-------------------------------------------------------------------------------
conf = Configure(env
	, custom_tests = { 
		'CheckGcc' : CheckGcc
		, 'CheckGccVisibility' : CheckGccVisibility
		, 'CheckCppUnitConfig' : CheckCppUnitConfig
		, 'CheckSoQt' : CheckSoQt
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

env['HAVE_SOQT'] = conf.CheckSoQt()

#-------------
# documentation

env.Append(SUBST_DICT= {
	'@VERSION@':version
	,'@CHANGELOG@':file("ChangeLog.txt").read()
})

env.SConscript('doc/SConscript',['env'])

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
	exA.3dd  exC.3dd  exE.3dd  exG.3dd  exI.3dd
	exB.3dd  exD.3dd  exF.3dd  exH.3dd  
""")

exampledir=Dir(env.subst("$INSTALL_ROOT$INSTALL_FRAMEDATA/examples"))
libdir=Dir(env.subst("$INSTALL_ROOT$INSTALL_LIB"))
bindir=Dir(env.subst("$INSTALL_ROOT$INSTALL_BIN"))
incdir=Dir(env.subst("$INSTALL_ROOT$INSTALL_INCLUDE"))
env.Install(exampledir,['examples/%s'%e for e in examples])
env.Install(exampledir,['test/truss.arc','test/bent-cantilever.arc'])

#------------
# install documentation

docdir = Dir(env.subst("$INSTALL_ROOT$INSTALL_DOC"))
env.Install(docdir,['doc/user-manual.html','doc/version.html'])
docimgdir = Dir(env.subst("$INSTALL_ROOT$INSTALL_DOC/img"))
env.Install(docimgdir,Glob("doc/img/*.jpg"))
env.Install(docimgdir,Glob("doc/img/*.png"))

#------------
# register the 'install' target

datadir=Dir(env.subst("$INSTALL_ROOT$INSTALL_FRAMEDATA"))
env['installdirs']+=[datadir, libdir, bindir, incdir, docdir, docimgdir]

env.Alias("install",env['installdirs'])

#------------
# create the RPM .spec file

specfile = env.SubstInFile('frame3dd.spec.in')
Depends(specfile,"ChangeLog.txt")

#------------
# create distribution zip-file

env['DISTTAR_FORMAT']='bz2'
env.Append(
	DISTTAR_EXCLUDEEXTS=['.o','.os','.so','.a','.dll','.cc','.cache','.pyc'
		,'.cvsignore','.dblite','.log', '.gz', '.bz2', '.zip', '.patch', '.mm'
		,'.out','.tmp','.swp','.db']
	,DISTTAR_EXCLUDEDIRS=['CVS','.svn','.sconf_temp', 'dist','build'
		,'development','buildings','gui','images','tests']
)

tar = env.DistTar("dist/frame3dd-"+version
        , [env.Dir('#')]
)

env.Depends(tar,["frame3dd.spec","doc/version.html"])

if platform.system()=="Windows":
	env.Append(NSISDEFINES={
		'OUTFILE':"#dist/"+env['WIN_INSTALLER_NAME']
		,'VERSION':version
	})
	
	installer = env.Installer('installer.nsi')
	Depends(installer,env['PROGS'])
	env.Alias('installer',installer)

#-------------
# debian.tar.gz for Debian packaging (Ubuntu,...)

import glob
deb_files = glob.glob('debian/*.install')
deb_files += glob.glob('debian/*.docs')
deb_files += glob.glob('debian/*.dirs')
deb_files += glob.glob('debian/*.man')
deb_files += glob.glob('debian/*.manpages')
deb_files += ['debian/%s' % s for s in ['rules','control','changelog','compat','copyright','dirs']]

deb_tar = env.Tar(
	'dist/debian.tar.gz'
	,deb_files
	,TARFLAGS = ['cz']
)

#-------

env.Default(env['PROGS'] + specfile)

# vim: set syntax=python :
