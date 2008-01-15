# CSconstruct for FRAME
# vim: set syntax=python:

version = '20071220'

env = Environment(
	tools=['default','disttar','substinfile']
	,toolpath=['scons']
)

opts = Options()

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

opts.Add("INSTALL_PREFIX","Install location prefix (usually /usr or /usr/local)","/usr/local")
opts.Add("INSTALL_BIN","Install location for binaries","$INSTALL_PREFIX/bin")
opts.Add("INSTALL_DATA","Install location for general data files","$INSTALL_PREFIX/share");
opts.Add("INSTALL_FRAMEDATA","Install location FRAME's data files","$INSTALL_DATA/frame");
opts.Add("INSTALL_ROOT","Install root (for building RPMs etc)","")

opts.Update(env)

env['CCFLAGS']=['-O']

if env.get('DEBUG'):
	env.Append(CCFLAGS=['-g'])

#-------------
# build the program

env['installdirs'] = []

env.BuildDir('build','src')
env.SConscript('build/SConscript',['env'])

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

env.Append(SUBST_DICT= {
	'@VERSION@':version
})

env.SubstInFile('frame.spec.in')

#------------
# create distribution zip-file

env['DISTTAR_FORMAT']='bz2'
env.Append(
	DISTTAR_EXCLUDEEXTS=['.o','.os','.so','.a','.dll','.cc','.cache','.pyc'
		,'.cvsignore','.dblite','.log', '.gz', '.bz2', '.zip', '.patch']
	,DISTTAR_EXCLUDEDIRS=['CVS','.svn','.sconf_temp', 'dist','build'
		,'development','buildings','gui','images','tests']
)

tar = env.DistTar("dist/frame-"+version
        , [env.Dir('#')]
)

env.Depends(tar,"frame.spec")


# vim: set syntax=python:
