import os, platform
from SCons.Script import *

def generate(env):
	"""
	Detect SOQT settings and add them to the environment.
	"""
	try:
		if platform.system()=="Windows":
			import _winreg
			x=_winreg.ConnectRegistry(None,_winreg.HKEY_LOCAL_MACHINE)
			y= _winreg.OpenKey(x,r"SOFTWARE\soqt")
			LIB,t = _winreg.QueryValueEx(y,"INSTALL_LIB")
			BIN,t = _winreg.QueryValueEx(y,"INSTALL_BIN")
			INCLUDE,t = _winreg.QueryValueEx(y,"INSTALL_INCLUDE")

			env['SOQT_CPPPATH'] = [INCLUDE]
			env['SOQT_LIBPATH'] = [LIB]
			env['SOQT_LIBS'] = ['SoQt']
		else:
			cmd = ['soqt-config','--cppflags','--ldflags','--libs']
			env1 = env.Copy()
			env1.ParseConfig(cmd)
			env['SOQT_CPPPATH'] = env1.get('CPPPATH')
			env['SOQT_LIBPATH'] = env1.get('LIBPATH')
			env['SOQT_LIBS'] = env1.get('LIBS')

		print "SOQT_LIBS =",env.get('SOQT_LIBS')
		print "SOQT_LIBPATH =",env.get('SOQT_LIBPATH')
		print "SOQT_CPPPATH =",env.get('SOQT_CPPPATH')

	except:
		print "FAILED TO SET UP SOQT"
		pass

def exists(env):
	"""
	Make sure this tool exists.
	"""
	if platform.system()=="Windows":
		try:
			import _winreg
			x=_winreg.ConnectRegistry(None,_winreg.HKEY_LOCAL_MACHINE)
			y= _winreg.OpenKey(x,r"SOFTWARE\soqt")
			INCLUDE,t = _winreg.QueryValueEx(y,'INSTALL_INCLUDE')
			return True
		except:
			return False
	else:
		if env.Which('soqt-config'):
			return True
		return False

