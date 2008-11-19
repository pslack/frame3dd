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
			env['SOQT_CPPDEFINES'] = ['SOQT_DLL']
			
			x=_winreg.ConnectRegistry(None,_winreg.HKEY_LOCAL_MACHINE)
			y= _winreg.OpenKey(x,r"SOFTWARE\coin3d")
			LIB,t = _winreg.QueryValueEx(y,"INSTALL_LIB")
			BIN,t = _winreg.QueryValueEx(y,"INSTALL_BIN")
			INCLUDE,t = _winreg.QueryValueEx(y,"INSTALL_INCLUDE")
			
			env.AppendUnique(
				SOQT_CPPPATH = [INCLUDE]
				,SOQT_LIBPATH = [LIB]
				,SOQT_LIBS = ['Coin']
				,SOQT_CPPDEFINES = ['COIN_DLL']
			)
			
			env.AppendUnique(
				QT_PREFIX = r'c:/Qt/4.3.3'
				,SOQT_CPPPATH = ["$QT_PREFIX/include","$QT_PREFIX/include/Qt"]
				,SOQT_LIBPATH = ["$QT_PREFIX/lib"]
			)
			
		else:
			cmd = ['soqt-config','--cppflags','--ldflags','--libs']
			env1 = env.Clone()
			env1.ParseConfig(cmd)
			env['SOQT_CPPPATH'] = env1.get('CPPPATH')
			env['SOQT_LIBPATH'] = env1.get('LIBPATH')
			env['SOQT_LIBS'] = env1.get('LIBS')
			env['SOQT_CPPDEFINES'] = env1.get('CPPDEFINES')

		print "SOQT_LIBS =",env.get('SOQT_LIBS')
		print "SOQT_LIBPATH =",env.get('SOQT_LIBPATH')
		print "SOQT_CPPPATH =",env.get('SOQT_CPPPATH')
		print "SOQT_CPPDEFINES =",env.get('SOQT_CPPDEFINES')

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

