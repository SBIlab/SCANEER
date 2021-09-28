# Definition of Path
import os, sys
#lib_path = os.path.dirname(__file__)
lib_path = os.path.dirname(os.path.abspath(__file__))

# Get Executable File Path of Platform-dependent Program
def get_executable_arch(exe_path):
	# Win32 Platform
	if sys.platform == 'win32':
		arch_name = 'win32'
		exe_path += '.exe'
		exe_path = exe_path.replace('/', '\\')

	# Linux Platform
	elif 'linux' in sys.platform:
		arch_name = 'linux'

	# POSIX Compatible Platform: Non-implemented
	else:
		print "ERROR: Not supported Platform!"
		return None

	# Return Executable Path
	return exe_path.replace('<arch>', arch_name)
