import socket
hostname = socket.gethostname()
if 'antares' in hostname:
	import matplotlib 
	matplotlib.use('Agg')