#!/usr/bin/env python
"""Created on 28-June-2017
by wojdyla_j, modified by basu_s
"""
import os
import pwd
from subprocess import Popen, PIPE, STDOUT
import shlex
import logging

def demote(user_uid, user_gid):
	def result():
		os.setgid(user_gid)
		os.setuid(user_uid)
	return result

def run_command(logger_name, directory_name, user, command, logname):
	"""Function to run command of choice with user credentials.
	Works only when original script is run as root"""
	logger = logging.getLogger('{}'.format(logger_name))
	tpath = directory_name
	try:
		pw_record = pwd.getpwnam(user)
		logger.info('pw_record:{}'.format(pw_record))
	except (OSError, TypeError, Exception) as e:
		logger.info('switch to local mode'.format(e))

	user_name      = pw_record.pw_name
	user_home_dir  = pw_record.pw_dir
	user_uid       = pw_record.pw_uid
	user_gid       = pw_record.pw_gid
	env = os.environ.copy()
	env[ 'HOME'     ]  = user_home_dir
	env[ 'LOGNAME'  ]  = user_name
	env[ 'PWD'      ]  = tpath
	env[ 'USER'     ]  = user_name
	os.chdir(tpath)
	logger.info('Running with command: {}'.format(command))
	logger.info('Directory name: {}'.format(directory_name))
	logger.info('User_uid:{}, user_gid: {}'.format(user_uid, user_gid))
	logger.info('Regaining root privilages to run command')
	os.seteuid(0)
	os.setegid(0)
	logger.info('gid:{}, uid:{}'.format(os.getegid(),os.geteuid()))
	logger.info('user_uid: {}, user_gui: {}, user_name: {}, user_home_dir: {}'.format(user_uid, user_gid, user_name, user_home_dir))
	output = Popen(shlex.split(command), preexec_fn=demote(user_uid, user_gid), env=env, stdout=PIPE, stderr=STDOUT).communicate()[0]
	os.setegid(pw_record.pw_gid)
	os.seteuid(pw_record.pw_uid)
	logger.info('Returning to user privilages')
	with open(os.path.join(directory_name, logname), 'a+') as log:
		log.write(output)
	return output
