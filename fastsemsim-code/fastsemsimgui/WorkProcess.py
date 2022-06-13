# -*- coding: iso-8859-1 -*-
'''
Copyright 2011 Marco Mina. All rights reserved.

This file is part of fastSemSim

fastSemSim is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fastSemSim is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with fastsemsim.  If not, see <http://www.gnu.org/licenses/>.
'''

'''
WorkProcess is the background process that handles data and evaluates the semantic similarity.
'''
import multiprocessing
import time 
import copy
import os
# import fastsemsim.SemSim import SemSimMeasures
import fastsemsim.SemSim
from fastsemsim.SemSim import ObjSemSim
from fastsemsim.SemSim import SemSimUtils
from fastsemsim.Ontology import GeneOntology
from fastsemsim.Ontology import AnnotationCorpus
import gzip

DEBUG_LEVEL = 0 # NOTE 0: no debug output. 1:level 1 debug output. 2: verbose debug output. In any case it does not affect the performance of the calculation step

# WorkProcess is designed as a state machine. Here are the possible status
STATUS_BASE = 0
STATUS_INIT = STATUS_BASE + 0 # performing reset.
STATUS_WAIT = STATUS_BASE + 1 # waiting for requests.
STATUS_RUN = STATUS_BASE + 2 # computation in progress.
STATUS_SET = STATUS_BASE + 3 # Loading data from the gui process.
STATUS_PAUSE = STATUS_BASE + 4 # Pause.
STATUS_SEND = STATUS_BASE + 5 # Sending data to gui process.
STATUS_KILL = STATUS_BASE + 6 # Killing process

# Commands handled by WorkProcess
CMD_BASE = 100
CMD_NONE = CMD_BASE
CMD_START = CMD_BASE + 1 # Start the computation
CMD_STOP = CMD_BASE + 2 # Stop the computation. Semantic similarity scores are cleared, but all the other data and settings are mantained
CMD_PAUSE = CMD_BASE + 3 # Pause the computation.
CMD_RESET = CMD_BASE + 4 # Reinizialize the class. Clear all the variables!
CMD_STATUS = CMD_BASE + 5 # Return status information
CMD_KILL =  CMD_BASE + 6 # Kill the process
CMD_OUTPUT = CMD_BASE + 10 # Output data
CMD_GET = CMD_BASE + 7 # Get
CMD_SET = CMD_BASE + 11 # Set


# Codes used in incoming set requests
CMD_LOAD_AC = CMD_SET + 0 # Load an Annotation Corpus
CMD_LOAD_GO = CMD_SET + 1 # Load a Gene Ontology
CMD_SET_QUERY = CMD_SET + 2 # Load a Query
CMD_SET_SS = CMD_SET + 3 # Initialize Semantic Similarity Engine
CMD_SET_OUTPUT = CMD_SET + 4 # Load Output parameters

# Codes used in incoming get requests
CMD_GET_AC = CMD_GET + 1
CMD_GET_AC_OBJECTS = CMD_GET_AC + 0
CMD_GET_AC_TERMS = CMD_GET_AC + 1
CMD_GET_AC_OBJECTS_NUMBER = CMD_GET_AC + 2
CMD_GET_AC_TERMS_NUMBER = CMD_GET_AC + 3

CMD_GET_PARAMS = CMD_GET + 100
CMD_GET_PARAMS_AC = CMD_GET_PARAMS + 0
CMD_GET_PARAMS_GO = CMD_GET_PARAMS + 1
CMD_GET_PARAMS_SS = CMD_GET_PARAMS + 2
CMD_GET_PARAMS_QUERY = CMD_GET_PARAMS + 3
CMD_GET_PARAMS_OUTPUT = CMD_GET_PARAMS + 4


# Codes used in outgoing messages
ANSWER_BASE = 100
ANSWER_IGNORED =  ANSWER_BASE + 1 # Previous message has been ignored 
ANSWER_PROCESSED =  ANSWER_BASE + 2 # Previous message has been accepted and completely processed
ANSWER_PROCESSING =  ANSWER_BASE + 3 # Previous message has been accepted but is still being process

# Results
RESULT_OK = True
RESULT_BAD = False

# Constants
LOAD_AC_BASE = 200
LOAD_AC_END = LOAD_AC_BASE + 1
LOAD_AC_STATUS = LOAD_AC_BASE + 2

LOAD_GO_BASE = 300 
LOAD_GO_END = LOAD_GO_BASE + 1
LOAD_GO_STATUS = LOAD_GO_BASE + 2

OUTPUT_TO_GUI = 0
OUTPUT_TO_FILE = 1

QUERY_FROM_GUI = 0
QUERY_FROM_AC = 2
QUERY_FROM_FILE = 1

QUERY_PAIRS = 0
QUERY_LIST = 1















	#-#-#-#-#-#-#-#-#-#-#-#
	# Class WorkProcess   #
	#-#-#-#-#-#-#-#-#-#-#-#

class WorkProcess(multiprocessing.Process):

	# some static settings
	print_output = False
	MAX_BUFFER_SIZE = 1000
	CURRENT_BUFFER_SIZE = MAX_BUFFER_SIZE
	CHECK_INTERVAL = 500

#### __init__ receives in input the 4 channels used for bidirectional communication with gui process
# WorkProcess communicates through 4 channels: 2 ingoing and 2 outgoing. In each direction a pipe and a queue are available.#
	def __init__(self, gui2ssprocess_queue, ssprocess2gui_queue, gui2ssprocess_pipe, ssprocess2gui_pipe):
		multiprocessing.Process.__init__(self)
		if DEBUG_LEVEL>0:
			print "WorkProcess: init()"
		self.output_queue = ssprocess2gui_queue
		self.output_pipe = ssprocess2gui_pipe
		self.input_queue = gui2ssprocess_queue
		self.input_pipe = gui2ssprocess_pipe
		self.status = STATUS_INIT

#### NOTE All the variables should be defined here!!!
	def reset(self):
		if DEBUG_LEVEL>0:
			print "WorkProcess: reset()"
# component ok flags. True if the corresponding component is correctly and completely configured and data are loaded.
		self.ok_ss = False
		self.ok_ac = False
		self.ok_go = False
		self.ok_query = False
		self.ok_output = False
# parameters ok flags. True if the parameters of the component are correctly configured. Data are not necessarily loaded though.
		self.ok_params_ss = False
		self.ok_params_ac = False
		self.ok_params_go = False
		self.ok_params_query = False
		self.ok_params_output = False
# GO parameters and datastructures
		self.go = None # Gene Ontology
		self.go_filename = None
		self.go_parameters = {}
# AC parameters and datastructures
		self.ac = None # Annotation Corpus
		self.ac_filename = None
		self.ac_filetype = None
		self.ac_filetypeparams = {}
# Query parameters and datastructures
		self.query = None # Query
		#### DANGER Fill missing variables
# SS parameters and datastructures
		self.ss = None # ss object implementing the desired semantic similarity measure
		self.ss_name = None
		self.ss_ontology = None
		self.ss_params = None
# Output parameters and datastructures
		self.output = None # 
		self.output_to = None
		self.output_file = None
		self.output_buffer = None
		#### DANGER Fill missing variables

		self.status = STATUS_WAIT
#

#-#-#-#-#-#-#-#-#-#-#-#	#-#-#-#-#-#-#-#-#-#-#-#	#-#-#-#-#-#-#-#-#-#-#-#	#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#	#-#-#-#-#-#-#-#-#-#-#-#	#-#-#-#-#-#-#-#-#-#-#-#	#-#-#-#-#-#-#-#-#-#-#













	#-#-#-#-#-#-#-#-#-#-#-#
	# main loop function  #
	#-#-#-#-#-#-#-#-#-#-#-#

	def run(self):
		if DEBUG_LEVEL>0:
			print "WorkProcess: run(). (Should appear once)"
		while self.control():
			pass
		if DEBUG_LEVEL>0:
			print "WorkProcess: end()"
		return 0
#





	def control(self):
		if self.status == STATUS_RUN:
			if DEBUG_LEVEL > 0:
				print "WorkProcess: control(). STATUS_RUN"
			try:
				if not self.input_queue.empty():
					data = self.input_queue.get(False)
					if data[0] == CMD_KILL: # kill the process
						self._kill()
						return False
					if data[0] == CMD_STOP:
						self._stop()
					#elif data[0] == CMD_PAUSE:
						#self.pause()
					elif data[0] == CMD_RESET: # reinitialize
						self.reset()
					elif data[0] == CMD_STATUS:
						self._status(data[1:len(data)])
						self._calculate()
					elif data[0] == CMD_START:
						self._calculate()
					else:
						print "Unhandled message in RUN status: " + str(data[0])
						self._calculate()
				else:
					self._calculate()
			except Exception:
				pass
		elif self.status == STATUS_INIT:
			self.reset()
		elif self.status == STATUS_WAIT:
			data = self.input_queue.get()
			if data[0] == CMD_KILL:
				self._kill()
				return False
			if data[0] == CMD_STOP: #### or data[0] == CMD_PAUSE: # No effect here
				pass
			elif data[0] == CMD_RESET:
				self.status = STATUS_INIT
			elif data[0] == CMD_SET:
				self.status = STATUS_SET
				self._set(data[1:len(data)])
			elif data[0] == CMD_GET:
				self._get(data[1:len(data)])
			elif data[0] == CMD_START:
				self._start(data[1:len(data)])
			elif data[0] == CMD_STATUS:
				self._status(data[1:len(data)])
			elif data[0] == CMD_NONE: # useless statement
				pass
			else: # any other message does not have any effect
				print "Unknown request. Ignoring request..."
				pass
		elif self.status == STATUS_SET:
			if not self.input_queue.empty():
				data = self.input_queue.get(False)
				if data[0] == CMD_STATUS:
					self._status(data[1:len(data)])
				elif data[0] == CMD_KILL: # kill the process
					self._kill()
					return False
				else:
					print "\"Set\" status. Ignoring request..."
#### Start To Remove 
		# in PAUSE status only STOP, PAUSE, RESET, STATUS and START commands are accepted. The other are ignored.
		# get is casted in blocking mode since no computation has to be performed
		# If no messages have been received the computation continues
		#elif self.status == STATUS_PAUSE:
			##print "STATUS: Pause"
			#data = self.gui2ssprocess_queue.get()
			#if data[0] == CMD_KILL:
				#return False
			#if data[0] == CMD_STOP:
				#self.stop()
			#elif data[0] == CMD_RESET:
				#self.status = STATUS_INIT
			#elif data[0] == CMD_STATUS:
				#self._status()
			#elif data[0] == CMD_START:
				#self._restart()
			#elif data[0] == CMD_PAUSE:
				#pass
			#else:
				#print "Unhandled message in PAUSE status: " + str(data[0])
				## send ignored message?
				#pass
#### End To Remove
		else:
			if not self.input_queue.empty():
				data = self.input_queue.get(False)
				if data[0] == CMD_STATUS:
					self._status(data[1:len(data)])
				elif data[0] == CMD_KILL: # kill the process
					self._kill()
					return False
			print "Invalid status. Reinitializing..."
			self.status = STATUS_INIT
		return True
#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-














	#-#-#-#-#-#-#-#-#-#-#-#-#
	#  main loop functions  #
	#-#-#-#-#-#-#-#-#-#-#-#-#

	def send(self, message):
		if DEBUG_LEVEL > 2:
			print "WorkProcess: send()"
		self.output_queue.put(message)
#
	def send_data(self, message):
		if DEBUG_LEVEL > 2:
			print "WorkProcess: send_data()"
		#print str(message[0]) + "\t"  + str(message[1])
		#self.output_pipe.send(message)
		self.send((CMD_OUTPUT, message))
		#if DEBUG_LEVEL > 0:
			#print "WorkProcess: send_data(). done"
#





	def _kill(self):
		# WARNING Should close open files!
		self.status == STATUS_KILL
		self._status()
		self.send((CMD_KILL, True))
#





	def _status(self, req = None):
		if DEBUG_LEVEL > 0:
			print "WorkProcess: _status()"
		if req == None or req[0] == None:
			data = {}
			
			data['status'] = self.status
			data['go'] = {}
			data['ac'] = {}
			data['query'] = {}
			data['ss'] = {}
			data['output'] = {}
			
			data['go']['ok'] = self.ok_go
			data['go']['ok_params'] = self.ok_params_go
			data['ac']['ok'] = self.ok_ac
			data['ac']['ok_params'] = self.ok_params_ac
			data['query']['ok'] = self.ok_query
			data['query']['ok_params'] = self.ok_params_query
			data['output']['ok'] = self.ok_output
			data['output']['ok_params'] = self.ok_params_output
			data['ss']['ok'] = self.ok_ss
			data['ss']['ok_params'] = self.ok_params_ss
			self.send((CMD_STATUS, data))
#
		
		#if self.print_output:
			#print "----------GO----------"
			#if self.ok_go:
				#print "GO is ok"
				#print "GO loaded: " + str(self.go_filename)
				#print "Nodes: " + str(self.go.node_num()) + ". Edges: " + str(self.go.edge_num())
			#else:
				#print "GO is not ok"
			#print "----------AC----------"
			#if self.ok_ac:
				#print "AC is ok"
				#print "AC loaded: " + str(self.ac_filename)
				#print "Nodes: " + str(len(self.ac.annotations))
			#else:
				#print "AC is not ok"
			#print "----------Query----------"
			#if self.ok_query:
				#print "Query is ok"
				#if self.query_from == QUERY_FROM_GUI:
					#print "Query loaded from gui"
				#elif self.query_from == QUERY_FROM_FILE:
					#print "Query loaded from file " + str(self.query_filename)
				#elif  self.query_from == QUERY_FROM_AC:
					#print "Query loaded from ac"
			#else:
				#print "Query is not ok"
			#print "----------Output----------"
			#if self.ok_output:
				#print "Output is ok"
				#if self.query_from == OUTPUT_TO_GUI:
					#print "Output to gui"
				#elif self.query_from == OUTPUT_TO_FILE:
					#print "Output to file " + str(self.output_filename)
			#else:
				#print "Output is not ok"
			#print "----------SS----------"
			#if self.ok_ss:
				#print "SS is ok"
				#print "SS: " + str(self.ss_name)
				#print "MS: " + str(self.ms_name)
				#print "Ontology selected: " + str(self.ss_ontology)
			#else:
				#print "SS is not ok"
			#print "--------------------"
			#print "Current status"
#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-











	#-#-#-#-#-#-#-#-#-#-#-#-#
	# data parsing routines #
	#-#-#-#-#-#-#-#-#-#-#-#-#
	
	def _set(self, data):
		if DEBUG_LEVEL>0:
			print "WorkProcess: _set()"
		if data[0] == CMD_LOAD_GO:
			self.set_GO(data[1:len(data)])
		elif data[0] == CMD_LOAD_AC:
			self.set_AC(data[1:len(data)])
		elif data[0] == CMD_SET_QUERY:
			self.set_query(data[1:len(data)])
		elif data[0] == CMD_SET_OUTPUT:
			self.set_output_params(data[1:len(data)])
		elif data[0] == CMD_SET_SS:
			self.set_SS(data[1:len(data)])
		else:
			print "Unknown set request. Ignoring request..."
		self.status = STATUS_WAIT
#








	####-#-#-#-#-#-#-#-#-#
	#### load_AC routine #
	####-#-#-#-#-#-#-#-#-#

	def set_AC(self, data): #### data format: (filename, other) other = (file format, file format params)
		if DEBUG_LEVEL>1:
			print "WorkProcess: set_AC()"
		self.send((CMD_LOAD_AC, ANSWER_PROCESSING))
		self.ok_ac = False
		self.ok_params_ac = False
		
		data = data[0] #### WARNING: Why?
		
		self.ac_filename = data['filename']
		self.ac_filetype = data['type']
		self.ac_filetypeparams = {}
		if 'params' in data:
			self.ac_filetypeparams = data['params']
		self.ok_params_ac = True
		#self.ac = None # can remove
		#self.ok_ss = False # can remove
		#self.ok_query = False # can remove
		self.load_AC()
#





	def load_AC(self):
		if DEBUG_LEVEL>1:
			print "WorkProcess: load_AC()"
		if self.ok_ac:
			return True
		self.ok_query = False # important
		self.ok_ss = False # important
		self.ac = None # important
		if self.ok_params_ac:
			#try:
			self.ac = AnnotationCorpus.AnnotationCorpus(self.go)
			if self.ac.parse(str(self.ac_filename), self.ac_filetype, self.ac_filetypeparams):
				if self.go == None:
					self.ok_ac = False
				elif self.ac.sanitize():
					self.ok_ac = True
			#except Exception:
				#pass
		if self.ok_ac:
			self.send((CMD_LOAD_AC, ANSWER_PROCESSED, True, len(self.ac.annotations), len(self.ac.reverse_annotations)))
			return True
		else:
			self.send((CMD_LOAD_AC, ANSWER_PROCESSED, False))
			return False
#







	####-#-#-#-#-#-#-#-#-#
	#### load_GO routine #
	####-#-#-#-#-#-#-#-#-#
	
	def set_GO(self, data): #### data format: (filename, other) other = (file format, file format params)
		if DEBUG_LEVEL>1:
			print "WorkProcess: set_GO()"
		self.output_queue.put((CMD_LOAD_GO, ANSWER_PROCESSING, True))
		self.ok_go = False
		self.ok_params_go = False
		self.go = None # can remove
		
		data = data[0] #### WARNING: Why?
		
		if data == None or len(data) == 0:
			self.go_filename = None
			self.go_parameters = {}
		else:
			self.go_filename = data['filename']
			self.go_parameters = data
		self.ok_params_go = True
		self.ok_ac = False # can remove
		self.ok_ss = False # can remove
		self.ok_query = False # can remove
		self.load_GO()
#



	def load_GO(self):
		if DEBUG_LEVEL>1:
			print "WorkProcess: load_GO()"
		if self.ok_go:
			return True

		self.ok_ac = False # important
		self.ok_ss = False # important
		self.ok_query = False # important
		self.go = None
		if self.ok_params_go:
			try:
				#add "iffilename is none... then use the std file! Or is it already included?
				#program_dir = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
				#program_dir = '.' # use this with py2exe to build a working binary
				#go_file = program_dir + "/data/GO_2012-02-24.obo-xml.gz"
				self.go = GeneOntology.load(str(self.go_filename), 	self.go_parameters) # convert unicode to str!!!
			except Exception:
				#print "Load Gene Ontology. Exception loading " + str(self.go_filename)
				pass
			if not self.go == None:
				self.ok_go = True
				self.send((CMD_LOAD_GO, ANSWER_PROCESSED, True, self.go.node_num(), self.go.edge_num()))
				return True
		self.ok_go = False
		self.send((CMD_LOAD_GO, ANSWER_PROCESSED, False))
		return False
#







	####-#-#-#-#-#-#-#-#-#
	#### set_query       #
	####-#-#-#-#-#-#-#-#-#

	def set_query(self, data): #### data format: (query from, query_params) query_params: none (fom ac), (type, filename) (from file), type (from gui)
		if DEBUG_LEVEL>0:
			print "WorkProcess: set_query()"

		self.ok_query = False
		self.ok_params_query = False
		self.query = None
		
		if len(data) > 1:
			#print "WorkProcess: set_query() using data[1]"
			self.query = data[1]
			self.ok_query = True

		data = data[0] #### WARNING: Why?

		self.query_from = data['source']
		self.query_type = data['type']

		if self.query == None and 'query' in data:
			#print "WorkProcess: set_query() using data['query']"
			self.query = data['query'] #### NOTE 'query' should not be part of parameters
			self.ok_query = True

		if self.query_from == QUERY_FROM_AC:
			#print "Query from AC selected."
			self.ok_params_query = True
		elif self.query_from == QUERY_FROM_FILE:
			#print "Query from FILE selected."
			self.query_filename = data['filename']
			#self.query_filetype = data[3]
			#self.query_filetypeparams = data[4]
			self.ok_params_query = True
		elif self.query_from == QUERY_FROM_GUI:
				##print "Query from GUI selected."
				self.ok_params_query = True
		else:
			pass
#






	def load_query(self):
		if DEBUG_LEVEL>0:
			print "WorkProcess: load_query()"
		if self.ok_query:
			if DEBUG_LEVEL>0:
				print "WorkProcess: load_query() already ok"
			return True
		self.query = None
		if not self.ok_params_query:
			if DEBUG_LEVEL>0:
				print "WorkProcess: load_query() params not ok"
			return False
		if self.query_from == QUERY_FROM_GUI:
			if DEBUG_LEVEL>0:
				print "WorkProcess: load_query() Wrong"
			self.query = self.start_data[0]; #### DANGER expecting query as start parameters!
		elif self.query_from == QUERY_FROM_FILE:
			if DEBUG_LEVEL>0:
				print "WorkProcess: load_query() load from file"
			self.build_query_from_file()
		elif self.query_from == QUERY_FROM_AC:
			if DEBUG_LEVEL>0:
				print "WorkProcess: load_query() load from AC"
			self.build_query_from_ac()
		if not self.query == None:
			if DEBUG_LEVEL>0:
				print "WorkProcess: load_query() fixing pairs_number"
			if self.query_type == QUERY_PAIRS:
				self.query_pairs_number = len(self.query) 
			elif self.query_type == QUERY_LIST:
				self.query_pairs_number = len(self.query) * (len(self.query)-1) / 2
			self.ok_query = True
			return True
		if DEBUG_LEVEL>0:
			print "WorkProcess: load_query() Returning false"
		return False
#








	####-#-#-#-#-#-#
	#### set_SS    #
	####-#-#-#-#-#-#

	def set_SS(self, data): #### data format: (ss measure, ss measure params, mixing strat., mixing strat. params, ontology)
		if DEBUG_LEVEL>1:
			print "WorkProcess: set_SS()"
			print "NOTE. WorkProcess: Should check whether SS measure has changed or not. If so, rebuild the SS object."
		self.ok_ss = False
		self.ok_params_ss = False
		
		data = data[0] #### WARNING Why??
		self.ss_name = data['measure']
		self.ms_name = data['mixing_strategy']
		self.ss_ontology = data['ontology']

		self.ok_params_ss = True
		if self.ms_name == None:
			self.ok_params_ss = False
		if self.ss_name == None:
			self.ok_params_ss = False
		if self.ss_ontology == None:
			self.ok_params_ss = False

		if self.ok_params_ss:
			self.send((CMD_SET_SS, True))
		else:
			self.send((CMD_SET_SS, False))
#



	
	def _get(self, data):
		if DEBUG_LEVEL>0:
			print "WorkProcess: _get()"
		if data[0] == CMD_GET_AC:
			if DEBUG_LEVEL>0:
				print "WorkProcess: _get() CMD_GET_AC"
			if data[1] == CMD_GET_AC_OBJECTS:
				results = None
				if not self.ac == None:
					results = self.ac.obj_set.keys()
				self.send((CMD_GET, CMD_GET_AC, CMD_GET_AC_OBJECTS, results))
			elif data[1] == CMD_GET_AC_OBJECTS_NUMBER:
				if DEBUG_LEVEL>0:
					print "WorkProcess: _get() CMD_GET_AC_OBJECTS_NUMBER"
				results = None
				if not self.ac == None:
					results = len(self.ac.obj_set.keys())
				self.send((CMD_GET, CMD_GET_AC, CMD_GET_AC_OBJECTS_NUMBER, results))
			elif data[1] == CMD_GET_AC_TERMS:
				results = None
				if not self.ac == None:
					results = self.ac.term_set.keys()
				self.send((CMD_GET, CMD_GET_AC, CMD_GET_AC_TERMS, results))
			elif data[1] == CMD_GET_AC_TERMS_NUMBER:
				if DEBUG_LEVEL>0:
					print "WorkProcess: _get() CMD_GET_AC_TERMS_NUMBER"
				results = None
				if not self.ac == None:
					results = len(self.ac.term_set.keys())
				self.send((CMD_GET, CMD_GET_AC, CMD_GET_AC_TERMS_NUMBER, results))
		elif data[0] == CMD_GET_PARAMS:
			if data[1] == CMD_GET_PARAMS_AC:
				results = {}
				results['filename'] = self.ac_filename
				results['type'] = self.ac_filetype
				results['params'] = self.ac_filetypeparams
			elif data[1] == CMD_GET_PARAMS_GO:
				results = self.go_parameters
				results['filename'] = self.go_filename
			elif data[1] == CMD_GET_PARAMS_SS:
				results = {}
				results['ontology'] = self.ss_ontology
				results['measure'] = self.ss_name
				results['mixing_strategy'] = self.ms_name
			elif data[1] == CMD_GET_PARAMS_OUTPUT:
				results = {}
				results['to'] = self.output_to
				results['filename'] = self.output_filename
				results['filetype'] = self.output_filetype
				results['params'] = self.output_fileparams
				results['fileparams'] = None
			elif data[1] == CMD_GET_PARAMS_QUERY:
				results = {}
				results['source'] = self.query_from
				results['type'] = self.query_type
				results['filename'] = self.query_filename
			else:
				print "Unknown get params request. Ignoring request: " + str(data[1])
				return
			self.send((CMD_GET, data[0], data[1], results))
		else:
			print "Unknown get request. Ignoring request: " + str(data[0])
#



	def build_ss(self): # NOTE better if built when start is called. This way the parameters can be varied multiple times without rebuilding the object (potentially expensive!)
		if DEBUG_LEVEL>1:
			print "WorkProcess: build_ss()"
		if self.ok_ss:
			return True
		self.ss = None
		if self.ok_params_ss and self.ok_ac and self.ok_go:
			try:
				ssu = SemSimUtils.SemSimUtils(self.ac, self.go)
				ssu.det_IC_table()
				self.ss = ObjSemSim.ObjSemSim(self.ac, self.go, self.ss_name, self.ms_name, ssu)
				self.ss.TSS.format_and_check_data = False # NOTE Enhance performances disabling internal checks
				self.ok_ss = True
				return True
			except Exception:
				pass
		return False
#








	####-#-#-#-#-#-#-#-#-#-#-#
	#### set_output_params   #
	####-#-#-#-#-#-#-#-#-#-#-#
	
	def set_output_params(self, data): #### data format: (output_to, output_params) output_params: for file: (filename, type, params), for gui: (how many scores)
		if DEBUG_LEVEL>1:
			print "WorkProcess: set_output_params()"
			
		data = data[0] #### WARNING Why??
		
		self.ok_output = False
		self.ok_params_output = False
		self.output_fileparams = None
		self.output_params = None
		self.output_filters = None
		self.filter_0 = False
		self.filter_less = 0
		self.filter_None = False
		
		self.output_to = data['to']
		if 'params' in data:
			self.output_fileparams = data['params']
			self.output_params = data['params']
		if 'filter' in data:
			self.output_filters = data['filter']
			if 'None' in self.output_filters:
				self.filter_None = self.output_filters['None']
			if '0' in self.output_filters:
				self.filter_0 = self.output_filters['0']
			if 'less' in self.output_filters and self.output_filters['less']:
				self.filter_less = self.output_filters['value']
			
		if self.output_to == OUTPUT_TO_FILE:
			self.output_filename = data['filename']
		elif  self.output_to == OUTPUT_TO_GUI:
			self.output_filename = None
			#self.output_params = data['params']
		self.ok_params_output = True
		if self.ok_params_output:
			self.send((CMD_SET_OUTPUT, True))
		else:
			self.send((CMD_SET_OUTPUT, False))
#






	def init_output(self):
		if DEBUG_LEVEL>1:
			print "WorkProcess: init_output()"
		self.ok_output = False # do this to force reinit. of output module
		#if self.ok_output:
			#return
		if self.ok_params_output:
			if self.output_to == OUTPUT_TO_FILE:
				self.output_file = open(self.output_filename, 'w')
			elif self.output_to == OUTPUT_TO_GUI:
				self.output_buffer = [[]] * self.CURRENT_BUFFER_SIZE
				self.query_pairs_saved = 0
			else:
				return False
			self.ok_output = True
			return True
		else:
			print "WorkProcess: init_output() NO PARAMS OK!!!"
		return False
#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-








	#-#-#-#-#-#-#-#-#-#-#-#
	# utility functions   #
	#-#-#-#-#-#-#-#-#-#-#-#
	
	def build_query_from_file(self):
		h = open(self.query_filename,'r')
		self.query = []
		for line in h:
			line = line.rstrip('\n')
			line = line.rstrip('\r')
			if self.query_type == 0: #pairs
				line = line.rsplit('\t')
				self.query.append((line[0], line[1]))
			else: # list
				self.query.append(line)
		h.close()
#

	def build_query_from_ac(self):
		self.query = None
		if not self.ok_ac:
			return False
		self.query = []
		for line in self.ac.annotations:
			self.query.append(line)
		return True
#



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

	def _init_structures(self):
		if DEBUG_LEVEL>1:
			print "WorkProcess: _init_structures()"
		if not self.ok_params_go or not self.ok_params_ac or not self.ok_params_ss or not self.ok_params_query or not self.ok_params_output:
			if DEBUG_LEVEL>1:
				print "WorkProcess: _init_structures(). Something wrong in the parameters."
			return False
		if not self.load_GO():
			if DEBUG_LEVEL>1:
				print "WorkProcess: _init_structures(). Something wrong in the GO component."
			return False
		if not self.load_AC():
			if DEBUG_LEVEL>1:
				print "WorkProcess: _init_structures(). Something wrong in the AC component."
			return False
		if not self.build_ss():
			if DEBUG_LEVEL>1:
				print "WorkProcess: _init_structures(). Something wrong in the SS component."
			return False
		if not self.load_query():
			if DEBUG_LEVEL>1:
				print "WorkProcess: _init_structures(). Something wrong in the query component."
			return False
		if not self.init_output():
			if DEBUG_LEVEL>1:
				print "WorkProcess: _init_structures(). Something wrong in the output component."
			return False
		if not self.ok_go or not self.ok_ac or not self.ok_ss or not self.ok_query or not self.ok_output:
			if DEBUG_LEVEL>1:
				print self.ok_go
				print self.ok_ac
				print self.ok_ss
				print self.ok_query
				print self.ok_output
				print "WorkProcess: _init_structures(). Something still wrong."
			return False
		return True
#


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#	MAIN ROUTINES	MAIN ROUTINES	MAIN ROUTINES	MAIN ROUTINES	MAIN ROUTINES	MAIN ROUTINES	MAIN ROUTINES








	#-#-#-#-#-#-#-#-#-#-#
	# _start function   #
	#-#-#-#-#-#-#-#-#-#-#

	def _start(self, data):
		if DEBUG_LEVEL>0:
			print "WorkProcess: _start()"
		self.status = STATUS_RUN
		self.start_data = data
		self.send((CMD_START, ANSWER_PROCESSING, RESULT_OK, {'go':self.ok_go , 'ac':self.ok_ac, 'ss':self.ok_ss, 'query':self.ok_query, 'output':self.ok_output}))
		# initialize modules
		if not self._init_structures():
			if DEBUG_LEVEL>0:
				print "WorkProcess: _start(). Bad Parameters"
			self.send((CMD_START, ANSWER_PROCESSED, RESULT_BAD, {'go':self.ok_go , 'ac':self.ok_ac, 'ss':self.ok_ss, 'query':self.ok_query, 'output':self.ok_output}))
			self._stop()
		
		else:
			self.send((CMD_START, ANSWER_PROCESSED, RESULT_OK))
			self.obj1_pos = 0
			self.obj2_pos = 1
			self.obj_pos = 0
			self.pairs_done = 0
			#print "First Params:\tobj1_pos: " + str(self.obj1_pos) + "\tobj2_pos: " + str(self.obj2_pos) + "\tpairs_done: " + str(self.pairs_done) + "\tquery_pairs_saved: " + str(self.query_pairs_saved)
			self.CURRENT_BUFFER_SIZE = 100
			self._calculate()
#





	#-#-#-#-#-#-#-#-#-#-#-#
	# _restart function   #
	#-#-#-#-#-#-#-#-#-#-#-#
	
	def _restart(self):
		self.status = STATUS_RUN
		self.send((CMD_START, ANSWER_PROCESSED, True, self.query_pairs_number))
		self._calculate()
#









	#-#-#-#-#-#-#-#-#-#-#
	# _stop function    #
	#-#-#-#-#-#-#-#-#-#-#
	
	# go, ac, query, output parameters and ss should be mantained as valid
	# output data should be flushed, as well as temporary variables.
	# output_update should be set to true to reinitialize output structures
	def _stop(self):
		if DEBUG_LEVEL>0:
			print "WorkProcess: _stop()"
		self.send((CMD_STOP, ANSWER_PROCESSING))
		if self.output_to == OUTPUT_TO_FILE and not self.output_file == None:
			self.output_file.close()
		self.output_file = None
		self.output_buffer = None
		self.output_update = True
		self.status = STATUS_WAIT
		self.send((CMD_STOP, ANSWER_PROCESSED, RESULT_OK))
#








	#-#-#-#-#-#-#-#-#-#-#
	# _pause function   #
	#-#-#-#-#-#-#-#-#-#-#
	
	#def pause(self):
		#self.ssprocess2gui_queue.put((CMD_PAUSE, True))
		#self.status = STATUS_PAUSE
#









	#-#-#-#-#-#-#-#-#-#-#-#
	# _calculate function #
	#-#-#-#-#-#-#-#-#-#-#-#

	def _calculate(self):
		if DEBUG_LEVEL>1:
			print "WorkProcess: _calculate()"
		if self.query_type == QUERY_LIST:
			for self.obj1_pos in range(self.obj1_pos, len(self.query)):
				#print str(self.query[self.obj1_pos])
				for self.obj2_pos in range(self.obj2_pos, len(self.query)):
					self._temp_ss = self.ss.SemSim(self.query[self.obj1_pos],self.query[self.obj2_pos], self.ss_ontology)
					self._dispatch_output()
					self.pairs_done += 1
					
					#if True: # self.obj2_pos < self.obj1_pos + 10:
						#print "\tobj1_pos: " + str(self.obj1_pos) + "\tobj2_pos: " + str(self.obj2_pos) + "\tpairs_done: " + str(self.pairs_done) + "\tquery_pairs_saved: " + str(self.query_pairs_saved) + "\tobj1: " + str(self.query[self.obj1_pos]) + "\tobj2: "  + str(self.query[self.obj2_pos]) + "\tscore: " + str(self._temp_ss) 

					#print self.pairs_done
					if self.obj2_pos == len(self.query) - 1:
						self.obj2_pos = self.obj1_pos + 2
						break
					if self.pairs_done%self.CHECK_INTERVAL == 0:
						self.obj2_pos = self.obj2_pos + 1
						return

		elif self.query_type == QUERY_PAIRS:
			for self.obj_pos in range(self.obj_pos,len(self.query)):
				self._temp_ss = self.ss.SemSim(self.query[self.obj_pos][0], self.query[self.obj_pos][1], self.ss_ontology)
				self._dispatch_output()
				self.pairs_done += 1
				#print self.pairs_done
				if self.pairs_done%self.CHECK_INTERVAL == 0:
					return

		self._flush_output()
		self._stop()
#











	#-#-#-#-#-#-#-#-#-#-#-#
	# _dispatch_output f. #
	#-#-#-#-#-#-#-#-#-#-#-#
	
	def _dispatch_output(self):
		# add here switch to ignore None if required
		if type(self._temp_ss) == float:
			if self.filter_0 and self._temp_ss == float(0):
				return
			if self.filter_less > 0 and self._temp_ss < self.filter_less:
				return
			self._temp_ss = str('%.4f' %self._temp_ss)
		else:
			if self.filter_None and self._temp_ss == None:
				return

		if self.query_type == QUERY_LIST:
			if self.output_to == OUTPUT_TO_GUI:
				self._send_output(str(self.query[self.obj1_pos]), str(self.query[self.obj2_pos]), str(self._temp_ss))
			else:
				self._write_output(str(self.query[self.obj1_pos]), str(self.query[self.obj2_pos]), str(self._temp_ss))

		elif self.query_type == QUERY_PAIRS:
			if self.output_to == OUTPUT_TO_GUI:
				self._send_output(str(self.query[self.obj_pos][0]), str(self.query[self.obj_pos][1]), str(self._temp_ss))
			else:
				self._write_output(str(self.query[self.obj_pos][0]), str(self.query[self.obj_pos][1]), str(self._temp_ss))
#






	#-#-#-#-#-#-#-#-#-#-#-#
	# _send_output func.  #
	#-#-#-#-#-#-#-#-#-#-#-#
	
	def _send_output(self, a, b, c):
		#print (str(a),str(b),str(c))
		#if self.store_all:
			#self.last_added += 1
			#self.buffer[self.last_added] = ((str(a),str(b),str(c)))
			#if self.last_added - (self.last_sent + 1) > self.MAX_CACHE_SIZE:
				#self.ssprocess2gui_pipe.send(self.buffer[self.last_sent+1: self.last_added+1])
				#self.last_sent = self.last_added
		#else:
		self.output_buffer[self.query_pairs_saved] = (a,b,c)
		self.query_pairs_saved += 1
		if self.query_pairs_saved == len(self.output_buffer):
			if DEBUG_LEVEL>1:
				print "WorkProcess: _send_output(). Cycle"
			self.send_data(self.output_buffer)
			self.CURRENT_BUFFER_SIZE = self.MAX_BUFFER_SIZE
			self.output_buffer = [[]] * self.CURRENT_BUFFER_SIZE
				#self.query_pairs_saved = 0
			self.query_pairs_saved = 0
#






	#-#-#-#-#-#-#-#-#-#-#-#
	# _flush_output func. #
	#-#-#-#-#-#-#-#-#-#-#-#
	
	def _flush_output(self):
		if DEBUG_LEVEL>0:
				print "WorkProcess: _flush_output"
		#if self.store_all:
			#if self.last_sent == self.last_added:
				#return
			#self.ssprocess2gui_pipe.send(self.buffer[self.last_sent+1: self.last_added+1])
			#self.last_sent = self.last_added
		#else:
		if not self.output_to == OUTPUT_TO_GUI:
			return
		if self.query_pairs_saved == 0:
			return
		self.send_data(self.output_buffer[0: self.query_pairs_saved])
		self.output_buffer = [[]] * self.CURRENT_BUFFER_SIZE
		self.query_pairs_saved = 0
#






	#-#-#-#-#-#-#-#-#-#-#-#
	# _write_output func. #
	#-#-#-#-#-#-#-#-#-#-#-#
	
	def _write_output(self, a, b, c):
		self.output_file.write(str(a) + "\t" + str(b) + "\t" + str(c) + "\n")
		return True
#
