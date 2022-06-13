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

#### Import basic packages
import sys
import os
import copy

#### Import fastSemSim modules
from fastsemsim.Ontology import GeneOntology
from fastsemsim.Ontology import AnnotationCorpus
import fastsemsim.SemSim
from fastsemsim.SemSim import ObjSemSim

#### Include this to work with processes and threads
import threading
import multiprocessing 
import WorkProcess

#### Import wxWidgets
import wx
import wx.lib.newevent

#### Import modules
from GOGui import *
from ACGui import *
from QueryGui import *
from SSGui import *
from OutputGui import *
from ConfigGui import *

#--------------------------------------------------#

DEBUG_LEVEL = 0



# WorkProcess is designed as a state machine. Here are the possible status
STATUS_BASE = 55
STATUS_INIT = STATUS_BASE + 0 # performing reset.
STATUS_WAIT = STATUS_BASE + 1 # waiting for requests.
STATUS_RUN = STATUS_BASE + 2 # computation in progress.
STATUS_EXIT = STATUS_BASE + 3
#







class fastSemSimGui(wx.Frame):

# control & status flags
	status_gui = STATUS_INIT # current GUI status
	status = None # most recently known process status
#


# params
	config_file = None
	
	params_GO = {'filename':None, 'type':None, 'ignore':{}}
	params_AC = {'filename':None, 'type':None, 'params':{'filter':{}}}
	params_SS = {'ontology':None, 'measure':None, 'mixing_strategy':None}
	params_query = {'type':None, 'source':None}
	params_output = {'to':None, 'params':{}, 'filter':{}, 'filename':None}

	GO_status = False
	AC_status = False
	SS_status = False
	query_status = False
	output_status = False
#


# multiprocessing data
	running = False # whether the background process is running
	ssprocess = None
#









#-------------------------------------------------------------------------------#
#     This methods initialize the main GUI and the subprocesses/subthreads      #
#-------------------------------------------------------------------------------#

	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "fastSemSimGui v.2 - Marco Mina", pos = wx.DefaultPosition, size = wx.Size( 652,571 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: init()"

# Set constants. Load images
		self.programdirectory = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
		#self.programdirectory = 'images' # use this with py2exe to build a working binary

		self.icons_folder 	= self.programdirectory + '/icons/32x32/'

		self.Ok_pic 		= self.icons_folder + 'V.png'
		self.Warning_pic 	= self.icons_folder + 'W.png'
		self.Advanced_pic 	= self.icons_folder + 'advanced.png'
		self.work_pic 		= self.icons_folder + 'work.png'
		self.query_pic 		= self.icons_folder + 'query.png'
		self.GO_ok_pic 		= self.icons_folder + 'GO_ok.png'
		self.GO_warn_pic 	= self.icons_folder + 'GO_warn.png'
		self.output_pic 	= self.icons_folder + 'output.png'
		self.SS_pic 		= self.icons_folder + 'tweak.png'
	
		self.font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
		self.font.SetPointSize(9)

		self.status_images = wx.ImageList(32,32)

		self.status_images.Add(wx.Bitmap(self.Ok_pic)) # 0
		self.status_images.Add(wx.Bitmap(self.Warning_pic)) # 1
		self.status_images.Add(wx.Bitmap(self.Advanced_pic)) # 2
		self.status_images.Add(wx.Bitmap(self.work_pic)) # 3
		self.status_images.Add(wx.Bitmap(self.query_pic)) # 4
		self.status_images.Add(wx.Bitmap(self.output_pic)) # 5
		self.status_images.Add(wx.Bitmap(self.SS_pic)) # 6
		self.status_images.Add(wx.Bitmap(self.GO_ok_pic)) # 7
		self.status_images.Add(wx.Bitmap(self.GO_warn_pic)) # 8

		self._init_communication()
		self._init_handles()
		
# Build main GUI
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		fastSemSim_sizer_1 = wx.BoxSizer( wx.VERTICAL )
		self.fastSemSim_listbook = wx.Listbook( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LB_DEFAULT )
		self.fastSemSim_listbook.SetImageList(self.status_images)
		
		fastSemSim_sizer_1.Add( self.fastSemSim_listbook, 1, wx.EXPAND |wx.ALL, 5 ) #### This section was in the end

		self.fastSemSim_statusbar = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
		self.InitMenu()

		self.SetSizer( fastSemSim_sizer_1 )
		self.Layout()
		self.Centre( wx.BOTH )

		#### Panel Controls
		#self.query_panel = ControlsPanel(self, self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		#self.fastSemSim_listbook.AddPage( self.query_panel, u"Controls", True )
		#self.fastSemSim_listbook.SetPageImage(0, 3)

		#### Panel GO
		self.GO_panel = GOPanel(self, self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.fastSemSim_listbook.AddPage( self.GO_panel, u"Gene Ontology", False )
		self.fastSemSim_listbook.SetPageImage(0, 8)
		#### Panel AC
		self.AC_panel = ACPanel(self, self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.fastSemSim_listbook.AddPage( self.AC_panel, u"Annotation Corpus", False )
		self.fastSemSim_listbook.SetPageImage(1, 8)
		#### Panel SS
		self.SS_panel = SSPanel(self, self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.fastSemSim_listbook.AddPage( self.SS_panel, u"Semantic Similarity", False )
		self.fastSemSim_listbook.SetPageImage(2, 6)
		#### Panel Query
		self.query_panel = QueryPanel(self, self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.fastSemSim_listbook.AddPage( self.query_panel, u"Query", False )
		self.fastSemSim_listbook.SetPageImage(3, 4)
		#### Panel Output Ctrl
		self.output_ctrl_panel = OutputCtrlPanel(self, self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.fastSemSim_listbook.AddPage( self.output_ctrl_panel, u"Output Settings", False )
		self.fastSemSim_listbook.SetPageImage(4, 5)
		#### Panel Output
		#self.output_window = OutputPanel(self, self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		#self.fastSemSim_listbook.AddPage( self.output_window, u"Output", False )
		##self.fastSemSim_listbook.SetPageImage(5, 5)

		#self.SetStatus(STATUS_WAIT)
		#self.activateGoCmd()
		#self.OnAnyUpdate()
		# Init background process and communication module

		self._init_process()
		self._reset()
		#self.reset() # should restart gui side section. Do no touch communication channels or background process
#



	def InitMenu(self):
		#self.ID_NEW=wx.NewId()
		self.ID_OPEN=wx.NewId()
		self.ID_QUIT=wx.NewId()
		self.ID_SAVE=wx.NewId()
		self.ID_SAVE_AS=wx.NewId()
		self.MenuBar = wx.MenuBar()
		self.FileMenu = wx.Menu()
		#self.FileMenu.Append(self.ID_NEW, 'New', 'New')
		self.FileMenu.Append(self.ID_OPEN, 'Open Configuration...', 'Open Configuration')
		self.FileMenu.Append(self.ID_SAVE, 'Save Configuration', 'Save Configuration')
		self.FileMenu.Append(self.ID_SAVE_AS, 'Save Configuration As...', 'Save Configuration As')
		self.FileMenu.Append(self.ID_QUIT, 'Quit', 'Quit application')
		self.MenuBar.Append(self.FileMenu,'&File')
		self.SetMenuBar(self.MenuBar)
		self.Bind(wx.EVT_MENU, self.OnMenuQuit, id=self.ID_QUIT)
		self.Bind(wx.EVT_MENU, self.OnMenuOpen, id=self.ID_OPEN)
		self.Bind(wx.EVT_MENU, self.OnMenuSave, id=self.ID_SAVE)
		self.Bind(wx.EVT_MENU, self.OnMenuSaveAs, id=self.ID_SAVE_AS)
		#self.Bind(wx.EVT_MENU, self.OnMenuNew, id=self.ID_NEW)
#



	#def OnMenuNew(self, event):
		#print "OnMenuNew still to be implemented."
#





	def OnMenuQuit(self, event):
		self.OnQuit(None)
#






	def OnMenuOpen(self, event):
		dialog = wx.FileDialog(None, style = wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			self.fastSemSim_statusbar.SetStatusText("Loading config file: " + str(dialog.GetPath()))
			#if not self.config_file == None:
				#self.OnMenuNew(None)
			self.config_file = dialog.GetPath()
			self.loadConfigGui = LoadConfigGui(self)
#









	def OnMenuSave(self, event):
		if self.config_file == None:
			return self.OnMenuSaveAs(event)
		self.saveConfigGui = SaveConfigGui(self)
#




	def OnMenuSaveAs(self, event):
		dialog = wx.FileDialog(None, style = wx.SAVE|wx.OVERWRITE_PROMPT)
		if dialog.ShowModal() == wx.ID_OK:
			self.config_file = dialog.GetPath()
			self.OnMenuSave(event)
#






	def _reset(self):
		self.update()
#


	def _init_process(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: _init_process()"
		if not self.running:
			if DEBUG_LEVEL>1:
				print "Starting process"
			self.ssprocess = []
			self.ssprocess.append(None)
			#self._init_communication()
			self.ssprocess[0] = WorkProcess.WorkProcess(self.gui2ssprocess_queue, self.ssprocess2gui_queue, self.ssprocess2gui_pipe, self.gui2ssprocess_pipe)
			self.running = True
			self.ssprocess[0].start()
			return True
		else:
			if DEBUG_LEVEL>1:
				print "Process already in execution"
			return False
#



	def _init_communication(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: _init_communication()"
		self.process_busy = False
		self.process_busy_lock = multiprocessing.Lock()
		self.process_busy_event = multiprocessing.Event()
		self.gui2ssprocess_queue = multiprocessing.Queue()
		self.ssprocess2gui_queue = multiprocessing.Queue()
		self.ssprocess2gui_pipe = None
		self.gui2ssprocess_pipe = None
		#self.gui2ssprocess_pipe, self.ssprocess2gui_pipe = multiprocessing.Pipe()
		self._init_thread()
#



	def _init_thread(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: _init_thread()"
		self.communication_thread = CommunicationThread(self)
		self.communication_thread.start()
#




	def _init_handles(self):
		#self.start_check_handle = None
		self.status_handle = None
		#if not self.start_check_handle == None:
			#communication_thread.unregister_callback(self.start_check_handle)
			#self.start_check_handle = None
		#if not self.status_handle == None:
			#communication_thread.unregister_callback(self.status_handle)
			#self.status_handle = None

		self.status_handle = self.communication_thread.register_callback(self.EVT_CUSTOM_STATUS, self.OnStatusEvent)
		self.generic_handle = self.communication_thread.register_callback(self.EVT_CUSTOM_GENERIC, self.OnGenericEvent)
		self.start_handle = self.communication_thread.register_callback(self.EVT_CUSTOM_START, self.OnStartEvent)
		self.stop_handle = self.communication_thread.register_callback(self.EVT_CUSTOM_STOP, self.OnStopEvent)
		self.output_handle = self.communication_thread.register_callback(self.EVT_CUSTOM_OUTPUT, self.OnOutputEvent)
		self.kill_handle = self.communication_thread.register_callback(self.EVT_CUSTOM_KILL, self.OnKillEvent)

		self.Bind(wx.EVT_CLOSE, self.OnQuit, id=self.GetId())
#




#######################################################
#######################################################
#################           EVENT HANDLING         ####
#######################################################
#######################################################


# When the GUI raises a Quit Request
	def OnQuit(self, event):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: OnQuit()"
		self.status_gui = STATUS_EXIT # WARNING Should be done in a dedicated function, perhaps disabling all the commands
		self._kill_process()
		if not event == None:
			event.Veto()
#






# When the background process terminates. If the gui is in termination phase, then exit.
# Otherwise reinitialize the background process
	def OnKillEvent(self, event):
		self.running = False
		if DEBUG_LEVEL>1:
			print "fastSemSimGui: OnKillEvent()"
		if self.status_gui == STATUS_EXIT:
			if DEBUG_LEVEL>1:
				print "fastSemSimGui: OnKillEvent(). status is STATUS_EXIT. Quitting."
			self.Destroy()
			if not self.communication_thread == None:
				self.ssprocess2gui_queue.put((-1,))
				self.communication_thread.join()
				self.communication_thread = None
			sys.exit()
		else:
			if DEBUG_LEVEL>1:
				print "fastSemSimGui: OnKillEvent(). status is NOT STATUS_EXIT. Respawning background process..."
			self._init_process() # WARNING what about updating the status?
			self.update()
#









	def OnStatusEvent(self, event):
		if DEBUG_LEVEL>1:
			print "fastSemSimGui: OnStatusEvent()"
			#print str(event.data)
		data = event.data
		self.GO_status = data[1]['go']['ok']
		self.AC_status = data[1]['ac']['ok']
		self.ss_status = data[1]['ss']['ok_params']
		self.query_status = data[1]['query']['ok_params']
		self.output_status = data[1]['output']['ok_params']
		self.status = data[1]['status']
		self._update()
		event.Skip()
#



	def OnGenericEvent(self, event):
		if DEBUG_LEVEL>1:
			print "fastSemSimGui: OnStatusEvent()"
		#event.Skip()
#



	def OnStartEvent(self, event):
		if DEBUG_LEVEL>1:
			print "fastSemSimGui: OnStartEvent()"
		#print len(event.data)
		data = event.data
		if data[0] == WorkProcess.CMD_START:
			if data[1] == WorkProcess.ANSWER_PROCESSED:
				if data[2] == WorkProcess. RESULT_OK:
					if DEBUG_LEVEL>0:
						print "start_outcome: Computation start."
					if self.params_output['to'] == WorkProcess.OUTPUT_TO_FILE:
						self.query_panel.controls_log.AppendText("Output data redirected to the selected output file.\n")
					else:
						self.output_window = OutputWindow(self)
						self.output_window.Show(True)
				else:
					if DEBUG_LEVEL>0:
						print "start_outcome: Computation not started."
					if not data[3]['go']:
						self.query_panel.controls_log.AppendText("Check Gene Ontology.\nAborted\n")
					if not data[3]['ac']:
						self.query_panel.controls_log.AppendText("Check Annotation Corpus.\nAborted.\n")
					if not data[3]['query']:
						self.query_panel.controls_log.AppendText("Check query.\nAborted.\n")
					if not data[3]['output']:
						self.query_panel.controls_log.AppendText("Check output parameters.\nAborted.\n")
					if not data[3]['ss']:
						self.query_panel.controls_log.AppendText("Check operation parameters.\nAborted.\n")

			elif data[1] == WorkProcess.ANSWER_PROCESSING:
				if DEBUG_LEVEL>0:
					print "start_outcome: Start request is being processed."
				self.query_panel.controls_log.AppendText("Starting computation...\n")
				return
			elif data[1] == WorkProcess.ANSWER_IGNORED:
				if DEBUG_LEVEL>0:
					print "start_outcome: Start request ignored."
			else:
				if DEBUG_LEVEL>0:
					print "start_outcome: Unknown answer."
			self.update()
			#self.query_panel.controls_start_button.Enable()
			#self.query_panel.controls_start_button.SetLabel("Stop")
#






	def OnStopEvent(self, event):
		if DEBUG_LEVEL>1:
			print "fastSemSimGui: OnStopEvent()"
		data = event.data
		if data[0] == WorkProcess.CMD_STOP:
			if data[1] == WorkProcess.ANSWER_PROCESSED:
				if data[2] == WorkProcess. RESULT_OK:
					if DEBUG_LEVEL>0:
						print "stop_outcome: Computation stopped."
				else:
					if DEBUG_LEVEL>1:
						print "stop_outcome: Computation not stopped."
			elif data[1] == WorkProcess.ANSWER_PROCESSING:
				if DEBUG_LEVEL>1:
					print "stop_outcome: Stop request is being processed."
				return
			elif data[1] == WorkProcess.ANSWER_IGNORED:
				if DEBUG_LEVEL>1:
					print "stop_outcome: Stop request ignored."
			else:
				if DEBUG_LEVEL>1:
					print "stop_outcome: Unknown answer."
			self.query_panel.controls_log.AppendText("Computation terminated.\n-------------------------\n")
			self.update()
			#event.Skip()
#





	def OnOutputEvent(self, event):
		if DEBUG_LEVEL>1:
			print "fastSemSimGui: OnStatusEvent()"
		data = event.data
		if data[0] == WorkProcess.CMD_OUTPUT:
			data = data[1]
			#print "OUTPUT: " + str(data[0]) + "\t" + str(data[1]) + "\t" + str(data[2])
			#temp = ("\n".join([str(`num`) for num in data]))
			temp = "\n".join(["\t".join([str(b) for b in num]) for num in data])
			temp += "\n"
			self.output_window.output_text.AppendText(temp)
			#event.Skip()
#




	#def lock(self, wait = True):
		#if DEBUG_LEVEL>0:
			#print "fastSemSimGui: lock()"
		#self.process_busy_lock.acquire()
		#if not wait:
			#if self.process_busy:
				#self.process_busy_lock.release()
				#return False
		#else:
			#while self.process_busy:
				#self.process_busy_lock.release()
				#self.process_busy_event.wait()
				#self.process_busy_lock.acquire()
		#self.process_busy = True
		#self.process_busy_event.clear()
		#self.process_busy_lock.release()
		#return True
##

	#def unlock(self):
		#if DEBUG_LEVEL>0:
			#print "fastSemSimGui: unlock()"
		#self.process_busy_lock.acquire()
		#self.process_busy = False
		#self.process_busy_event.clear()
		#self.process_busy_event.set()
		#self.process_busy_lock.release()
		#return True
##


	def OnCheckProcessData(self,event):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: OnCheckProcessData()"
		#print "timer elapsed"
		data = event.data
		if data[0] == WorkProcess.CMD_OUTPUT:
			if DEBUG_LEVEL>0:
				print "fastSemSimGui: OnCheckProcessData(). CMD_OUTPUT received"
			#if self.ssprocess2gui_pipe.poll():
				#incoming_data = self.ssprocess2gui_pipe.recv()
			incoming_data = data[1]
			self.OnPipeData(incoming_data)
#


	#def OnCheckPipes(self,event):
		##print "timer elapsed"
		#if self.ssprocess2gui_pipe.poll():
			##print "data!"
			#tmp = self.ssprocess2gui_pipe.recv()
			#self.OnOutputData(tmp)
			#gaugerange = self.progress.GetRange()
			#self.progress.SetValue((float(self.superconta)/self.pairs_to_process)*gaugerange)
			
		#if self.lock(wait = False):
			##print "THIS CHECK"
			#if not self.ssprocess2gui_queue.empty():
				#data = self.ssprocess2gui_queue.get()
				##print data
				#if data[0] == WorkProcess.CMD_STOP:
					#if data[1]:
						#self.completed()
					#else:
						#self.completed()
			#self.unlock()
		#else:
			#pass

##




#################################

#'''
#Following routines handle communications with other processes 
#'''



			#print "Busy"
		
		#self.OnProgress()
		#self.OnCompleted()

	#def OnProgress(self):
		##print "Updating progress bar..."
		#gaugerange = self.progress.GetRange()
		#self.progress.SetValue(self.sspdone.value*gaugerange)
		##self.OnUpdateDone()

	#def OnCompleted(self):
		#if self.sscompleted.value == 1:
			#self.log_field.AppendText("Finished.\n")
			#print self.superconta
			##self.running = False
			#self.stop()


	def OnLogData(self, msg):
		t = msg.data
		self.log_field.AppendText(t)
		#self.OnUpdateDone()
#






	def OnAnyUpdate(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: OnAnyUpdate()"
		#self.statusgridbox.Fit(self)
		#self.gostatsbox.Fit(self)
		#self.commandbox.Fit(self)
		self.mainbox.Fit(self)
		#w,h = self.commandbox.GetSizeTuple()
		#self.mainbox.SetItemMinSize(self.commandbox, w,h)
		#self.GetBestSize()
		self.mainbox.Layout()
		















	def OnUpdateDone(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: OnUpdateDone()"
		pass
		#print "Update done"
		#self.update_event.set()
#






	def OnDefaultAction(self, event):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: OnDefaultAction()"
#











#CustomEvent_Set, EVT_CUSTOM_SET = wx.lib.newevent.NewEvent()
#CustomEvent_LoadAC, EVT_CUSTOM_LOAD_AC = wx.lib.newevent.NewEvent()
#CustomEvent_LoadGO, EVT_CUSTOM_LOAD_GO = wx.lib.newevent.NewEvent()
#CustomEvent_SetSS, EVT_CUSTOM_SET_SS = wx.lib.newevent.NewEvent()
#CustomEvent_SetOutput, EVT_CUSTOM_SET_OUTPUT = wx.lib.newevent.NewEvent()
#CustomEvent_SetQuery, EVT_CUSTOM_SET_QUERY = wx.lib.newevent.NewEvent()


	def OnPipeData(self, t):
		if DEBUG_LEVEL>1:
			print "fastSemSimGui: OnPipeData()"
		#print t
		
		frase = ""
		for i in t: #range(t[0], t[1] + 1):
			##frase = self.ssprocess[0].buffer[i][0] + "\t" + self.ssprocess[0].buffer[i][1] + "\t" + self.ssprocess[0].buffer[i][2] + "\n"
			#self.superconta += 1
			frase += str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\n"
			#self.OutputGui.output_field.AppendText(frase)
		self.output_window.output_text.AppendText(frase)
		#frase = ""
		#for i in t: #range(t[0], t[1] + 1):
			#frase = self.ssprocess[0].buffer[i][0] + "\t" + self.ssprocess[0].buffer[i][1] + "\t" + self.ssprocess[0].buffer[i][2] + "\n"
			#self.superconta += 1
			#frase = frase + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\n"
		#self.OutputGui.output_field.AppendText(frase)
		#print "Output data printed"
		#self.OnUpdateDone()
#



	def OnQueueData(self, data):
		if DEBUG_LEVEL>1:
			print "fastSemSimGui: OnQueueData()"
		if data[0] == WorkProcess.CMD_STOP:
			if data[1]:
				#def stop(self):
					self.completed()
					self.stop_process()
			else:
				self.SetStatus(1)
				self.completed()
		return True
#







	def _update(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: _update()"

		if self.status == WorkProcess.STATUS_WAIT:
			#print "STATUS_WAIT"
			self._unfreeze()

		if self.GO_status:
			self.fastSemSim_listbook.SetPageImage(1, 7)
			self.AC_panel._unfreeze()
		else:
			self.fastSemSim_listbook.SetPageImage(1, 8)
			self.AC_panel._freeze()
		if self.AC_status:
			self.fastSemSim_listbook.SetPageImage(2, 7)
		else:
			self.fastSemSim_listbook.SetPageImage(2, 8)

		self.GO_panel._update()
		self.AC_panel._update()
		self.SS_panel._update()
		self.query_panel._update()
		self.output_ctrl_panel._update()
		
		if not self.status == WorkProcess.STATUS_WAIT:
			#print "not STATUS_WAIT"
			self._freeze()
			
		if not self.GO_status:
			self.fastSemSim_statusbar.SetStatusText("Load a Gene Ontology (Gene Ontology panel)")
		elif not self.AC_status:
			self.fastSemSim_statusbar.SetStatusText("Load an Annotation Corpus (Annotation Corpus panel)")
		elif self.status == WorkProcess.STATUS_RUN:
			self.fastSemSim_statusbar.SetStatusText("Calculation in progress. Press Stop to abort (Query panel)")
		elif self.status == WorkProcess.STATUS_WAIT:
			self.fastSemSim_statusbar.SetStatusText("Use the Query panel to start the computation.")
#



	def update(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: update()"
		self.gui2ssprocess_queue.put((WorkProcess.CMD_STATUS, None))
#









	def _kill_process(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: _kill_process()"

		self.gui2ssprocess_queue.put((WorkProcess.CMD_KILL, None))
		
		#if not self.ssprocess[0] == None:
			#self.ssprocess[0].terminate() # DANGER abruptly terminates computation process. Should not use this
			#self.ssprocess[0] = None
		#self.running = False
#






	#def _start_process(self):
		#if DEBUG_LEVEL>0:
			#print "fastSemSimGui: _start_process()"
		##self.progress.SetValue(float(0))
		#if self.status == 1:
			#if DEBUG_LEVEL>1:
				#print "Job already running."
			#return False
		#if self.status == -1:
			#if DEBUG_LEVEL>1:
				#print "Process not initialized."
			#return False

		##self.lock()
		#self.gui2ssprocess_queue.put((WorkProcess.CMD_LOAD_SS, self.param_ss_measure, None, self.param_mixing_strategy, None, self.param_selected_GO))
		##data = self.ssprocess2gui_queue.get()
		#self.gui2ssprocess_queue.put((WorkProcess.CMD_STATUS, None))
		##data = self.ssprocess2gui_queue.get()
		##if self.query_from == WorkProcess.QUERY_FROM_GUI:
			##self.loadFromField()
		#self.start_outcome_handle = self.communication_thread.register_callback(self.start_outcome)
		#self.gui2ssprocess_queue.put((WorkProcess.CMD_START, self.param_query))

		##data = self.ssprocess2gui_queue.get()
		##self.unlock()
		##if data[0] == WorkProcess.CMD_START and data[1]:
			##self.pairs_to_process = data[2]
			##self.SetStatus(1)
		##else:
			##return False
		#return True
##
















	def SetStatus(self, status):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: SetStatus()"
		prev = self.status
		self.status = status
		if DEBUG_LEVEL>1:
			print "Status switched from " + str(prev) + " to " + str(self.status)
		
		if prev==-1 and not self.status==-1:
			#self.timer.Start(1000)
			pass
		
		#if self.show_pics:
			#if self.status == -1: # Fatal error
				#if not self.start_cmd is None:
					#self.start_cmd.Disable()
					#self.start_cmd.SetLabel("Start")
			#if self.status == 0: # Not running. Not Ready
				#if not self.start_cmd is None:
					#self.start_cmd.Disable()
					#self.start_cmd.SetLabel("Start")
			#if self.status == 1: # Not Running. Ready
				#if not self.start_cmd is None:
					#self.start_cmd.SetLabel("Start")
					#self.start_cmd.Enable()
			#if self.status == 2: # Running. (Completed)
				#if not self.start_cmd is None:
					#self.start_cmd.SetLabel("Stop")
					##self.statuspicture.SetBitmap(wx.Bitmap(self.s0_pic))
					#self.start_cmd.Enable()
#




























###################################################################################

	def loadFromField(self):
		h = self.QueryGui.text_query.GetValue()
		self.query= []
		h = h.splitlines()
		for line in h:
			if self.query_type == WorkProcess.QUERY_PAIRS:
				line = line.rsplit(' ')
				self.query.append((str(line[0]), str(line[1])))
			else:
				self.query.append(str(line))
		#print self.query
		return True
#








	def completed(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: completed()"
		self.progress.SetValue(self.progress.GetRange())
		self.running = False
		self.log_field.AppendText("Task completed.\n")
		self.SetStatus(1)
#



	def _freeze(self):
		self.GO_panel._freeze()
		self.AC_panel._freeze()
		self.query_panel._freeze()
		self.output_ctrl_panel._freeze()
		self.SS_panel._freeze()
#

	def _unfreeze(self):
		self.GO_panel._unfreeze()
		self.AC_panel._unfreeze()
		self.query_panel._unfreeze()
		self.output_ctrl_panel._unfreeze()
		self.SS_panel._unfreeze()
#




	def start(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: start()"

		self._freeze()
		
		if self.SS_panel.SS_to_send:
			self.SS_panel.SS_to_send = False
			self.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_SET_SS, self.params_SS))
		self.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_SET_OUTPUT, self.params_output))
		if self.query_panel.query_to_update:
			self.query_panel.OnLoadFromGui(None)
		if len(self.query) == 0:
			wx.MessageBox('Empty or Invalid query. Perhaps the wrong type (list/pairs) has been selected?', 'Invalid Query', 
        wx.OK | wx.ICON_ERROR)
			self._unfreeze()
			return
		self.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_SET_QUERY, self.params_query, self.query))

		#self.gui2ssprocess_queue.put((WorkProcess.CMD_STATUS, None))
		#print self.params_query
		#print self.params_output
		#print self.params_SS
		#print self.params_GO
		#print self.params_AC
		#print self.query

		self.query_panel.controls_start_button.Disable()
		self.gui2ssprocess_queue.put((WorkProcess.CMD_START, None))
#









	def stop(self):
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: stop()"
		self.query_panel.controls_log.AppendText("Stopping computation...\n")
		self.gui2ssprocess_queue.put((WorkProcess.CMD_STOP, None))
		self.query_panel.controls_start_button.Disable()
#








#############################################################################################
#############################################################################################
#################           GUI EVENTS HANDLERS          ####################################
#############################################################################################
#############################################################################################







	def activateGoCmd(self): # should not be called if status is 2
		if DEBUG_LEVEL>0:
			print "fastSemSimGui: activateGoCmd()"
		if self.status == -1:
			return
		if self.debug or (self.GO_status and self.AC_status and self.query_status and self.outputctrl_ok and self.operation_ok):
			self.SetStatus(1) # set to 0
		else:
			self.SetStatus(0)
			
	#def SetStatus(self, status):
		#self.status = status
		#return None ## To remove when complete
		#if self.show_pics:
			#if self.status == -1: # Fatal error
				#if not self.start_cmd is None:
					#self.start_cmd.Disable()
					#self.start_cmd.SetLabel("Start")
			#if self.status == 0: # Not running. Not Ready
				#if not self.start_cmd is None:
					#self.start_cmd.Disable()
					#self.start_cmd.SetLabel("Start")
			#if self.status == 1: # Not Running. Ready
				#if not self.start_cmd is None:
					#self.start_cmd.SetLabel("Start")
					#self.start_cmd.Enable()
			#if self.status == 2: # Running. (Completed)
				#if not self.start_cmd is None:
					#self.start_cmd.SetLabel("Stop")
					##self.statuspicture.SetBitmap(wx.Bitmap(self.s0_pic))
					#self.start_cmd.Enable()
#




##----------------------------------------------------------------------------------------------------------------------------------
	##'''
	##Menubar
	##'''
	#def InitMenu(self):
		#self.ID_NEW=wx.NewId()
		#self.ID_OPEN=wx.NewId()
		#self.ID_QUIT=wx.NewId()
		#self.ID_SAVE=wx.NewId()
		#self.ID_SAVE_AS=wx.NewId()
		#self.MenuBar = wx.MenuBar()
		#self.FileMenu = wx.Menu()
		#self.FileMenu.Append(self.ID_NEW, 'New', 'New')
		#self.FileMenu.Append(self.ID_OPEN, 'Open Configuration...', 'Open Configuration')
		#self.FileMenu.Append(self.ID_SAVE, 'Save Configuration', 'Save Configuration')
		#self.FileMenu.Append(self.ID_SAVE_AS, 'Save Configuration As...', 'Save Configuration As')
		#self.FileMenu.Append(self.ID_QUIT, 'Quit', 'Quit application')
		#self.MenuBar.Append(self.FileMenu,'&File')
		#self.SetMenuBar(self.MenuBar)
		#self.Bind(wx.EVT_MENU, self.OnMenuQuit, id=self.ID_QUIT)
		#self.Bind(wx.EVT_MENU, self.OnMenuOpen, id=self.ID_OPEN)
		#self.Bind(wx.EVT_MENU, self.OnMenuSave, id=self.ID_SAVE)
		#self.Bind(wx.EVT_MENU, self.OnMenuSaveAs, id=self.ID_SAVE_AS)
		#self.Bind(wx.EVT_MENU, self.OnMenuNew, id=self.ID_NEW)
		#return True
	
	#def OnMenuNew(self, event):
		#print "OnMenuNew still to be implemented."

	#def OnMenuQuit(self, event):
		#print "OnMenuQuit still to be implemented correctly."
		#sys.exit()

	#def OnMenuOpen(self, event):
		#print "OnMenuOpen still to be implemented."
		#dialog = wx.FileDialog(None, style = wx.OPEN)
		#if dialog.ShowModal() == wx.ID_OK:
			#print 'Loading: ', dialog.GetPath()
			#if not self.config_file == None:
				#self.OnMenuNew(None)
			#self.config_file = dialog.GetPath()
			#self.loadConfigGui = LoadConfigGui(self)

	#def OnMenuSave(self, event):
		#if self.config_file == None:
			#return self.OnMenuSaveAs(event)
		#self.saveConfigGui = SaveConfigGui(self)

	#def OnMenuSaveAs(self, event):
		#dialog = wx.FileDialog(None, style = wx.SAVE|wx.OVERWRITE_PROMPT)
		#if dialog.ShowModal() == wx.ID_OK:
			#self.config_file = dialog.GetPath()
			#self.OnMenuSave(event)
#






# routines to manage flags, status, and variables in general

	#def SetGoOk(self, status):
		#self.GO_status = status
		#if self.show_pics:
			#if self.GO_status:
				#self.StatusGui.GO_status_pic.SetBitmap(wx.Bitmap(self.Ok_pic))
			#else:
				#self.StatusGui.GO_status_pic.SetBitmap(wx.Bitmap(self.Warning_pic))
		#self.activateGoCmd()
		
	#def SetAcOk(self, status):
		#self.AC_status = status
		#if self.show_pics:
			#if self.AC_status:
				#self.StatusGui.AC_status_pic.SetBitmap(wx.Bitmap(self.Ok_pic))
			#else:
				#self.StatusGui.AC_status_pic.SetBitmap(wx.Bitmap(self.Warning_pic))
		#self.activateGoCmd()
		
	#def SetQueryOk(self, status):
		#self.query_status = status
		#if self.show_pics:
			#if self.query_status:
				#self.StatusGui.query_status_pic.SetBitmap(wx.Bitmap(self.Ok_pic))
			#else:
				#self.StatusGui.query_status_pic.SetBitmap(wx.Bitmap(self.Warning_pic))
		#self.activateGoCmd()
		
	#def SetOutputCtrlOk(self, status):
		#self.outputctrl_ok = status
		#if self.show_pics:
			#if self.outputctrl_ok:
				#self.StatusGui.output_status_pic.SetBitmap(wx.Bitmap(self.Ok_pic))
			#else:
				#self.StatusGui.output_status_pic.SetBitmap(wx.Bitmap(self.Warning_pic))
		#self.activateGoCmd()

	#def SetOperationOk(self, status):
		#self.operation_ok = status
		#if self.show_pics:
			#if self.operation_ok:
				#self.StatusGui.operation_status_pic.SetBitmap(wx.Bitmap(self.Ok_pic))
			#else:
				#self.StatusGui.operation_status_pic.SetBitmap(wx.Bitmap(self.Warning_pic))
		#self.activateGoCmd()
##


#----------------------------------------------------------------------------------------------------------------#








	CustomEvent_Generic, EVT_CUSTOM_GENERIC = wx.lib.newevent.NewEvent()
	CustomEvent_Status, EVT_CUSTOM_STATUS = wx.lib.newevent.NewEvent()
	CustomEvent_Start, EVT_CUSTOM_START = wx.lib.newevent.NewEvent()
	CustomEvent_Stop, EVT_CUSTOM_STOP = wx.lib.newevent.NewEvent()
	CustomEvent_Output, EVT_CUSTOM_OUTPUT = wx.lib.newevent.NewEvent()
	CustomEvent_Set, EVT_CUSTOM_SET = wx.lib.newevent.NewEvent()
	CustomEvent_Get, EVT_CUSTOM_GET = wx.lib.newevent.NewEvent()
	CustomEvent_LoadAC, EVT_CUSTOM_LOAD_AC = wx.lib.newevent.NewEvent()
	CustomEvent_LoadGO, EVT_CUSTOM_LOAD_GO = wx.lib.newevent.NewEvent()
	CustomEvent_SetSS, EVT_CUSTOM_SET_SS = wx.lib.newevent.NewEvent()
	CustomEvent_SetOutput, EVT_CUSTOM_SET_OUTPUT = wx.lib.newevent.NewEvent()
	CustomEvent_SetQuery, EVT_CUSTOM_SET_QUERY = wx.lib.newevent.NewEvent()
	CustomEvent_Kill, EVT_CUSTOM_KILL = wx.lib.newevent.NewEvent()







#----------------------------------------------------------------------------------------------------------------#
##################################################################################################################
##################################################################################################################
#################          COMMUNICATION CLASS - THREAD                    #######################################
##################################################################################################################
##################################################################################################################


class CommunicationThread(threading.Thread):

	def __init__(self, gui):
		threading.Thread.__init__ ( self)
		if DEBUG_LEVEL>2:
			print "CommunicationThread: init()"
		self.gui = gui
		self.uc = 0
		self.callbacks = {}
		self.lock = multiprocessing.Lock()
		
		self.connection = {}
		self.connection[WorkProcess.CMD_STATUS] = self.gui.CustomEvent_Status
		self.connection[WorkProcess.CMD_START] = self.gui.CustomEvent_Start
		self.connection[WorkProcess.CMD_STOP] = self.gui.CustomEvent_Stop
		self.connection[WorkProcess.CMD_OUTPUT] = self.gui.CustomEvent_Output
		self.connection[WorkProcess.CMD_SET] = self.gui.CustomEvent_Set
		self.connection[WorkProcess.CMD_GET] = self.gui.CustomEvent_Get
		self.connection[WorkProcess.CMD_LOAD_AC] = self.gui.CustomEvent_LoadAC
		self.connection[WorkProcess.CMD_LOAD_GO] = self.gui.CustomEvent_LoadGO
		self.connection[WorkProcess.CMD_SET_SS] = self.gui.CustomEvent_SetSS
		self.connection[WorkProcess.CMD_SET_OUTPUT] = self.gui.CustomEvent_SetOutput
		self.connection[WorkProcess.CMD_SET_QUERY] = self.gui.CustomEvent_SetQuery
		self.connection[WorkProcess.CMD_KILL] = self.gui.CustomEvent_Kill
#




	def _lock(self):
		self.lock.acquire()
#





	def _unlock(self):
		self.lock.release()
#



	def run(self):
		if DEBUG_LEVEL>2:
			print "CommunicationThread: run(). (Should appear once.)"
		while True:
			data = self.gui.ssprocess2gui_queue.get(True)
			if DEBUG_LEVEL>3:
				print "CommunicationThread: run(). message received from computation process"
			if data[0]==-1: # NOTE used to kill the thread
				return
			if data == None:
				return
			if data[0] in self.connection:
				ne = self.connection[data[0]](data=data)
			else:
				ne = CustomEvent_Generic(data=data)
			wx.PostEvent(self.gui.GetEventHandler(), ne) # WARNING is it safe tu run it without a lock?
#





	def register_callback(self, what, func):
		if DEBUG_LEVEL>2:
			print "CommunicationThread: register_callback()"
		
		self._lock()
		self.gui.Bind(what, func)
		self.callbacks[self.uc] = (what, func)
		self.uc += 1
		self._unlock()
		return self.uc - 1
#





	def unregister_callback(self, func_id):
		if func_id not in self.callbacks:
			if DEBUG_LEVEL>2:
				print "CommunicationThread: unregister_callback() error"
		else:
			self._lock()
			#print "A: " + str(self.callbacks[func_id][0])
			#print "B: " + str(self.callbacks[func_id][1])
			self.gui.Unbind(event=self.callbacks[func_id][0], handler=self.callbacks[func_id][1])
			del self.callbacks[func_id]
			self._unlock()
#



##################################################################################
#### start the main gui. Should be called from outside. Why? ####
##################################################################################

def _start():
	multiprocessing.freeze_support()
	app = wx.App()
	window = fastSemSimGui(None)
	window.Centre()
	window.Show()  
	app.MainLoop()
#

if __name__=='__main__':
	_start()
#
