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
along with fastSemSim.  If not, see <http://www.gnu.org/licenses/>.
'''

import wx
import time
import WorkProcess

DEBUG_LEVEL = 0

class ConfigFileParser():
	tags = {}

	separator = '\t'
	
	def __init__(self, filename):
		self.filename = str(filename)

	def load(self):
		if type(self.filename) == str:
			stream = open(self.filename,'r')
		else:
			stream = self.filename
			
		lines_counter = 0
		
		for line in stream:
			line = str(line.rstrip('\n'))
			line = line.rstrip('\r')
			line = line.split(self.separator) # should be there at least 3 tokens
			if len(line) < 2:
				continue
			
			cur = self.tags
			precur = None

			for j in range(0,len(line)):
				if j == (len(line) - 1):
					tow = line[j]
					if len(tow) > 0 and not tow[0] == "\"":
						try:
							tow = int(tow)
						except Exception:
							try:
								tow = float(tow)
							except Exception:
								pass
					else:
						tow = tow[1:len(tow)-1]
					precur[line[j-1]] = tow
				else:
					if not line[j] in cur:
						cur[line[j]] = {}
					precur = cur
					cur = cur[line[j]]
				lines_counter += 1
		if type(self.filename) == str:
			stream.close()
		return self.tags
#



	def save(self):
		if type(self.filename) == str:
			stream = open(self.filename,'w')
		else:
			stream = self.filename
		self.save_branch(self.tags, '', stream)
		if type(self.filename) == str:
			stream.close()
		return self.tags
#



	def save_branch(self, branch, prefix, h):
		for j in branch:
			if branch[j] == None:
				pass
			elif type(branch[j]) == dict:
				self.save_branch(branch[j], prefix + str(j) + self.separator, h)
			elif type(branch[j]) == str or type(branch[j]) == unicode:
				h.write(prefix + j + "\t\"" + str(branch[j]) + "\"\n")
			elif type(branch[j]) == bool:
				h.write(prefix + j + "\t" + str(int(branch[j])) + "\n")
			else:
				h.write(prefix + j + "\t" + str(branch[j]) + "\n")
#





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


class DummyEvent(wx.PyEvent):
	def __init__(self, data):
		wx.PyEvent.__init__(self)
		self.SetEventType(data)
		self.data = None







class LoadConfigGui(wx.Dialog):

	def __init__(self, parent):
		self.real_parent = parent
		super(LoadConfigGui, self).__init__(self.real_parent, title="Loading Configuration", size=(500,200))
		self.InitUI()
		self.Process()
	
	def InitUI(self):
		self.panel = wx.Panel(self, -1)
		self.mainbox = wx.BoxSizer(wx.VERTICAL)
		self.filebox = wx.BoxSizer(wx.HORIZONTAL)
		self.file_label = wx.StaticText(self.panel, label='Configuration file: ')
		self.filename_label = wx.StaticText(self.panel, label='')
		self.filebox.Add(self.file_label)
		self.filebox.Add(self.filename_label)
		self.mainbox.Add(self.filebox, flag=wx.EXPAND|wx.LEFT|wx.RIGHT|wx.TOP, border=10)
		
		self.go_label = wx.StaticText(self.panel, label='Load Gene Ontology...')
		self.ac_label = wx.StaticText(self.panel, label='Load Annotation Corpus...')
		self.ss_label = wx.StaticText(self.panel, label='Load Semantic Similarity parameters...')
		self.query_label = wx.StaticText(self.panel, label='Load Query...')
		self.output_label = wx.StaticText(self.panel, label='Load Output parameters...')
		self.go_label.SetFont(self.real_parent.font)
		self.ac_label.SetFont(self.real_parent.font)
		self.ss_label.SetFont(self.real_parent.font)
		self.query_label.SetFont(self.real_parent.font)
		self.output_label.SetFont(self.real_parent.font)
		
		self.go_label.Disable()
		self.ac_label.Disable()
		self.ss_label.Disable()
		self.query_label.Disable()
		self.output_label.Disable()
		
		self.go_pic = wx.StaticBitmap(self.panel)
		self.ac_pic = wx.StaticBitmap(self.panel)
		self.query_pic = wx.StaticBitmap(self.panel)
		self.ss_pic = wx.StaticBitmap(self.panel)
		self.output_pic = wx.StaticBitmap(self.panel)
		
		self.descbox = wx.FlexGridSizer(rows = 5, cols = 2, vgap = 6, hgap = 20)
		self.descbox.AddMany([(self.go_label), (self.go_pic), (self.ac_label), (self.ac_pic), (self.ss_label), (self.ss_pic), (self.query_label), (self.query_pic), (self.output_label), (self.output_pic)])
		
		self.mainbox.Add(self.descbox, flag=wx.EXPAND|wx.LEFT|wx.RIGHT|wx.TOP, border=10)
		#self.progress = wx.Gauge(self.panel, -1, 50, size=(250, 25))
		#self.mainbox.Add(self.progress, flag=wx.TOP, border = 5)
		self.panel.SetSizerAndFit(self.mainbox)
		
		#self.EVT_GO_ID = wx.NewId()
		#self.EVT_AC_ID = wx.NewId()
		#self.EVT_SS_ID = wx.NewId()
		self.EVT_QUERY_ID = wx.NewId()
		#self.EVT_OUTPUT_ID = wx.NewId()
		#self.EVT_DONE_ID = wx.NewId()
		#self.Connect(-1, -1, self.EVT_GO_ID, self.ProcessGO)
		#self.Connect(-1, -1, self.EVT_AC_ID, self.ProcessAC)
		#self.Connect(-1, -1, self.EVT_SS_ID, self.ProcessSS)
		self.Connect(-1, -1, self.EVT_QUERY_ID, self.ProcessQuery)
		#self.Connect(-1, -1, self.EVT_OUTPUT_ID, self.ProcessOutput)
		#self.Connect(-1, -1, self.EVT_DONE_ID, self.ProcessDone)
		
		self.Show()

	def Process(self):
		#load configuration
		if self.real_parent.config_file == None:
			return
		self.filename_label.SetLabel(self.real_parent.config_file)
		self.config = ConfigFileParser(self.real_parent.config_file)
		self.tags = self.config.load()
		#print self.tags
		self.ProcessGO()
#






	def ProcessGO(self):
		temp = self.tags['go']
		self.real_parent.params_GO = temp.copy()
		self.GO_load_outcome_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_LOAD_GO, self.OnProcessGODone)
		self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_LOAD_GO, temp))
		self.go_label.Enable()
#




	def OnProcessGODone(self, event):
		if DEBUG_LEVEL>0:
			print "ConfigGui: OnProcessGODone()"
		data = event.data

		if data[0] == WorkProcess.CMD_LOAD_GO:
			if data[1] == WorkProcess.ANSWER_PROCESSED:
				if data[2] == WorkProcess. RESULT_OK:
					self.go_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
				else:
					print "OnProcessGODone: Errore in loading!"
					self.go_pic.SetBitmap(wx.Bitmap(self.real_parent.Warning_pic))

				self.real_parent.communication_thread.unregister_callback(self.GO_load_outcome_handle)
				event.Skip()
				#event.Allow()
				self.ProcessAC()

			elif data[1] == WorkProcess.ANSWER_PROCESSING:
				pass
			elif data[1] == WorkProcess.ANSWER_IGNORED:
				if DEBUG_LEVEL>2:
					print "GO load outcome: Load request ignored."
				self.real_parent.communication_thread.unregister_callback(self.GO_load_outcome_handle)
				event.Skip()
				#event.Allow()
				self.ProcessAC()
			else:
				if DEBUG_LEVEL>2:
					print "GO load outcome: Unknown answer."
				self.real_parent.communication_thread.unregister_callback(self.GO_load_outcome_handle)
				event.Skip()
				#event.Allow()
				self.ProcessAC()
#








	def ProcessAC(self):
		self.ac_label.Enable()
		if not 'ac' in self.tags:
			self.ac_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
			#self.ProcessQuery()
			self.real_parent.update()
			wx.PostEvent(self.GetEventHandler(), DummyEvent(self.EVT_QUERY_ID))
			return
		temp = self.tags['ac']
		self.real_parent.params_AC = temp.copy()
		self.AC_load_outcome_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_LOAD_AC, self.OnProcessACDone)
		self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_LOAD_AC, temp))
#






	def OnProcessACDone(self, event):
		if DEBUG_LEVEL>0:
			print "ConfigGui: OnProcessACDone()"
		data = event.data

		if data[0] == WorkProcess.CMD_LOAD_AC:
			if data[1] == WorkProcess.ANSWER_PROCESSED:
				if data[2] == WorkProcess. RESULT_OK:
					self.ac_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
				else:
					self.ac_pic.SetBitmap(wx.Bitmap(self.real_parent.Warning_pic))

				self.real_parent.communication_thread.unregister_callback(self.AC_load_outcome_handle)
				event.Skip()
				#event.Allow()
				#self.ProcessQuery()
				self.real_parent.update()
				wx.PostEvent(self.GetEventHandler(), DummyEvent(self.EVT_QUERY_ID))

			elif data[1] == WorkProcess.ANSWER_PROCESSING:
				pass
			elif data[1] == WorkProcess.ANSWER_IGNORED:
				if DEBUG_LEVEL>2:
					print "AC load outcome: Load request ignored."
				self.real_parent.communication_thread.unregister_callback(self.AC_load_outcome_handle)
				event.Skip()
				#event.Allow()
				self.ProcessQuery()
				self.real_parent.update()
				wx.PostEvent(self.GetEventHandler(), DummyEvent(self.EVT_QUERY_ID))
				
			else:
				if DEBUG_LEVEL>2:
					print "AC load outcome: Unknown answer."
				self.real_parent.communication_thread.unregister_callback(self.AC_load_outcome_handle)
				event.Skip()
				#event.Allow()
				#self.ProcessQuery()
				self.real_parent.update()
				wx.PostEvent(self.GetEventHandler(), DummyEvent(self.EVT_QUERY_ID))
#






	def ProcessQuery(self, event):
		self.query_label.Enable()
		if not 'query' in self.tags:
			self.query_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
			self.ProcessSS()
			return
		self.real_parent.params_query = self.tags['query'].copy()
		self.query_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
		
		if self.real_parent.params_query['source'] == WorkProcess.QUERY_FROM_AC:
			self.real_parent.query_panel.OnLoadFromAC(None)
		elif self.real_parent.params_query['source'] == WorkProcess.QUERY_FROM_FILE:
			print "Batch load query from file still to be implemented"
		elif self.real_parent.params_query['source'] == WorkProcess.QUERY_FROM_GUI:
			self.real_parent.query_panel.OnReset(None)
			
		#self.real_parent.query_file = self.tags['query']['query_filename']
		#self.real_parent.query_type = int(self.tags['query']['query_format'])
		#self.real_parent.query_from = int(self.tags['query']['query_from'])
		#if not self.tags['query']['query_format'] == None:
			#self.real_parent.QueryGui.set_query_format(int(self.tags['query']['query_format']))
		#if not self.tags['query']['query_from'] == None:
			#self.real_parent.QueryGui.set_query_from(int(self.tags['query']['query_from']))
		self.ProcessSS()
#





	def ProcessSS(self):
		self.ss_label.Enable()
		if not 'ss' in self.tags:
			self.ss_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
			self.ProcessOutput()
			return
		
		self.real_parent.params_SS = self.tags['ss'].copy()
		#self.real_parent.selectedGO = self.tags['ss']['ss_go_tree']
		#self.real_parent.mixingstrategy = self.tags['ss']['ss_mixing_strategy']
		#self.real_parent.ssmeasure = self.tags['ss']['ss_measure']
		#if not self.tags['ss']['ss_measure'] == None:
			#self.real_parent.OperationGui.set_ss(self.tags['ss']['ss_measure'])
		#if not self.tags['ss']['ss_mixing_strategy'] == None:
			#self.real_parent.OperationGui.set_ms(self.tags['ss']['ss_mixing_strategy'])
		#if not self.tags['ss']['ss_go_tree'] == None:
			#self.real_parent.OperationGui.set_go(self.tags['ss']['ss_go_tree'])
		self.ss_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
		#self.ss_pic.SetBitmap(wx.Bitmap(self.real_parent.Warning_pic))
		self.ProcessOutput()
#






	def ProcessOutput(self):
		self.output_label.Enable()
		if not 'output' in self.tags:
			self.output_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
			self.ProcessDone()
			return
		self.real_parent.params_output = self.tags['output'].copy()
		#self.real_parent.output_type = int(self.tags['output']['output_to'])
		#self.real_parent.output_filename = self.tags['output']['output_filename']
		#if not self.tags['output']['output_to'] == None:
			#self.real_parent.OutputCtrlGui.set_output_type(int(self.tags['output']['output_to']))
		self.output_pic.SetBitmap(wx.Bitmap(self.real_parent.Ok_pic))
		#self.output_pic.SetBitmap(wx.Bitmap(self.real_parent.Warning_pic))
		self.ProcessDone()
#






	def ProcessDone(self):
		self.real_parent.update()
		#self.Quit()
		pass
#





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-



class SaveConfigGui(wx.Dialog):
	
	def __init__(self, parent):
		self.real_parent = parent
		super(SaveConfigGui, self).__init__(self.real_parent, title="Saving Configuration", size=(400,200))
		#self.InitUI()
		self.CollectData()
		self.SaveData()
		
	def CollectData(self):
		self.OnCollectData(None)
		#self.temp_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_STATUS, self.OnCollectData)
		#self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_LOAD_GO, self.param_filename, {"ignore_part_of":self.param_GO_ignore_haspart, "ignore_regulates":self.param_GO_ignore_regulates}))

	def OnCollectData(self, event):
		self.config = ConfigFileParser(self.real_parent.config_file)
		self.config.tags['go'] = self.real_parent.params_GO.copy()
		self.config.tags['ac'] = self.real_parent.params_AC.copy()
		self.config.tags['ss'] = self.real_parent.params_SS.copy()
		self.config.tags['query'] = self.real_parent.params_query.copy()
		self.config.tags['output'] = self.real_parent.params_output.copy()
#


	def SaveData(self):
		self.config.save()
#
