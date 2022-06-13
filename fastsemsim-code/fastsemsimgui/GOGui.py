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
from fastsemsim.Ontology import GeneOntology
import WorkProcess
import os


DEBUG_LEVEL = 0 # NOTE 0: no debug output. 1:level 1 debug output. 2: verbose debug output. In any case it does not affect the performance of the calculation step

class GOPanel(wx.Panel):
	
	node_info = None
	edge_info = None
	
	def __init__( self, real_parent, parent, id, pos, size, style):
		wx.Panel.__init__ ( self, parent, id, pos, size, style)
		
		if DEBUG_LEVEL>1:
			print "GOPanel: init()"

		self.GO_panel = self # temporary workaround
		self.real_parent = real_parent

		

		GO_panel_sizer_1 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer3 = wx.StaticBoxSizer( wx.StaticBox( self.GO_panel, wx.ID_ANY, u"Source" ), wx.HORIZONTAL )
		
		self.GO_source_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"No file loaded", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		self.GO_source_label.Wrap( -1 )
		sbSizer3.Add( self.GO_source_label, 1, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		self.GO_load_button = wx.Button( self.GO_panel, wx.ID_ANY, u"Load ...", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_load_button.SetDefault() 
		sbSizer3.Add( self.GO_load_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		GO_panel_sizer_1.Add( sbSizer3, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer5 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer5 = wx.StaticBoxSizer( wx.StaticBox( self.GO_panel, wx.ID_ANY, u"Statistics" ), wx.HORIZONTAL )
		
		gSizer2 = wx.GridSizer( 6, 2, 0, 0 )
		
		self.m_staticText2 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"Terms", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText2.Wrap( -1 )
		gSizer2.Add( self.m_staticText2, 0, wx.ALL, 5 )
		
		self.GO_terms_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_terms_label.Wrap( -1 )
		gSizer2.Add( self.GO_terms_label, 0, wx.ALL, 5 )
		
		self.m_staticText4 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"Edges", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText4.Wrap( -1 )
		gSizer2.Add( self.m_staticText4, 0, wx.ALL, 5 )
		
		self.GO_edges_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_edges_label.Wrap( -1 )
		gSizer2.Add( self.GO_edges_label, 0, wx.ALL, 5 )
		
		self.m_staticText6 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"Categories", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText6.Wrap( -1 )
		gSizer2.Add( self.m_staticText6, 0, wx.ALL, 5 )
		
		self.GO_categories_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_categories_label.Wrap( -1 )
		gSizer2.Add( self.GO_categories_label, 0, wx.ALL, 5 )
		
		sbSizer5.Add( gSizer2, 1, wx.ALIGN_CENTER|wx.EXPAND|wx.SHAPED, 2 )
		
		bSizer5.Add( sbSizer5, 0, wx.ALL|wx.EXPAND, 5 )
		
		sbSizer6 = wx.StaticBoxSizer( wx.StaticBox( self.GO_panel, wx.ID_ANY, u"Relationships" ), wx.VERTICAL )
		
		gSizer22 = wx.GridSizer( 6, 2, 0, 0 )
		
		self.m_staticText24 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"is_a", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText24.Wrap( -1 )
		gSizer22.Add( self.m_staticText24, 0, wx.ALL, 5 )
		
		self.GO_isa_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_isa_label.Wrap( -1 )
		gSizer22.Add( self.GO_isa_label, 0, wx.ALL, 5 )
		
		self.m_staticText42 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"part_of", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText42.Wrap( -1 )
		gSizer22.Add( self.m_staticText42, 0, wx.ALL, 5 )
		
		self.GO_partof_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_partof_label.Wrap( -1 )
		gSizer22.Add( self.GO_partof_label, 0, wx.ALL, 5 )
		
		self.m_staticText62 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"has_part", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText62.Wrap( -1 )
		gSizer22.Add( self.m_staticText62, 0, wx.ALL, 5 )
		
		self.GO_haspart_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_haspart_label.Wrap( -1 )
		gSizer22.Add( self.GO_haspart_label, 0, wx.ALL, 5 )
		
		self.m_staticText59 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"regulates", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText59.Wrap( -1 )
		gSizer22.Add( self.m_staticText59, 0, wx.ALL, 5 )
		
		self.GO_regulates_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_regulates_label.Wrap( -1 )
		gSizer22.Add( self.GO_regulates_label, 0, wx.ALL, 5 )
		
		sbSizer6.Add( gSizer22, 1, wx.EXPAND, 5 )
		
		bSizer5.Add( sbSizer6, 0, wx.ALL, 5 )
		
		GO_panel_sizer_1.Add( bSizer5, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.GO_status_label = wx.StaticText( self.GO_panel, wx.ID_ANY, "", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_status_label.Wrap( -1 )
		self.GO_status_label.SetFont( wx.Font( wx.NORMAL_FONT.GetPointSize(), 70, 90, 92, False, wx.EmptyString ) )
		
		GO_panel_sizer_1.Add( self.GO_status_label, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_TOP|wx.ALL, 10 )
		
		self.GO_panel.SetSizer( GO_panel_sizer_1 )
		self.GO_panel.Layout()
		GO_panel_sizer_1.Fit( self.GO_panel )
		
		
		#GO_panel_sizer_1 = wx.BoxSizer( wx.VERTICAL )
		#sbSizer3 = wx.StaticBoxSizer( wx.StaticBox( self.GO_panel, wx.ID_ANY, u"Source" ), wx.HORIZONTAL )
		#self.GO_source_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"No file loaded", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		#self.GO_source_label.Wrap( -1 )
		#sbSizer3.Add( self.GO_source_label, 1, wx.ALIGN_CENTER|wx.ALL, 5 )
		#self.GO_load_button = wx.Button( self.GO_panel, wx.ID_ANY, u"Load ...", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.GO_load_button.SetDefault() 
		#sbSizer3.Add( self.GO_load_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		#GO_panel_sizer_1.Add( sbSizer3, 0, wx.ALL|wx.EXPAND, 5 )
		#bSizer5 = wx.BoxSizer( wx.HORIZONTAL )
		#sbSizer5 = wx.StaticBoxSizer( wx.StaticBox( self.GO_panel, wx.ID_ANY, u"Statistics" ), wx.HORIZONTAL )
		#gSizer2 = wx.GridSizer( 6, 2, 0, 0 )
		#self.m_staticText2 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"Terms", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.m_staticText2.Wrap( -1 )
		#gSizer2.Add( self.m_staticText2, 0, wx.ALL, 5 )
		#self.GO_terms_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.GO_terms_label.Wrap( -1 )
		#gSizer2.Add( self.GO_terms_label, 0, wx.ALL, 5 )
		#self.m_staticText4 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"Edges", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.m_staticText4.Wrap( -1 )
		#gSizer2.Add( self.m_staticText4, 0, wx.ALL, 5 )
		#self.GO_edges_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.GO_edges_label.Wrap( -1 )
		#gSizer2.Add( self.GO_edges_label, 0, wx.ALL, 5 )
		##self.m_staticText6 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"Categories", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.m_staticText6.Wrap( -1 )
		##gSizer2.Add( self.m_staticText6, 0, wx.ALL, 5 )
		##self.GO_categories_label = wx.StaticText( self.GO_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.GO_categories_label.Wrap( -1 )
		##gSizer2.Add( self.GO_categories_label, 0, wx.ALL, 5 )
		#sbSizer5.Add( gSizer2, 1, wx.ALIGN_CENTER|wx.SHAPED, 2 )
		#bSizer5.Add( sbSizer5, 0, wx.ALL, 5 )
		#sbSizer6 = wx.StaticBoxSizer( wx.StaticBox( self.GO_panel, wx.ID_ANY, u"Additional information" ), wx.VERTICAL )
		#bSizer51 = wx.BoxSizer( wx.VERTICAL )
		#self.m_staticText8 = wx.StaticText( self.GO_panel, wx.ID_ANY, "", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.m_staticText8.Wrap( -1 )
		#bSizer51.Add( self.m_staticText8, 0, wx.ALL, 5 )
		#sbSizer6.Add( bSizer51, 1, wx.EXPAND, 5 )
		#bSizer5.Add( sbSizer6, 0, wx.ALL, 5 )
		#GO_panel_sizer_1.Add( bSizer5, 0, wx.EXPAND, 5 )
		#bSizer4 = wx.BoxSizer( wx.HORIZONTAL )
		#GO_panel_sizer_1.Add( bSizer4, 0, wx.ALIGN_RIGHT, 5 )
		#self.GO_panel.SetSizer( GO_panel_sizer_1 )
		#self.GO_panel.Layout()
		#GO_panel_sizer_1.Fit( self.GO_panel )
		
		self.go_load_gui = GO_load_gui(self, self.real_parent)
		
		self.Bind(wx.EVT_BUTTON, self.OnLoadButton, id=self.GO_load_button.GetId())
		self.GO_load_outcome_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_LOAD_GO, self.OnLoadDone)
		#self.GO_params_handle = None
		self.GO_params_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_GET, self.OnGetParams)
		self._reset()
#






	def _reset(self):
		self._update()
#







	def _freeze(self):
		if DEBUG_LEVEL>2:
			print "GOPanel: _freeze()"
		self.GO_load_button.Disable()
#






	def _unfreeze(self):
		if DEBUG_LEVEL>2:
			print "GOPanel: _unfreeze()"
		self.GO_load_button.Enable()
#





	def OnLoadButton(self, event):
		if DEBUG_LEVEL>2:
			print "GOPanel: OnLoadButton()"
		self.go_load_gui.ShowModal()
#





	def _update(self):
		if DEBUG_LEVEL>1:
			print "GOPanel: _update()"
		#print "DEBUG 3 " + str(self.node_info)
		if self.real_parent.GO_status:
			if not self.real_parent.params_GO['filename']==None:
				self.GO_source_label.SetLabel(os.path.basename(self.real_parent.params_GO['filename']))
			else:
				self.GO_source_label.SetLabel("[built-in version]")

			self.GO_terms_label.SetLabel(str(self.node_info))
			self.GO_edges_label.SetLabel(str(self.edge_info))
			self.GO_categories_label.SetLabel(str(3))
			self.GO_status_label.SetLabel("Gene Ontology correcty loaded")
			self.GO_isa_label.SetLabel("considered")
			self.GO_haspart_label.SetLabel("considered")
			self.GO_regulates_label.SetLabel("considered")
			self.GO_partof_label.SetLabel("considered")
			if "ignore" in self.real_parent.params_GO:
				if "is_a" in self.real_parent.params_GO['ignore'] and self.real_parent.params_GO['ignore']['is_a']:
					self.GO_isa_label.SetLabel("ignored")
				if "has_part" in self.real_parent.params_GO['ignore'] and self.real_parent.params_GO['ignore']['has_part']:
					self.GO_haspart_label.SetLabel("ignored")
				if "regulates" in self.real_parent.params_GO['ignore'] and self.real_parent.params_GO['ignore']['regulates']:
					self.GO_regulates_label.SetLabel("ignored")
				if "part_of" in self.real_parent.params_GO['ignore'] and self.real_parent.params_GO['ignore']['part_of']:
					self.GO_partof_label.SetLabel("ignored")
		else:
			self.GO_source_label.SetLabel("-")
			self.GO_terms_label.SetLabel("-")
			self.GO_edges_label.SetLabel("-")
			self.GO_categories_label.SetLabel("-")
			self.GO_isa_label.SetLabel("-")
			self.GO_haspart_label.SetLabel("-")
			self.GO_regulates_label.SetLabel("-")
			self.GO_partof_label.SetLabel("-")
			self.GO_status_label.SetLabel("No Gene Ontology currently loaded")
			#self.m_staticText8.SetLabel("No Gene Ontology is currently loaded.")
#




	def OnLoadDone(self, event):
		if DEBUG_LEVEL>0:
			print "GOPanel: OnLoadDone()"
		data = event.data

		if data[0] == WorkProcess.CMD_LOAD_GO:
			if data[1] == WorkProcess.ANSWER_PROCESSED:
				if data[2] == WorkProcess. RESULT_OK:
					if DEBUG_LEVEL>1:
						print "GO load outcome: Load successful."
					self.node_info = data[3]
					self.edge_info = data[4]
					self.go_load_gui.Hide()
					
				else:
					if DEBUG_LEVEL>2:
						print "GO load outcome: Load Fail."
					self.node_info = None
					self.edge_info = None
					wx.MessageBox('Gene Ontology not loaded. Perhaps the wrong file has been selected?', 'Invalid Gene Ontology file', 
        wx.OK | wx.ICON_ERROR)
			elif data[1] == WorkProcess.ANSWER_PROCESSING:
				if DEBUG_LEVEL>2:
					print "GO load outcome: Load in progress."
				return
			elif data[1] == WorkProcess.ANSWER_IGNORED:
				if DEBUG_LEVEL>2:
					print "GO load outcome: Load request ignored."
			else:
				if DEBUG_LEVEL>2:
					print "GO load outcome: Unknown answer."
		#self.unfreeze()
		self.go_load_gui.unfreeze()
		#if self.GO_params_handle == None:
			#self.GO_params_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_GET, self.OnGetParams)
		self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_GET, WorkProcess.CMD_GET_PARAMS, WorkProcess.CMD_GET_PARAMS_GO))
		self.real_parent.update()
		event.Skip()
#






	def OnGetParams(self, event):
		if DEBUG_LEVEL>0:
			print "GOPanel: OnGetParams()"
		data = event.data

		if data[0] == WorkProcess.CMD_GET:
			if data[1] == WorkProcess.CMD_GET_PARAMS:
				if data[2] == WorkProcess. CMD_GET_PARAMS_GO:
					data = data[3]
					#print data
					self.real_parent.params_GO = data
					self._update()
					#self.real_parent.params_GO['filename'] = self.param_filename
					#self.real_parent.params_GO['ignore']['has_part'] = self.param_GO_ignore_haspart
					#self.real_parent.params_GO['ignore']['regulates'] = self.param_GO_ignore_regulates
					#self.real_parent.params_GO['ignore']['is_a'] = self.param_GO_ignore_isa
					#self.real_parent.params_GO['ignore']['part_of'] = self.param_GO_ignore_partof
					#if not self.GO_params_handle == None:
						#self.real_parent.communication_thread.unregister_callback(self.GO_params_handle)
						#self.GO_params_handle = None
					#self.real_parent.update()
		event.Skip()
#








	#def _data_to_main(self):
		#if DEBUG_LEVEL>1:
			#print "GO_load_gui: _data_to_main()"

		#self.real_parent.GO_status = self.status
		#self.parent.edge_info = self.edge_info
		#self.parent.node_info = self.node_info
		#print "DEBUG 2 " + str(self.parent.node_info)
##



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# GO load window #


class GO_load_gui ( wx.Dialog ):
	
	param_GO_ignore_haspart = False
	param_GO_ignore_isa = False
	param_GO_ignore_partof = False
	param_GO_ignore_regulates = False
	param_filename = None
	
	def __init__( self, parent , real_parent):
		wx.Dialog.__init__ ( self, parent, id = wx.ID_ANY, title = u"Load Gene Ontology", pos = wx.DefaultPosition, size = wx.Size( 650,410 ), style = wx.DEFAULT_DIALOG_STYLE )
		
		if DEBUG_LEVEL>2:
			print "GO_load_gui: init()"
		self.parent = parent
		self.real_parent = real_parent
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		bSizer3 = wx.BoxSizer( wx.VERTICAL )
		bSizer22 = wx.BoxSizer( wx.HORIZONTAL )
		sbSizer3 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Source file" ), wx.HORIZONTAL )
		self.GO_load_source_label = wx.StaticText( self, wx.ID_ANY, u"No file selected. The built-in version will be loaded.", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_LEFT)
		self.GO_load_source_label.Wrap( -1 )
		sbSizer3.Add( self.GO_load_source_label, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL|wx.EXPAND, 5 )
		bSizer22.Add( sbSizer3, 1, wx.LEFT|wx.RIGHT, 5 )
		bSizer20 = wx.BoxSizer( wx.HORIZONTAL )
		self.GO_load_select_button = wx.Button( self, wx.ID_ANY, u"Select file...", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer20.Add( self.GO_load_select_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		self.GO_load_default_button = wx.Button( self, wx.ID_ANY, u"Select built-in GO", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer20.Add( self.GO_load_default_button, 0, wx.ALIGN_CENTER, 5 )
		bSizer22.Add( bSizer20, 0, wx.ALIGN_RIGHT|wx.RIGHT|wx.TOP, 5 )
		bSizer3.Add( bSizer22, 0, wx.ALIGN_CENTER|wx.EXPAND|wx.TOP, 5 )
		bSizer5 = wx.BoxSizer( wx.VERTICAL )
		bSizer23 = wx.BoxSizer( wx.VERTICAL )
		sbSizer5 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Options" ), wx.HORIZONTAL )
		fgSizer1 = wx.FlexGridSizer( 4, 2, 0, 0 )
		fgSizer1.SetFlexibleDirection( wx.BOTH )
		fgSizer1.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		self.GO_load_ignore_haspart_check = wx.CheckBox( self, wx.ID_ANY, u"Ignore has_part", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer1.Add( self.GO_load_ignore_haspart_check, 1, wx.ALL, 5 )
		self.m_staticText59 = wx.StaticText( self, wx.ID_ANY, u"Check this to ignore 'has_part' relationships", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText59.Wrap( -1 )
		fgSizer1.Add( self.m_staticText59, 1, wx.ALIGN_RIGHT|wx.ALL, 5 )
		self.GO_load_ignore_regulates_check = wx.CheckBox( self, wx.ID_ANY, u"Ignore regulates", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer1.Add( self.GO_load_ignore_regulates_check, 1, wx.ALL, 5 )
		self.m_staticText591 = wx.StaticText( self, wx.ID_ANY, u"Check this to ignore 'regulates' relationships", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText591.Wrap( -1 )
		fgSizer1.Add( self.m_staticText591, 1, wx.ALIGN_RIGHT|wx.ALL, 5 )
		
		self.GO_load_ignore_partof_check = wx.CheckBox( self, wx.ID_ANY, u"Ignore part_of", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer1.Add( self.GO_load_ignore_partof_check, 1, wx.ALL, 5 )
		self.m_staticText59 = wx.StaticText( self, wx.ID_ANY, u"Check this to ignore 'part_of' relationships", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText59.Wrap( -1 )
		fgSizer1.Add( self.m_staticText59, 1, wx.ALIGN_RIGHT|wx.ALL, 5 )
		
		self.GO_load_ignore_isa_check = wx.CheckBox( self, wx.ID_ANY, u"Ignore is_a", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer1.Add( self.GO_load_ignore_isa_check, 1, wx.ALL, 5 )
		self.m_staticText59 = wx.StaticText( self, wx.ID_ANY, u"Check this to ignore 'is_a' relationships", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText59.Wrap( -1 )
		fgSizer1.Add( self.m_staticText59, 1, wx.ALIGN_RIGHT|wx.ALL, 5 )
		
		
		sbSizer5.Add( fgSizer1, 1, wx.EXPAND, 5 )
		bSizer23.Add( sbSizer5, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		bSizer5.Add( bSizer23, 1, wx.EXPAND, 5 )
		bSizer3.Add( bSizer5, 0, wx.ALL|wx.EXPAND, 5 )
		bSizer4 = wx.BoxSizer( wx.HORIZONTAL )
		self.GO_load_load_button = wx.Button( self, wx.ID_ANY, u"Load", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_load_load_button.SetDefault() 
		bSizer4.Add( self.GO_load_load_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		self.GO_load_cancel_button = wx.Button( self, wx.ID_ANY, u"Cancel", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer4.Add( self.GO_load_cancel_button, 0, wx.ALL, 5 )
		bSizer3.Add( bSizer4, 0, wx.ALIGN_CENTER, 5 )
		sbSizer18 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Help" ), wx.VERTICAL )
		self.m_staticText31 = wx.StaticText( self, wx.ID_ANY, u"The Gene Ontology file must be encoded in obo-xml format. It can also be gzipped. Please refer to www.geneontology.org to retrieve the most updated Gene Ontology. If not provided, the Gene Ontology provided with fastSemSim will be loaded.", wx.DefaultPosition, wx.Size( 500,150 ), 0 )
		self.m_staticText31.Wrap( -1 )
		self.m_staticText31.SetMaxSize( wx.Size( 500,-1 ) )
		sbSizer18.Add( self.m_staticText31, 0, wx.ALL, 5 )
		bSizer3.Add( sbSizer18, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		self.SetSizer( bSizer3 )
		self.Layout()
		self.Centre( wx.BOTH )
# Bind events to controls
		self.Bind(wx.EVT_BUTTON, self.OnFileBrowse, id=self.GO_load_select_button.GetId())
		self.Bind(wx.EVT_BUTTON, self.OnCancel, id=self.GO_load_cancel_button.GetId())
		self.Bind(wx.EVT_BUTTON, self.OnLoad, id=self.GO_load_load_button.GetId())
		self.Bind(wx.EVT_BUTTON, self.OnDefaultButton, id=self.GO_load_default_button.GetId())
		self.Bind(wx.EVT_CHECKBOX, self.ignore_what, id=self.GO_load_ignore_haspart_check.GetId())
		self.Bind(wx.EVT_CHECKBOX, self.ignore_what, id=self.GO_load_ignore_regulates_check.GetId())
		self.Bind(wx.EVT_CHECKBOX, self.ignore_what, id=self.GO_load_ignore_partof_check.GetId())
		self.Bind(wx.EVT_CHECKBOX, self.ignore_what, id=self.GO_load_ignore_isa_check.GetId())
		
		self.OnDefaultButton(None)
#



	def ignore_what(self, event):
		self.param_GO_ignore_haspart = self.GO_load_ignore_haspart_check.GetValue()
		self.param_GO_ignore_regulates = self.GO_load_ignore_regulates_check.GetValue()
		self.param_GO_ignore_isa = self.GO_load_ignore_isa_check.GetValue()
		self.param_GO_ignore_partof = self.GO_load_ignore_partof_check.GetValue()
#




	def OnCancel(self, event):
		self.Hide()
#



	def OnDefaultButton(self, event):
		self._set_file_name(self.real_parent.programdirectory + "/../fastSemSim/data/GO_2012-10-22.obo-xml.gz")
		#self._set_file_name("data/GO_2012-10-22.obo-xml.gz")
#



	def OnFileBrowse(self, event):
		dialog = wx.FileDialog(None, style = wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			self._set_file_name(dialog.GetPath())
#


	def _set_file_name(self, fn):
		self.param_filename = fn
		if fn == None:
			self.GO_load_source_label.SetLabel(u"No GO selected.")
		else:
			self.GO_load_source_label.SetLabel(os.path.basename(self.param_filename))
#


	def freeze(self):
		self.GO_load_cancel_button.Disable()
		self.GO_load_load_button.Disable()
		self.GO_load_default_button.Disable()
		self.GO_load_select_button.Disable()
		self.GO_load_ignore_haspart_check.Disable()
		self.GO_load_ignore_regulates_check.Disable()
		self.GO_load_ignore_isa_check.Disable()
		self.GO_load_ignore_partof_check.Disable()
#

	def unfreeze(self):
		self.GO_load_cancel_button.Enable()
		self.GO_load_load_button.Enable()
		self.GO_load_default_button.Enable()
		self.GO_load_select_button.Enable()
		self.GO_load_ignore_haspart_check.Enable()
		self.GO_load_ignore_regulates_check.Enable()
		self.GO_load_ignore_isa_check.Enable()
		self.GO_load_ignore_partof_check.Enable()
#

	def OnLoad(self, data):
		if DEBUG_LEVEL>1:
			print "GO_load_gui: OnLoad()"
		self.freeze()
		
		temp = {}
		temp['filename'] = self.param_filename
		temp['type'] = 'obo-xml'
		temp['ignore'] = {"has_part":self.param_GO_ignore_haspart, "regulates":self.param_GO_ignore_regulates, 'part_of':self.param_GO_ignore_partof, 'is_a':self.param_GO_ignore_isa}
		
		self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_LOAD_GO, temp))
#
