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

import wx
        
from fastsemsim.Ontology import GeneOntology
import WorkProcess
import os

from fastsemsim.Ontology import GeneOntology
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.Ontology import PlainAnnotationCorpus
from fastsemsim.Ontology import GAF2AnnotationCorpus

GAF2_FILE_TYPE = 'GAF 2.0 (GOA)'
GAF2_FILE_TYPE_REF = 'gaf-2.0'
PLAIN_FILE_TYPE = 'Plain'
PLAIN_FILE_TYPE_REF = 'plain'
AC_FILE_TYPES = { GAF2_FILE_TYPE : GAF2_FILE_TYPE_REF, PLAIN_FILE_TYPE : PLAIN_FILE_TYPE_REF }

DEBUG_LEVEL = 0

class ACPanel(wx.Panel):
	def __init__( self, real_parent, parent, id, pos, size, style):
		wx.Panel.__init__ ( self, parent, id, pos, size, style)
		self.AC_panel = self # temporary workaround
		self.real_parent = real_parent
		
		
		bSizer31 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer28 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer31 = wx.StaticBoxSizer( wx.StaticBox( self.AC_panel, wx.ID_ANY, u"Source file" ), wx.HORIZONTAL )
		
		self.AC_source_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.AC_source_label.Wrap( -1 )
		sbSizer31.Add( self.AC_source_label, 1, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		self.AC_load_button = wx.Button( self.AC_panel, wx.ID_ANY, u"Load ...", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.AC_load_button.SetDefault() 
		sbSizer31.Add( self.AC_load_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		bSizer28.Add( sbSizer31, 6, wx.ALL|wx.EXPAND, 5 )
		
		bSizer31.Add( bSizer28, 0, wx.EXPAND, 5 )
		
		bSizer52 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer51 = wx.StaticBoxSizer( wx.StaticBox( self.AC_panel, wx.ID_ANY, u"Statistics" ), wx.HORIZONTAL )
		
		gSizer21 = wx.GridSizer( 6, 2, 0, 0 )
		
		self.m_staticText21 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"Terms", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText21.Wrap( -1 )
		gSizer21.Add( self.m_staticText21, 0, wx.ALL, 5 )
		
		self.terms_number_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.terms_number_label.Wrap( -1 )
		gSizer21.Add( self.terms_number_label, 0, wx.ALL, 5 )
		
		self.m_staticText41 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"Objects", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText41.Wrap( -1 )
		gSizer21.Add( self.m_staticText41, 0, wx.ALL, 5 )
		
		self.objects_number_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.objects_number_label.Wrap( -1 )
		gSizer21.Add( self.objects_number_label, 0, wx.ALL, 5 )
		
		sbSizer51.Add( gSizer21, 1, wx.ALIGN_CENTER|wx.SHAPED, 2 )
		
		bSizer52.Add( sbSizer51, 0, wx.ALL|wx.EXPAND, 5 )
		
		sbSizer512 = wx.StaticBoxSizer( wx.StaticBox( self.AC_panel, wx.ID_ANY, u"Information" ), wx.HORIZONTAL )
		
		gSizer211 = wx.GridSizer( 6, 2, 0, 0 )
		
		self.m_staticText211 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"File Type", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText211.Wrap( -1 )
		gSizer211.Add( self.m_staticText211, 0, wx.ALL, 5 )
		
		self.file_type_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.file_type_label.Wrap( -1 )
		gSizer211.Add( self.file_type_label, 0, wx.ALL, 5 )
		
		self.m_staticText411 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"IEA annotations", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText411.Wrap( -1 )
		gSizer211.Add( self.m_staticText411, 0, wx.ALL, 5 )
		
		self.IEA_annotations_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.IEA_annotations_label.Wrap( -1 )
		gSizer211.Add( self.IEA_annotations_label, 0, wx.ALL, 5 )
		
		sbSizer512.Add( gSizer211, 1, wx.ALIGN_CENTER|wx.SHAPED, 2 )
		
		bSizer52.Add( sbSizer512, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer31.Add( bSizer52, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.AC_status_label = wx.StaticText( self.AC_panel, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		self.AC_status_label.Wrap( -1 )
		self.AC_status_label.SetFont( wx.Font( wx.NORMAL_FONT.GetPointSize(), 70, 90, 92, False, wx.EmptyString ) )
		
		bSizer31.Add( self.AC_status_label, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 10 )
		
		self.AC_panel.SetSizer( bSizer31 )
		self.AC_panel.Layout()
		bSizer31.Fit( self.AC_panel )
		
		
		
		#bSizer31 = wx.BoxSizer( wx.VERTICAL )
		
		#bSizer28 = wx.BoxSizer( wx.HORIZONTAL )
		
		#sbSizer31 = wx.StaticBoxSizer( wx.StaticBox( self.AC_panel, wx.ID_ANY, u"Source file" ), wx.HORIZONTAL )
		
		#self.AC_source_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.AC_source_label.Wrap( -1 )
		#sbSizer31.Add( self.AC_source_label, 1, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		#self.AC_load_button = wx.Button( self.AC_panel, wx.ID_ANY, u"Load ...", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.AC_load_button.SetDefault() 
		#sbSizer31.Add( self.AC_load_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		#bSizer28.Add( sbSizer31, 6, wx.ALL|wx.EXPAND, 5 )
		
		#bSizer31.Add( bSizer28, 0, wx.EXPAND, 5 )
		
		#bSizer52 = wx.BoxSizer( wx.HORIZONTAL )
		
		#sbSizer51 = wx.StaticBoxSizer( wx.StaticBox( self.AC_panel, wx.ID_ANY, u"Statistics" ), wx.HORIZONTAL )
		
		#gSizer21 = wx.GridSizer( 3, 2, 0, 0 )
		
		#self.m_staticText21 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"Terms", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.m_staticText21.Wrap( -1 )
		#gSizer21.Add( self.m_staticText21, 0, wx.ALL, 5 )
		
		#self.terms_number_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.terms_number_label.Wrap( -1 )
		#gSizer21.Add( self.terms_number_label, 0, wx.ALL, 5 )

		#self.m_staticText51 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"Objects", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.m_staticText51.Wrap( -1 )
		#gSizer21.Add( self.m_staticText51, 0, wx.ALL, 5 )
		
		#self.objects_number_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.objects_number_label.Wrap( -1 )
		#gSizer21.Add( self.objects_number_label, 0, wx.ALL, 5 )
		
		#self.m_staticText18 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"File type", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.m_staticText18.Wrap( -1 )
		#gSizer21.Add( self.m_staticText18, 0, wx.ALL, 5 )
		
		#self.file_type_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		#self.file_type_label.Wrap( -1 )
		#gSizer21.Add( self.file_type_label, 0, wx.ALL, 5 )
		
		##self.m_staticText61 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"Categories", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.m_staticText61.Wrap( -1 )
		##gSizer21.Add( self.m_staticText61, 0, wx.ALL, 5 )
		
		##self.m_staticText71 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"3", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.m_staticText71.Wrap( -1 )
		##gSizer21.Add( self.m_staticText71, 0, wx.ALL, 5 )
		
		#sbSizer51.Add( gSizer21, 1, wx.ALIGN_CENTER|wx.SHAPED, 2 )
		
		#bSizer52.Add( sbSizer51, 0, wx.ALL, 5 )
		
		##sbSizer61 = wx.StaticBoxSizer( wx.StaticBox( self.AC_panel, wx.ID_ANY, u"Additional information" ), wx.VERTICAL )
		
		##bSizer511 = wx.BoxSizer( wx.VERTICAL )
		
		##gbSizer1 = wx.GridBagSizer( 0, 0 )
		##gbSizer1.SetFlexibleDirection( wx.BOTH )
		##gbSizer1.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		##self.m_staticText18 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"File type", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.m_staticText18.Wrap( -1 )
		##gbSizer1.Add( self.m_staticText18, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		##self.file_type_label = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.file_type_label.Wrap( -1 )
		##gbSizer1.Add( self.file_type_label, wx.GBPosition( 0, 1 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		##self.m_staticText20 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"MyLabel", wx.Point( -1,-1 ), wx.DefaultSize, 0 )
		##self.m_staticText20.Wrap( -1 )
		##gbSizer1.Add( self.m_staticText20, wx.GBPosition( 1, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		##self.m_staticText211 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.m_staticText211.Wrap( -1 )
		##gbSizer1.Add( self.m_staticText211, wx.GBPosition( 1, 1 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		##self.m_staticText22 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"MyLabel", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.m_staticText22.Wrap( -1 )
		##gbSizer1.Add( self.m_staticText22, wx.GBPosition( 2, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		##self.m_staticText23 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		##self.m_staticText23.Wrap( -1 )
		##gbSizer1.Add( self.m_staticText23, wx.GBPosition( 2, 1 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		##bSizer511.Add( gbSizer1, 1, wx.EXPAND, 5 )
		
		##sbSizer61.Add( bSizer511, 1, wx.EXPAND, 5 )
		
		##bSizer52.Add( sbSizer61, 0, wx.ALL|wx.EXPAND, 5 )
		
		#bSizer31.Add( bSizer52, 0, wx.EXPAND, 5 )
		
		#self.AC_panel.SetSizer( bSizer31 )
		#self.AC_panel.Layout()
		#bSizer31.Fit( self.AC_panel )
		
		self.ac_load_gui = AC_load_gui(self, self.real_parent)
		self.Bind(wx.EVT_BUTTON, self.OnLoadButton, id=self.AC_load_button.GetId())
		self.AC_load_outcome_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_LOAD_AC, self.OnLoadDone)
		self.AC_get_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_GET, self.OnGet)
		self.objects_number = None
		self.terms_number = None
		self._reset()
#






	def _reset(self):
		self._update()
#







	def _freeze(self):
		self.AC_load_button.Disable()
#






	def _unfreeze(self):
		self.AC_load_button.Enable()
#




	def update(self):
		if DEBUG_LEVEL>2:
			print "ACPanel: update()"
		self.real_parent.update()
#






	def OnLoadButton(self, event):
		self.ac_load_gui.ShowModal()
#






	def _update(self):
		if DEBUG_LEVEL>1:
			print "ACPanel: _update()"
			
		if self.real_parent.AC_status:
			self.AC_source_label.SetLabel(os.path.basename(self.real_parent.params_AC['filename']))
			self.objects_number_label.SetLabel(str(self.objects_number))
			self.terms_number_label.SetLabel(str(self.terms_number))
			self.file_type_label.SetLabel(str(self.real_parent.params_AC['type']))
			if str(self.real_parent.params_AC['type']) == 'gaf-2.0':
				ook = False
				if 'params' in self.real_parent.params_AC:
					if 'filter' in self.real_parent.params_AC['params']:
						if 'EC' in self.real_parent.params_AC['params']['filter']:
							if 'EC' in self.real_parent.params_AC['params']['filter']['EC']:
								if 'IEA' in self.real_parent.params_AC['params']['filter']['EC']['EC']:
									if 'inclusive' in self.real_parent.params_AC['params']['filter']['EC']:
										if not self.real_parent.params_AC['params']['filter']['EC']['inclusive']:
											self.IEA_annotations_label.SetLabel("ignored")
											ook = True
				if not ook:
					self.IEA_annotations_label.SetLabel("considered")
			else:
				self.IEA_annotations_label.SetLabel("-")
			self.AC_status_label.SetLabel("Annotation Corpus correctly loaded")
		else:
			if not self.real_parent.GO_status:
				self.AC_status_label.SetLabel("No Annotation Corpus currently loaded.\nA Gene Ontology should be loaded before")
			else:
				self.AC_status_label.SetLabel("No Annotation Corpus currently loaded")
			self.AC_source_label.SetLabel(u"No file loaded.")
			self.objects_number_label.SetLabel("-")
			self.terms_number_label.SetLabel("-")
			self.file_type_label.SetLabel("-")
			self.IEA_annotations_label.SetLabel("-")
#





	def OnLoadDone(self, event):
		if DEBUG_LEVEL>0:
			print "ACPanel: OnLoadDone()"
		data = event.data

		if data[0] == WorkProcess.CMD_LOAD_AC:
			if data[1] == WorkProcess.ANSWER_PROCESSED:
				if data[2] == WorkProcess. RESULT_OK:
					if DEBUG_LEVEL>1:
						print "AC load outcome: Load successful."
					#self.node_info = data[3]
					#self.edge_info = data[4]
					self.ac_load_gui.Hide()
					
				else:
					if DEBUG_LEVEL>2:
						print "AC load outcome: Load Fail."
					#self.node_info = None
					#self.edge_info = None
			elif data[1] == WorkProcess.ANSWER_PROCESSING:
				if DEBUG_LEVEL>2:
					print "AC load outcome: Load in progress."
				return
			elif data[1] == WorkProcess.ANSWER_IGNORED:
				if DEBUG_LEVEL>2:
					print "AC load outcome: Load request ignored."
			else:
				if DEBUG_LEVEL>2:
					print "AC load outcome: Unknown answer."
		#self.unfreeze()
		self.ac_load_gui.unfreeze()
		#if self.AC_get_handle == None:
			#self.AC_get_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_GET, self.OnGetParams)
		self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_GET, WorkProcess.CMD_GET_PARAMS, WorkProcess.CMD_GET_PARAMS_AC))
		self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_GET, WorkProcess.CMD_GET_AC, WorkProcess.CMD_GET_AC_OBJECTS_NUMBER))
		self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_GET, WorkProcess.CMD_GET_AC, WorkProcess.CMD_GET_AC_TERMS_NUMBER))
		#print "##############################################################################################################"
		#print str((WorkProcess.CMD_GET, WorkProcess.CMD_GET_AC, WorkProcess.CMD_GET_AC_TERMS_NUMBER))
		#print "##############################################################################################################"
		self.real_parent.update()
		event.Skip()
#




	def OnGet(self, event):
		if DEBUG_LEVEL>0:
			print "ACPanel: OnGetParams()"
		data = event.data
		if data[0] == WorkProcess.CMD_GET:
			if data[1] == WorkProcess.CMD_GET_PARAMS:
				if data[2] == WorkProcess. CMD_GET_PARAMS_AC:
					data = data[3]
					self.real_parent.params_AC = data
			elif data[1] == WorkProcess.CMD_GET_AC:
				if DEBUG_LEVEL>3:
					print "ACPanel: OnGetParams() GET_AC"
				if data[2] == WorkProcess. CMD_GET_AC_OBJECTS_NUMBER:
					self.objects_number = data[3]
				if data[2] == WorkProcess. CMD_GET_AC_TERMS_NUMBER:
					self.terms_number = data[3]
				self._update()
		event.Skip()
#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# AC load window #

class AC_load_gui ( wx.Dialog ):

	param_filename = None
	param_filetype = None
	params = {}
	
	def __init__( self, parent , real_parent):
		wx.Dialog.__init__ ( self, parent, id = wx.ID_ANY, title = u"Load Annotation Corpus", pos = wx.DefaultPosition, size = wx.Size( 685,442 ), style = wx.DEFAULT_DIALOG_STYLE )
		
		self.parent = parent
		self.real_parent = real_parent
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		bSizer3 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer22 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer3 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Source file" ), wx.HORIZONTAL )
		
		self.AC_load_source_label = wx.StaticText( self, wx.ID_ANY, u"No file selected", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.AC_load_source_label.Wrap( -1 )
		sbSizer3.Add( self.AC_load_source_label, 1, wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT, 10 )
		
		self.AC_load_select_button = wx.Button( self, wx.ID_ANY, u"Select file...", wx.DefaultPosition, wx.DefaultSize, 0 )
		sbSizer3.Add( self.AC_load_select_button, 0, wx.ALIGN_RIGHT|wx.ALL, 5 )
		
		bSizer22.Add( sbSizer3, 1, wx.ALL, 5 )
		
		bSizer3.Add( bSizer22, 0, wx.ALIGN_CENTER|wx.EXPAND, 5 )
		
		bSizer5 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer22 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Source type" ), wx.VERTICAL )
		
		self.AC_load_type_box = wx.Choicebook( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.CHB_DEFAULT )

# Gaf-2.0 file type

		self.m_panel9 = wx.Panel( self.AC_load_type_box, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer23 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer29 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel9, wx.ID_ANY, u"Options" ), wx.VERTICAL )
		
		gbSizer2 = wx.GridBagSizer( 0, 3 )
		gbSizer2.SetFlexibleDirection( wx.BOTH )
		gbSizer2.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.AC_load_ignore_IEA_check = wx.CheckBox( self.m_panel9, wx.ID_ANY, u"Ignore IEA", wx.DefaultPosition, wx.DefaultSize, 0 )
		gbSizer2.Add( self.AC_load_ignore_IEA_check, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.m_staticText59 = wx.StaticText( self.m_panel9, wx.ID_ANY, u"Check this to ignore IEA annotations (Inferred Electronically Annotations)", wx.DefaultPosition, wx.Size( 450,-1 ), 0 )
		self.m_staticText59.Wrap( -1 )
		self.m_staticText59.SetMaxSize( wx.Size( 450,-1 ) )
		
		gbSizer2.Add( self.m_staticText59, wx.GBPosition( 0, 1 ), wx.GBSpan( 1, 2 ), wx.ALL, 5 )
		
		self.AC_load_simplify_check = wx.CheckBox( self.m_panel9, wx.ID_ANY, u"Simplify", wx.DefaultPosition, wx.DefaultSize, 0 )
		gbSizer2.Add( self.AC_load_simplify_check, wx.GBPosition( 1, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.m_staticText591 = wx.StaticText( self.m_panel9, wx.ID_ANY, u"Check this to simplify the annotation corpus after loading. Removes informartion about annotation reliability, taxonomy, ...", wx.DefaultPosition, wx.Size( 450,-1 ), 0 )
		self.m_staticText591.Wrap( -1 )
		self.m_staticText591.SetMaxSize( wx.Size( 450,-1 ) )
		
		gbSizer2.Add( self.m_staticText591, wx.GBPosition( 1, 1 ), wx.GBSpan( 1, 2 ), wx.ALL, 5 )
		
		self.AC_load_tax_check = wx.CheckBox( self.m_panel9, wx.ID_ANY, u"Select taxonomy", wx.DefaultPosition, wx.DefaultSize, 0 )
		gbSizer2.Add( self.AC_load_tax_check, wx.GBPosition( 2, 0 ), wx.GBSpan( 1, 1 ), wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		bSizer51 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer23 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel9, wx.ID_ANY, u"Taxonomy Id" ), wx.VERTICAL )
		
		self.AC_load_tax_text = wx.TextCtrl( self.m_panel9, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		sbSizer23.Add( self.AC_load_tax_text, 0, wx.ALL, 5 )
		
		bSizer51.Add( sbSizer23, 0, wx.ALL|wx.EXPAND, 5 )
		
		gbSizer2.Add( bSizer51, wx.GBPosition( 2, 1 ), wx.GBSpan( 1, 1 ), wx.ALIGN_CENTER_VERTICAL|wx.BOTTOM|wx.EXPAND, 5 )
		
		self.m_staticText51 = wx.StaticText( self.m_panel9, wx.ID_ANY, u"Select a taxonomy to consider. All the objects not belonging to the selected taxonomy will be ignored.", wx.DefaultPosition, wx.Size( 300,-1 ), 0 )
		self.m_staticText51.Wrap( -1 )
		gbSizer2.Add( self.m_staticText51, wx.GBPosition( 2, 2 ), wx.GBSpan( 1, 1 ), wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		sbSizer29.Add( gbSizer2, 1, wx.BOTTOM|wx.EXPAND, 10 )
		
		bSizer23.Add( sbSizer29, 1, wx.EXPAND, 5 )
		
		self.m_panel9.SetSizer( bSizer23 )
		self.m_panel9.Layout()
		bSizer23.Fit( self.m_panel9 )
		self.AC_load_type_box.AddPage( self.m_panel9, GAF2_FILE_TYPE, True )

# Plain file type

		self.m_panel10 = wx.Panel( self.AC_load_type_box, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer231 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer30 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel10, wx.ID_ANY, u"Options" ), wx.VERTICAL )
		
		fgSizer8 = wx.FlexGridSizer( 3, 3, 15, 5 )
		fgSizer8.SetFlexibleDirection( wx.BOTH )
		fgSizer8.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.m_staticText30 = wx.StaticText( self.m_panel10, wx.ID_ANY, u"Separator", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText30.Wrap( -1 )
		fgSizer8.Add( self.m_staticText30, 0, wx.ALL, 5 )
		
		bSizer46 = wx.BoxSizer( wx.VERTICAL )
		
		self.AC_load_sep_tab_radio = wx.RadioButton( self.m_panel10, wx.ID_ANY, u"[tab]  '\\t'", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer46.Add( self.AC_load_sep_tab_radio, 0, wx.LEFT|wx.RIGHT|wx.TOP|wx.RB_GROUP, 5 )
		
		self.AC_load_sep_space_radio = wx.RadioButton( self.m_panel10, wx.ID_ANY, u"[space]", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer46.Add( self.AC_load_sep_space_radio, 0, wx.LEFT|wx.RIGHT, 5 )
		
		bSizer47 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.AC_load_sep_custom_radio = wx.RadioButton( self.m_panel10, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer47.Add( self.AC_load_sep_custom_radio, 0, wx.LEFT|wx.RIGHT, 5 )
		
		self.AC_sep_custom_text = wx.TextCtrl( self.m_panel10, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 35,-1 ), 0 )
		bSizer47.Add( self.AC_sep_custom_text, 0, wx.RIGHT, 5 )
		
		bSizer46.Add( bSizer47, 1, wx.EXPAND, 5 )
		
		fgSizer8.Add( bSizer46, 1, wx.EXPAND, 5 )
		
		self.m_staticText48 = wx.StaticText( self.m_panel10, wx.ID_ANY, u"Character used to separate fields within each row.", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText48.Wrap( -1 )
		fgSizer8.Add( self.m_staticText48, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		self.m_staticText49 = wx.StaticText( self.m_panel10, wx.ID_ANY, u"Row format", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText49.Wrap( -1 )
		fgSizer8.Add( self.m_staticText49, 0, wx.ALL, 5 )
		
		bSizer48 = wx.BoxSizer( wx.VERTICAL )
		
		self.AC_load_multi_check = wx.CheckBox( self.m_panel10, wx.ID_ANY, u"Multiple associations", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer48.Add( self.AC_load_multi_check, 0, wx.LEFT|wx.RIGHT, 5 )
		
		fgSizer8.Add( bSizer48, 1, wx.EXPAND, 5 )
		
		self.m_staticText50 = wx.StaticText( self.m_panel10, wx.ID_ANY, u"Check this option if each row contains more than one association between objects and GO Terms.", wx.DefaultPosition, wx.Size( 325,-1 ), 0 )
		self.m_staticText50.Wrap( -1 )
		fgSizer8.Add( self.m_staticText50, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		self.m_staticText511 = wx.StaticText( self.m_panel10, wx.ID_ANY, u"Data order", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText511.Wrap( -1 )
		fgSizer8.Add( self.m_staticText511, 0, wx.ALL, 5 )
		
		bSizer49 = wx.BoxSizer( wx.VERTICAL )
		
		self.AC_load_obj_first_radio = wx.RadioButton( self.m_panel10, wx.ID_ANY, u"Object -> GO Term(s)", wx.DefaultPosition, wx.DefaultSize, wx.RB_GROUP )
		bSizer49.Add( self.AC_load_obj_first_radio, 0, wx.LEFT|wx.TOP, 5 )
		
		self.AC_load_term_first_radio = wx.RadioButton( self.m_panel10, wx.ID_ANY, u"GO Term -> Object(s)", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer49.Add( self.AC_load_term_first_radio, 0, wx.BOTTOM|wx.LEFT, 5 )
		
		fgSizer8.Add( bSizer49, 1, wx.EXPAND, 5 )
		
		self.m_staticText52 = wx.StaticText( self.m_panel10, wx.ID_ANY, u"Select whether the first field of each row is an object or a GO Term.", wx.DefaultPosition, wx.Size( 325,-1 ), 0 )
		self.m_staticText52.Wrap( -1 )
		fgSizer8.Add( self.m_staticText52, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		sbSizer30.Add( fgSizer8, 0, wx.BOTTOM|wx.EXPAND, 5 )
		
		bSizer231.Add( sbSizer30, 1, wx.EXPAND, 5 )
		
		self.m_panel10.SetSizer( bSizer231 )
		self.m_panel10.Layout()
		bSizer231.Fit( self.m_panel10 )
		self.AC_load_type_box.AddPage( self.m_panel10, PLAIN_FILE_TYPE, False )
		sbSizer22.Add( self.AC_load_type_box, 1, wx.EXPAND |wx.ALL, 5 )
		
		bSizer5.Add( sbSizer22, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer3.Add( bSizer5, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer4 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.AC_load_load_button = wx.Button( self, wx.ID_ANY, u"Load", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.AC_load_load_button.SetDefault() 
		bSizer4.Add( self.AC_load_load_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		self.AC_load_cancel_button = wx.Button( self, wx.ID_ANY, u"Cancel", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer4.Add( self.AC_load_cancel_button, 0, wx.ALL, 5 )
		
		bSizer3.Add( bSizer4, 0, wx.ALIGN_CENTER, 5 )
		
		self.SetSizer( bSizer3 )
		self.Layout()
		
		self.Centre( wx.BOTH )

		self.Bind(wx.EVT_BUTTON, self.OnFileBrowse, id=self.AC_load_select_button.GetId())
		self.Bind(wx.EVT_BUTTON, self.OnCancel, id=self.AC_load_cancel_button.GetId())
		self.Bind(wx.EVT_BUTTON, self.OnLoad, id=self.AC_load_load_button.GetId())

		#self.Bind(wx.EVT_CHOICEBOOK_PAGE_CHANGED, self.OnSelectType, id=self.AC_load_type_box.GetId())
		#self.Bind(wx.EVT_CHECKBOX, self.collect_gaf2_file_parms, id=self.AC_load_simplify_check.GetId())
		#self.Bind(wx.EVT_CHECKBOX, self.collect_gaf2_file_parms, id=self.AC_load_ignore_IEA_check.GetId())
		self.Bind(wx.EVT_CHECKBOX, self.OnTax, id=self.AC_load_tax_check.GetId())
		#self.real_parent.Bind(wx.EVT_TEXT, self.collect_gaf2_file_parms, id=self.AC_load_tax_text.GetId())

		self._reset_()
		#self.Bind(wx.EVT_RADIOBUTTON, self.OnSep, id=self.output_sep_tab_radio.GetId())
		#self.Bind(wx.EVT_RADIOBUTTON, self.OnSep, id=self.output_sep_space_radio.GetId())
#


	def _reset_(self):
		param_filename = None
		param_filetype = None
		self.params = {}
		self._set_file_name(None)
		#self.params[PlainAnnotationCorpus.SEPARATOR] = '\t'
		#self.params[PlainAnnotationCorpus.SEPARATOR] = ' '
		#self.params[PlainAnnotationCorpus.SEPARATOR] = self.text_separator.GetValue()
		#self.params[PlainAnnotationCorpus.MANYASSPERROW] = self.check_manyperline.GetValue()
		#self.params[PlainAnnotationCorpus.TERMFIRST] = self.radius_rowformat_2.GetValue()
		#self.params[GAF2AnnotationCorpus.SIMPLIFY] = self.check_simplify.GetValue()
		#self.params[AnnotationCorpus.FILTER_PARAM] = {}
		#self.params[AnnotationCorpus.FILTER_PARAM]['EC'] = {'EC':'IEA'}
		#self.params[AnnotationCorpus.FILTER_PARAM]['taxonomy'] = int(self.text_filtertax.GetValue())
#


#

	def OnSelectType(self, event):
		filetype = self.AC_load_type_box.GetPageText(self.AC_load_type_box.GetSelection())
		if filetype in AC_FILE_TYPES:
			self.real_parent.params_AC['type'] = AC_FILE_TYPES[filetype]
#

	def OnCancel(self, event):
		self.Hide()
#

	def OnFileBrowse(self, event):
		dialog = wx.FileDialog(None, style = wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			self._set_file_name(dialog.GetPath())
#

	def OnTax(self, event):
		if self.AC_load_tax_check.GetValue():
			self.AC_load_tax_text.Enable()
		else:
			self.AC_load_tax_text.Disable()
#





	def collect_gaf2_file_parms(self, event):
		if 'params' in self.real_parent.params_AC:
			del self.real_parent.params_AC['params']
		self.real_parent.params_AC['params'] = {}

		if not 'filter' in self.real_parent.params_AC['params']:
			self.real_parent.params_AC['params']['filter'] = {}

		if self.AC_load_ignore_IEA_check.GetValue():
		 self.real_parent.params_AC['params']['filter']['EC'] = {}
		 self.real_parent.params_AC['params']['filter']['EC']['EC'] = {'IEA':None}
		 self.real_parent.params_AC['params']['filter']['EC']['inclusive'] = False
		elif 'EC' in self.real_parent.params_AC['params']['filter']:
			del self.real_parent.params_AC['params']['filter']['EC']
			
		if self.AC_load_simplify_check.GetValue():
		 self.real_parent.params_AC['params']['simplify'] = True
		else:
			self.real_parent.params_AC['params']['simplify'] = False

		if self.AC_load_tax_check.GetValue():
			self.AC_load_tax_text.Enable()
			self.real_parent.params_AC['params']['filter']['taxonomy'] = self.AC_load_tax_text.GetValue()
		else:
			self.AC_load_tax_text.Disable()
			if 'taxonomy' in self.real_parent.params_AC['params']['filter']:
				del self.real_parent.params_AC['params']['filter']['taxonomy']
#






	def collect_plain_file_parms(self, event):
		if 'params' in self.real_parent.params_AC:
			del self.real_parent.params_AC['params']
		self.real_parent.params_AC['params'] = {}

		if self.AC_load_sep_tab_radio.GetValue():
			self.real_parent.params_AC['params']['separator'] = "\t"
		if self.AC_load_sep_space_radio.GetValue():
			self.real_parent.params_AC['params']['separator'] = " "
		if self.AC_load_sep_custom_radio.GetValue():
			self.real_parent.params_AC['params']['separator'] = self.AC_sep_custom_text.GetValue()
			
		if self.AC_load_obj_first_radio.GetValue():
			self.real_parent.params_AC['params']['term first'] = False
		if self.AC_load_term_first_radio.GetValue():
			self.real_parent.params_AC['params']['term first'] = True
		if self.AC_load_multi_check.GetValue():
			self.real_parent.params_AC['params']['multiple'] = True
		else:
			self.real_parent.params_AC['params']['multiple'] = False
#




	def _set_file_name(self, fn):
		self.param_filename = fn
		if fn == None:
			self.AC_load_source_label.SetLabel(u"No file selected.")
			self.AC_load_load_button.Disable()
		else:
			self.AC_load_source_label.SetLabel(self.param_filename)
			self.AC_load_load_button.Enable()
#



	def freeze(self):
		self.AC_load_cancel_button.Disable()
		self.AC_load_load_button.Disable()
		self.AC_load_select_button.Disable()
		self.AC_load_type_box.Disable()
#

	def unfreeze(self):
		self.AC_load_cancel_button.Enable()
		self.AC_load_load_button.Enable()
		self.AC_load_select_button.Enable()
		self.AC_load_type_box.Enable()
#

	def OnLoad(self, data):
		self.freeze()
		#self.filetype = self.AC_load_type_box.GetPageText(self.AC_load_type_box.GetSelection())
		#if self.filetype in AC_FILE_TYPES:
			#self.param_filetype = AC_FILE_TYPES[self.filetype]
			#print AC_FILE_TYPES[self.filetype]
		#if self.param_filetype == GAF2_FILE_TYPE:
		#self.AC_load_ignore_haspart_check.Disable()
		#self.AC_load_ignore_regulates_check.Disable()
		self.real_parent.params_AC['filename'] = self.param_filename
		self.OnSelectType(None)
		if self.real_parent.params_AC['type'] == GAF2_FILE_TYPE_REF:
			self.collect_gaf2_file_parms(None)
		elif self.real_parent.params_AC['type'] == PLAIN_FILE_TYPE_REF:
			self.collect_plain_file_parms(None)

		#print self.real_parent.params_AC
		self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_LOAD_AC, self.real_parent.params_AC))
#

