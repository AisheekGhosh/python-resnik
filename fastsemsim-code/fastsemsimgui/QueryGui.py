import wx
import WorkProcess
import os 
import sys

DEBUG_LEVEL = 0

################################################################
	
class QueryPanel(wx.Panel):
	def __init__( self, real_parent, parent, id, pos, size, style):
		wx.Panel.__init__ ( self, parent, id, pos, size, style)
		self.query_panel = self # temporary workaround
		self.real_parent = real_parent

		bSizer50 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer512 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer312 = wx.StaticBoxSizer( wx.StaticBox( self.query_panel, wx.ID_ANY, u"Query" ), wx.VERTICAL )
		
		self.query_text = wx.TextCtrl( self.query_panel, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 200,200 ), wx.TE_MULTILINE|wx.TE_PROCESS_ENTER|wx.TE_PROCESS_TAB)

		self.query_text.SetMinSize( wx.Size( 200,270 ) )
		
		sbSizer312.Add( self.query_text, 0, wx.ALL, 5 )
		
		bSizer512.Add( sbSizer312, 0, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.BOTTOM|wx.LEFT, 5 )
		
		bSizer54 = wx.BoxSizer( wx.VERTICAL )
		
		
		sbSizer32 = wx.StaticBoxSizer( wx.StaticBox( self.query_panel, wx.ID_ANY, u"Query format" ), wx.VERTICAL )
		
		self.query_type_choices = query_types = [ "Pairs", "List" ]
		self.query_type_ref = {"Pairs":WorkProcess.QUERY_PAIRS, "List":WorkProcess.QUERY_LIST}
		
		self.query_type_box = wx.ComboBox( self.query_panel, wx.ID_ANY, "", wx.DefaultPosition, wx.Size( 120,-1 ), self.query_type_choices, wx.CB_READONLY )
		self.query_type_box.SetSelection(0)
		sbSizer32.Add( self.query_type_box, 0, wx.ALL, 5 )
		
		bSizer54.Add( sbSizer32, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		sbSizer33 = wx.StaticBoxSizer( wx.StaticBox( self.query_panel, wx.ID_ANY, u"Help" ), wx.VERTICAL )
		
		self.m_staticText53 = wx.StaticText( self.query_panel, wx.ID_ANY, u"Manually enter a query in the field on the left. Separate the entries with tabs. Spaces will not be considered as separators. You can also load the query from a file or use the whole annotation corpus.", wx.DefaultPosition, wx.Size( 225,150 ), 0 )
		self.m_staticText53.Wrap( -1 )
		sbSizer33.Add( self.m_staticText53, 0, wx.ALL, 5 )
		
		bSizer541 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer33.Add( bSizer541, 1, wx.EXPAND, 5 )
		
		bSizer54.Add( sbSizer33, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.BOTTOM, 5 )
		
		bSizer53 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer55 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.query_from_file_button = wx.Button( self.query_panel, wx.ID_ANY, u"Load from File...", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		bSizer55.Add( self.query_from_file_button, 0, wx.ALL, 5 )
		
		self.query_reset_button = wx.Button( self.query_panel, wx.ID_ANY, u"Clear", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer55.Add( self.query_reset_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer53.Add( bSizer55, 1, wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.query_from_AC_button = wx.Button( self.query_panel, wx.ID_ANY, u"Use the entire Annotation Corpus", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		bSizer53.Add( self.query_from_AC_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer54.Add( bSizer53, 0, wx.EXPAND, 5 )
		
		
		bSizer512.Add( bSizer54, 1, wx.BOTTOM|wx.TOP, 5 )
		
		bSizer50.Add( bSizer512, 0, wx.EXPAND | wx.BOTTOM, 5 )
		
		self.query_status_label = wx.StaticText( self.query_panel, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.query_status_label.Wrap( -1 )
		self.query_status_label.SetFont( wx.Font( wx.NORMAL_FONT.GetPointSize(), 70, 90, 92, False, wx.EmptyString ) )
		
		bSizer50.Add( self.query_status_label, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.m_staticline1 = wx.StaticLine( self.query_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LI_HORIZONTAL )
		bSizer50.Add( self.m_staticline1, 0, wx.EXPAND |wx.ALL, 5 )
		
		sbSizer30 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.controls_start_button = wx.Button( self.query_panel, wx.ID_ANY, u"Start", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.controls_start_button.SetDefault() 
		sbSizer30.Add( self.controls_start_button, 0, wx.ALL, 5 )
		
		sbSizer311 = wx.StaticBoxSizer( wx.StaticBox( self.query_panel, wx.ID_ANY, u"Log" ), wx.VERTICAL )
		
		self.controls_log = wx.TextCtrl( self.query_panel, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH2|wx.TE_WORDWRAP )
		self.controls_log.SetMinSize( wx.Size( -1,-1 ) )
		
		sbSizer311.Add( self.controls_log, 0, wx.ALL|wx.EXPAND, 3 )
		
		sbSizer30.Add( sbSizer311, 1, wx.LEFT, 5 )
		
		bSizer50.Add( sbSizer30, 0, wx.BOTTOM|wx.EXPAND|wx.LEFT|wx.RIGHT, 5 )
		
		self.query_panel.SetSizer( bSizer50 )
		self.query_panel.Layout()
		bSizer50.Fit( self.query_panel )

		
		self.query_load_gui = Query_load_gui(self, self.real_parent)
		
		self.real_parent.Bind(wx.EVT_BUTTON, self.OnStartCmd, id=self.controls_start_button.GetId())
		self.real_parent.Bind(wx.EVT_BUTTON, self.OnLoadFromFile, id=self.query_from_file_button.GetId())
		self.real_parent.Bind(wx.EVT_BUTTON, self.OnLoadFromAC, id=self.query_from_AC_button.GetId())
		self.real_parent.Bind(wx.EVT_BUTTON, self.OnReset, id=self.query_reset_button.GetId())
		self.real_parent.Bind(wx.EVT_COMBOBOX, self.OnSelectQueryType, id=self.query_type_box.GetId())
		self.real_parent.Bind(wx.EVT_TEXT, self.OnQueryUpdate, id=self.query_text.GetId())

		self._reset()
#



	def _reset(self):
		if DEBUG_LEVEL > 0:
			print "QueryGui: _reset()"
		self.query_text.Enable()
		self.query_text.SetValue('')
		self.query_type_box.Enable()
		self.query_type_box.SetSelection(0)
		
		
		self.real_parent.query = None
		self.real_parent.params_query['type'] = None
		self.OnSelectQueryType(None)
		self.real_parent.params_query['source'] = WorkProcess.QUERY_FROM_GUI

		self.query_to_update = True
		self._update()
#




	def update_start_button(self):
		if self.real_parent.status == WorkProcess.STATUS_WAIT:
			self.controls_start_button.SetLabel('Start')
			if self.real_parent.GO_status and self.real_parent.AC_status:
				self.controls_start_button.Enable()
			else:
				self.controls_start_button.Disable()
		elif self.real_parent.status == WorkProcess.STATUS_RUN:
			self.controls_start_button.SetLabel('Stop')
			self.controls_start_button.Enable()
		else:
			self.controls_start_button.Disable()

#




	def _update(self):
		self.update_start_button()

		if not self.real_parent.GO_status or not self.real_parent.AC_status:
			self.query_from_AC_button.Disable()
		elif self.real_parent.GO_status and self.real_parent.AC_status:
			self.query_from_AC_button.Enable()

		if self.real_parent.params_query['source'] == WorkProcess.QUERY_FROM_GUI:
			pass
		elif self.real_parent.params_query['source'] == WorkProcess.QUERY_FROM_AC:
			pass
		elif self.real_parent.params_query['source'] == WorkProcess.QUERY_FROM_FILE:
			pass
		else:
			self.real_parent.params_query['source'] = WorkProcess.QUERY_FROM_GUI
			return self._update()

		if 'type' in self.real_parent.params_query and not self.real_parent.params_query['type']==None:
			self.query_type_box.SetSelection(self.real_parent.params_query['type'])
			
		if not self.real_parent.GO_status:
			self.query_status_label.SetLabel("Please load a valid Gene Ontology and Annotation Corpus.")
		elif  not self.real_parent.AC_status:
			self.query_status_label.SetLabel("Please load a valid Annotation Corpus.")
		else:
			self.query_status_label.SetLabel("Specify a query.")
#








	def _freeze(self):
		self.query_from_AC_button.Disable()
		self.query_from_file_button.Disable()
		self.query_reset_button.Disable()
		self.query_text.Disable()
		self.query_type_box.Disable()
#






	def _unfreeze(self):
		self.query_from_AC_button.Enable()
		self.query_from_file_button.Enable()
		self.query_reset_button.Enable()
		self.query_text.Enable()
		self.query_type_box.Enable()
#






	def OnStartCmd(self, event):
		if self.real_parent.status == WorkProcess.STATUS_WAIT:
			self.real_parent.start()
		elif self.real_parent.status == WorkProcess.STATUS_RUN:
			self.real_parent.stop()
			
#







	def setBitmap(self, handle, value):
		if value:
			handle.SetBitmap(self.ok_bitmap)
		else:
			handle.SetBitmap(self.warning_bitmap)
#

	def SetGOStatus(self, status):
		self.go_ok = status
		#self.setBitmap(self.go_status_bitmap, status)
		
	def SetAcStatus(self, status):
		self.ac_ok = status
		#self.setBitmap(self.ac_status_bitmap, status)
		
	def SetQueryStatus(self, status):
		self.query_ok = status
		#self.setBitmap(self.query_status_bitmap, status)
		##self.activateGoCmd()
		
	def SetOutputCtrlStatus(self, status):
		self.outputctrl_ok = status
		#if self.show_pics:
		#self.setBitmap(self.output_status_bitmap, status)
		##self.activateGoCmd()

	def SetSSCtrlStatus(self, status):
		#if self.show_pics:
		#self.setBitmap(self.ss_status_bitmap, status)
		##self.activateGoCmd()
		pass
#


	def OnReset(self, event):
		self._reset()
#


	def OnLoadFromGui(self, event):
		temp_query = self.query_panel.query_text.GetValue()
		current_query = []
		temp_query = temp_query.split("\n")
		for i in temp_query:
			i = i.rstrip("\n")
			i = i.rstrip("\r")
			if len(i) == 0:
				continue
			#i = i.rstrip(" ")
			if self.real_parent.params_query['type'] == WorkProcess.QUERY_LIST:
				line = i.split("\t")
				for j in line:
					if not str(j) == "" and not str(j) == " ":
						current_query.append(j)
			elif self.real_parent.params_query['type'] == WorkProcess.QUERY_PAIRS:
				line = i.split("\t")
				if len(line) >= 2:
					if not str(line[0]) == "" and not str(line[0]) == " " and not str(line[1]) == "" and not str(line[1]) == " ":
						current_query.append((line[0],line[1]))

		self.real_parent.query = current_query
		self.query_to_update = False
#





	def OnLoadFromAC(self, event):
		if not self.real_parent.AC_status:
			self._update()
		else:
			self.query_get_AC_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_GET, self.OnLoadFromACDone)
			self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_GET, WorkProcess.CMD_GET_AC, WorkProcess.CMD_GET_AC_OBJECTS))


	def OnLoadFromACDone(self, event):
		if DEBUG_LEVEL > 0:
			print "OnLoadFromACDone"
		data = event.data
		if data[1] == WorkProcess.CMD_GET_AC and data[2] == WorkProcess.CMD_GET_AC_OBJECTS:
			data = data[3]
			self.real_parent.communication_thread.unregister_callback(self.query_get_AC_handle)
			
			self.real_parent.query = data
			
			temp = ("\n".join([str(num) for num in data]))
			temp += "\n"
			#self.query_text.Disable()
			self.query_text.SetValue(temp)
			#self.query_type_box.Disable()
			self.query_type_box.SetSelection(1)

			self.real_parent.params_query['type'] = self.query_type_ref[str(self.query_type_box.GetString(self.query_type_box.GetCurrentSelection()))]
			#self.real_parent.params_query['query'] = self.query
			self.real_parent.params_query['source'] = WorkProcess.QUERY_FROM_AC
			self.query_to_update = False
		event.Skip()
#



	def OnLoadFromFile(self, event):
		self.query_load_gui.ShowModal()
#





	def OnLoadFromFileDone(self, event):
		self.query_type_box.SetSelection(self.query_load_gui.query_file_type_box.GetSelection())
		#self.query_type_box.Disable()
		#self.query_text.Disable()
	
		#self.query = self.query_load_gui.current_query
		self.real_parent.query = self.query_load_gui.current_query
		self.query_text.SetValue(self.text_query)
		
		self.real_parent.params_query['type'] = self.query_type_ref[str(self.query_type_box.GetString(self.query_type_box.GetCurrentSelection()))]
		self.real_parent.params_query['source'] = WorkProcess.QUERY_FROM_FILE
		self.real_parent.params_query['filename'] = self.query_load_gui.filename
		self.query_to_update = False
#




	def OnSelectQueryType(self, event):
		self.real_parent.params_query['type'] = self.query_type_ref[str(self.query_type_box.GetString(self.query_type_box.GetCurrentSelection()))]
		self.query_to_update = True
#


	def OnQueryUpdate(self, event):
		self.query_to_update = True
		self.real_parent.params_query['source'] = WorkProcess.QUERY_FROM_GUI
		#print self.query_text.GetValue()

#








#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Query LOAD GUI

class Query_load_gui ( wx.Dialog ):
	
	def __init__( self, parent , real_parent):
		wx.Dialog.__init__ ( self, parent, id = wx.ID_ANY, title = u"Load Query from File", pos = wx.DefaultPosition, size = wx.Size( 480,415 ), style = wx.DEFAULT_DIALOG_STYLE )
		
		self.parent = parent
		self.real_parent = real_parent
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		bSizer521 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer22 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer34 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Source file" ), wx.HORIZONTAL )
		
		self.query_input_file_label = wx.StaticText( self, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.query_input_file_label.Wrap( -1 )
		sbSizer34.Add( self.query_input_file_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		bSizer22.Add( sbSizer34, 1, wx.LEFT|wx.RIGHT, 5 )
		
		bSizer20 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.query_pick_file_button = wx.Button( self, wx.ID_ANY, u"Select file...", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer20.Add( self.query_pick_file_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer22.Add( bSizer20, 0, wx.ALIGN_RIGHT|wx.RIGHT|wx.TOP, 5 )
		
		bSizer521.Add( bSizer22, 0, wx.ALIGN_CENTER|wx.EXPAND|wx.TOP, 5 )
		
		bSizer55 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer56 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer29 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Options" ), wx.VERTICAL )
		
		fgSizer8 = wx.FlexGridSizer( 3, 3, 15, 5 )
		fgSizer8.SetFlexibleDirection( wx.BOTH )
		fgSizer8.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.m_staticText30 = wx.StaticText( self, wx.ID_ANY, u"Separator", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText30.Wrap( -1 )
		fgSizer8.Add( self.m_staticText30, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		bSizer46 = wx.BoxSizer( wx.VERTICAL )
		
		self.query_sep_none_radio = wx.RadioButton( self, wx.ID_ANY, u"None", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer46.Add( self.query_sep_none_radio, 0, wx.LEFT|wx.RIGHT|wx.TOP, 5 )
		
		self.query_sep_tab_radio = wx.RadioButton( self, wx.ID_ANY, u"[tab]  '\\t'", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer46.Add( self.query_sep_tab_radio, 0, wx.LEFT|wx.RIGHT, 5 )
		
		self.query_sep_space_radio = wx.RadioButton( self, wx.ID_ANY, u"[space]", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer46.Add( self.query_sep_space_radio, 0, wx.LEFT|wx.RIGHT, 5 )
		
		bSizer47 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.query_sep_custom_radio = wx.RadioButton( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer47.Add( self.query_sep_custom_radio, 0, wx.LEFT|wx.RIGHT|wx.TOP, 5 )
		
		self.query_sep_custom_text = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 35,-1 ), 0 )
		bSizer47.Add( self.query_sep_custom_text, 0, wx.RIGHT|wx.TOP, 5 )
		
		bSizer46.Add( bSizer47, 1, wx.EXPAND, 5 )
		
		fgSizer8.Add( bSizer46, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 5 )
		
		self.m_staticText48 = wx.StaticText( self, wx.ID_ANY, u"Character used to separate fields within each row.", wx.DefaultPosition, wx.Size( 200,40 ), 0 )
		self.m_staticText48.Wrap( -1 )
		fgSizer8.Add( self.m_staticText48, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		self.m_staticText511 = wx.StaticText( self, wx.ID_ANY, u"Query type", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText511.Wrap( -1 )
		fgSizer8.Add( self.m_staticText511, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

		self.query_file_type_box = wx.ComboBox( self, wx.ID_ANY, '', wx.DefaultPosition, wx.Size( 100,-1 ), self.parent.query_type_choices, wx.CB_READONLY )
		self.query_file_type_box.SetSelection(0)
		fgSizer8.Add( self.query_file_type_box, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		self.m_staticText52 = wx.StaticText( self, wx.ID_ANY, u"Select whether the input is a list of objects or a set of pairs.", wx.DefaultPosition, wx.Size( 200,60 ), 0 )
		self.m_staticText52.Wrap( -1 )
		fgSizer8.Add( self.m_staticText52, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		self.m_staticText46 = wx.StaticText( self, wx.ID_ANY, u"Process mode", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText46.Wrap( -1 )
		fgSizer8.Add( self.m_staticText46, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		self.query_background_check = wx.CheckBox( self, wx.ID_ANY, u"In background", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		fgSizer8.Add( self.query_background_check, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		#self.m_staticText45 = wx.StaticText( self, wx.ID_ANY, u"Select this option to load the query without displaying it. Useful for huge queries.", wx.DefaultPosition, wx.Size( 200,50 ), 0 )
		self.m_staticText45 = wx.StaticText( self, wx.ID_ANY, u"", wx.DefaultPosition, wx.Size( 200,50 ), 0 )
		self.m_staticText45.Wrap( -1 )
		fgSizer8.Add( self.m_staticText45, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		sbSizer29.Add( fgSizer8, 1, wx.ALL|wx.EXPAND, 5 )
		
		bSizer56.Add( sbSizer29, 1, wx.EXPAND, 5 )
		
		bSizer55.Add( bSizer56, 1, wx.EXPAND, 5 )
		
		bSizer521.Add( bSizer55, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer83 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.query_load_button = wx.Button( self, wx.ID_ANY, u"Load", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer83.Add( self.query_load_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		self.query_cancel_button = wx.Button( self, wx.ID_ANY, u"Cancel", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer83.Add( self.query_cancel_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		
		bSizer521.Add( bSizer83, 0, wx.ALIGN_CENTER|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		self.SetSizer( bSizer521 )
		self.Layout()
		self.Centre( wx.BOTH )
		
# Bind events to controls
		self.Bind(wx.EVT_BUTTON, self.OnFileBrowse, id=self.query_pick_file_button.GetId())
		self.Bind(wx.EVT_BUTTON, self.OnCancel, id=self.query_cancel_button.GetId())
		self.Bind(wx.EVT_BUTTON, self.OnLoad, id=self.query_load_button.GetId())
		self.Bind(wx.EVT_COMBOBOX, self.OnSelectQueryType, id=self.query_file_type_box.GetId())
		self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectSep, id=self.query_sep_custom_radio.GetId())
		self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectSep, id=self.query_sep_tab_radio.GetId())
		self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectSep, id=self.query_sep_space_radio.GetId())
		self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectSep, id=self.query_sep_none_radio.GetId())
		self.Bind(wx.EVT_TEXT, self.OnSelectSep, id=self.query_sep_custom_text.GetId())
		self.Bind(wx.EVT_CHECKBOX, self.OnBackgroundCheck, id=self.query_background_check.GetId())
		#self.Bind(wx.EVT_BUTTON, self.OnFileBrowse, id=self.AC_load_select_button.GetId())
		#self.Bind(wx.EVT_BUTTON, self.OnCancel, id=self.AC_load_cancel_button.GetId())
		#self.Bind(wx.EVT_BUTTON, self.OnLoad, id=self.AC_load_load_button.GetId())
		#self.Bind(wx.EVT_CHECKBOX, self.ignore_what, id=self.GO_load_ignore_haspart_check.GetId())
		#self.Bind(wx.EVT_CHECKBOX, self.ignore_what, id=self.GO_load_ignore_regulates_check.GetId())
		#self.box_filetype = wx.ComboBox(self.panel, wx.ID_ANY, choices=AC_FILE_TYPES, style=wx.CB_READONLY, size=(-1,-1)) #size=(100,30)
		#self.Bind(wx.EVT_COMBOBOX, self.OnSelectType, id=self.box_filetype.GetId())
		self.filename = None
		self._set_file_name(None)
		self._reset()
#




	def _reset(self):
		self.query_type = None
		self.file_sep = None
		self.OnSelectSep(None)
		self.OnSelectQueryType(None)
		self.query_background_check.Disable()
		self.OnBackgroundCheck(None)
		self.is_ok()
		
#



	def OnBackgroundCheck(self, event):
		if DEBUG_LEVEL>1:
			print "QueryGui:OnBackgroundCheck()"
		self.query_background = self.query_background_check.GetValue()
#







	def OnSelectSep(self, event):
		if DEBUG_LEVEL>2:
			print "QueryGui: OnSelectSep()"
		if self.query_sep_custom_radio.GetValue():
			self.query_sep_custom_text.Enable()
			self.file_sep  = self.query_sep_custom_text.GetValue()
		elif self.query_sep_space_radio.GetValue():
			self.query_sep_custom_text.Disable()
			self.file_sep = " "
		elif self.query_sep_none_radio.GetValue():
			self.query_sep_custom_text.Disable()
			self.file_sep = None
		elif self.query_sep_tab_radio.GetValue():
			self.query_sep_custom_text.Disable()
			self.file_sep = "\t"
#





	def OnSelectQueryType(self, event):
		if DEBUG_LEVEL>1:
			print "QueryGui:OnSelectQueryType()"

		if self.query_file_type_box.GetSelection() >= 0:
			self.query_type = self.query_file_type_box.GetString(self.query_file_type_box.GetCurrentSelection())
		else:
			self.query_type = None
		
		if self.query_type == 'Pairs':
			self.query_sep_none_radio.Disable()
			if self.query_sep_none_radio.GetValue():
				self.query_sep_tab_radio.SetValue(True)
				self.file_sep = "\t"
				self.query_sep_none_radio.SetValue(False)
		if not self.query_type == 'Pairs':
			self.query_sep_none_radio.Enable()
#




	def OnCancel(self, event):
		self.Hide()
#







	def OnFileBrowse(self, event):
		dialog = wx.FileDialog(None, style = wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			self._set_file_name(dialog.GetPath())
#






	def _set_file_name(self, fn):
		self.filename = fn
		if fn == None:
			self.query_input_file_label.SetLabel(u"No query file selected.")
		else:
			self.query_input_file_label.SetLabel(os.path.basename(self.filename))
		self.is_ok()
#



	def is_ok(self):
		self.query_load_button.Enable()
		if self.filename == None:
			self.query_load_button.Disable()
#



	def freeze(self):
		self.query_pick_file_button.Disable()
		self.query_cancel_button.Disable()
		self.query_load_button.Disable()
		self.query_file_type_box.Disable()
		self.query_sep_tab_radio.Disable()
		self.query_sep_space_radio.Disable()
		self.query_sep_none_radio.Disable()
		self.query_sep_custom_radio.Disable()
		self.query_sep_custom_text.Disable()
#




	def unfreeze(self):
		self.query_pick_file_button.Enable()
		self.query_cancel_button.Enable()
		self.query_load_button.Enable()
		self.query_file_type_box.Enable()
		self.query_sep_tab_radio.Enable()
		self.query_sep_space_radio.Enable()
		self.query_sep_custom_radio.Enable()
		self.query_sep_custom_text.Enable()
		self.query_sep_none_radio.Enable()
		self._reset()
#





	def OnLoad(self, event):
		if DEBUG_LEVEL>1:
			print "query_load_gui: OnLoad()"
		self.freeze()

		if self.query_background:
			# process background
			#self.query_load_outcome_handle = self.real_parent.communication_thread.register_callback(self.real_parent.EVT_CUSTOM_LOAD_GO, self.OnLoadDone)
			#self.real_parent.gui2ssprocess_queue.put((WorkProcess.CMD_SET, WorkProcess.CMD_LOAD_GO, self.param_filename, {"ignore_part_of":self.param_GO_ignore_haspart, "ignore_regulates":self.param_GO_ignore_regulates}))
			pass
		else:
			self.current_query = []
			
			self.parent.text_query = ""
			h = open(self.filename, 'r')
			for i in h:
				i = i.rstrip("\n")
				i = i.rstrip("\r")
				i = i.rstrip(" ")
				
				if self.query_type == 'List':
					if self.file_sep == '' or self.file_sep == None:
						self.current_query.append(i)
						self.parent.text_query += (i+"\n")
					else:
						line = i.split(self.file_sep)
						for j in line:
							self.current_query.append(j)
							self.parent.text_query += (j+"\n")
				if self.query_type == 'Pairs':
					if self.file_sep == '' or self.file_sep == None:
						raise Exception
					else:
						line = i.split(self.file_sep)
						if len(line) >= 2:
							self.parent.text_query += (line[0]+"\t"+line[1]+"\n")
							self.current_query.append((line[0],line[1]))

			#print temp_query
			h.close()
			
			#print self.current_query
			#self._data_to_main()
			self.parent.OnLoadFromFileDone(None)
			self.unfreeze()
			self.OnCancel(None)

			#if data[2] == WorkProcess. RESULT_OK:
			#if DEBUG_LEVEL>1:
				#print "GO load outcome: Load successful."
				#print data
			#self.status = True
			#self.node_info = data[3]
			#self.edge_info = data[4]
			#self._data_to_main()
			#self.parent.update()
			#self.Hide()
			
		#else:
			#if DEBUG_LEVEL>2:
				#print "GO load outcome: Load Fail."
			#self.status = False
			#self.node_info = 0
			#self.edge_info = 0
			#self._data_to_main()
			#self.parent.update()

#
