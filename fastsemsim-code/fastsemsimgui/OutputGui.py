import wx
import WorkProcess

DEBUG_LEVEL = 0

class OutputCtrlPanel(wx.Panel):
	
	param_filename = None
	
	def __init__( self, real_parent, parent, id, pos, size, style):
		wx.Panel.__init__ ( self, parent, id, pos, size, style)
		self.output_ctrl_panel = self # temporary workaround
		self.real_parent = real_parent
		bSizer93 = wx.BoxSizer( wx.VERTICAL )
		sbSizer511 = wx.StaticBoxSizer( wx.StaticBox( self.output_ctrl_panel, wx.ID_ANY, u"Output settings" ), wx.VERTICAL )
		sbSizer291 = wx.StaticBoxSizer( wx.StaticBox( self.output_ctrl_panel, wx.ID_ANY, u"Filter Options" ), wx.VERTICAL )
		gbSizer9 = wx.GridBagSizer( 0, 0 )
		gbSizer9.SetFlexibleDirection( wx.BOTH )
		gbSizer9.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		self.m_staticText1081 = wx.StaticText( self.output_ctrl_panel, wx.ID_ANY, u"Remove pairs with Semantic Similarity", wx.DefaultPosition, wx.Size( 200,-1 ), 0 )
		self.m_staticText1081.Wrap( -1 )
		gbSizer9.Add( self.m_staticText1081, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 1 ), wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		bSizer111 = wx.BoxSizer( wx.VERTICAL )
		self.output_filter_none_check = wx.CheckBox( self.output_ctrl_panel, wx.ID_ANY, u"= 'None'", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer111.Add( self.output_filter_none_check, 0, 0, 5 )
		self.output_filter_0_check = wx.CheckBox( self.output_ctrl_panel, wx.ID_ANY, u"= 0.0", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer111.Add( self.output_filter_0_check, 0, 0, 5 )
		bSizer112 = wx.BoxSizer( wx.HORIZONTAL )
		self.output_filter_less_check = wx.CheckBox( self.output_ctrl_panel, wx.ID_ANY, u"<", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer112.Add( self.output_filter_less_check, 0, 0, 5 )
		self.output_filter_less_text = wx.TextCtrl( self.output_ctrl_panel, wx.ID_ANY, u"0.0", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.output_filter_less_text.SetMaxSize( wx.Size( 70,-1 ) )
		bSizer112.Add( self.output_filter_less_text, 0, 0, 5 )
		bSizer111.Add( bSizer112, 1, wx.EXPAND, 5 )
		gbSizer9.Add( bSizer111, wx.GBPosition( 0, 1 ), wx.GBSpan( 1, 1 ), wx.EXPAND, 5 )
		sbSizer291.Add( gbSizer9, 1, wx.EXPAND, 5 )
		sbSizer511.Add( sbSizer291, 0, wx.ALL|wx.EXPAND, 5 )
		sbSizer52 = wx.StaticBoxSizer( wx.StaticBox( self.output_ctrl_panel, wx.ID_ANY, u"Redirect output to" ), wx.VERTICAL )
		self.output_to_box = wx.Choicebook( self.output_ctrl_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.CHB_DEFAULT )
		self.m_panel29 = wx.Panel( self.output_to_box, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer95 = wx.BoxSizer( wx.VERTICAL )
		self.m_panel29.SetSizer( bSizer95 )
		self.m_panel29.Layout()
		bSizer95.Fit( self.m_panel29 )
		self.output_to_box.AddPage( self.m_panel29, u"Output Window", True )
		self.m_panel30 = wx.Panel( self.output_to_box, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer94 = wx.BoxSizer( wx.VERTICAL )
		bSizer321 = wx.BoxSizer( wx.VERTICAL )
		bSizer221 = wx.BoxSizer( wx.HORIZONTAL )
		sbSizer341 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel30, wx.ID_ANY, u"Destination file" ), wx.HORIZONTAL )
		self.output_file_label = wx.StaticText( self.m_panel30, wx.ID_ANY, u"No file selected.", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.output_file_label.Wrap( -1 )
		sbSizer341.Add( self.output_file_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		bSizer221.Add( sbSizer341, 1, wx.LEFT|wx.RIGHT, 5 )
		bSizer201 = wx.BoxSizer( wx.HORIZONTAL )
		self.output_file_select_button = wx.Button( self.m_panel30, wx.ID_ANY, u"Select file...", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer201.Add( self.output_file_select_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		bSizer221.Add( bSizer201, 0, wx.ALIGN_RIGHT|wx.RIGHT|wx.TOP, 5 )
		bSizer321.Add( bSizer221, 0, wx.ALIGN_CENTER|wx.EXPAND|wx.TOP, 5 )
		bSizer551 = wx.BoxSizer( wx.VERTICAL )
		bSizer561 = wx.BoxSizer( wx.VERTICAL )
		bSizer551.Add( bSizer561, 1, wx.EXPAND, 5 )
		bSizer321.Add( bSizer551, 0, wx.ALL|wx.EXPAND, 5 )
		sbSizer2911 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel30, wx.ID_ANY, u"Options" ), wx.VERTICAL )
		gbSizer91 = wx.GridBagSizer( 0, 0 )
		gbSizer91.SetFlexibleDirection( wx.BOTH )
		gbSizer91.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		self.m_staticText3011 = wx.StaticText( self.m_panel30, wx.ID_ANY, u"Separator", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText3011.Wrap( -1 )
		gbSizer91.Add( self.m_staticText3011, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		bSizer4611 = wx.BoxSizer( wx.HORIZONTAL )
		self.output_sep_tab_radio = wx.RadioButton( self.m_panel30, wx.ID_ANY, u"[tab]  '\\t'", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer4611.Add( self.output_sep_tab_radio, 0, wx.ALL, 5 )
		self.output_sep_space_radio = wx.RadioButton( self.m_panel30, wx.ID_ANY, u"[space]", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer4611.Add( self.output_sep_space_radio, 0, wx.ALL, 5 )
		bSizer4711 = wx.BoxSizer( wx.HORIZONTAL )
		#self.output_sep_custom_radio = wx.RadioButton( self.m_panel30, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		#bSizer4711.Add( self.output_sep_custom_radio, 0, wx.ALL|wx.LEFT|wx.RIGHT, 5 )
		#self.output_sep_custom_text = wx.TextCtrl( self.m_panel30, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 35,-1 ), 0 )
		#bSizer4711.Add( self.output_sep_custom_text, 0, wx.ALL|wx.RIGHT, 5 )
		bSizer4611.Add( bSizer4711, 1, wx.EXPAND, 5 )
		gbSizer91.Add( bSizer4611, wx.GBPosition( 0, 1 ), wx.GBSpan( 1, 1 ), wx.EXPAND, 5 )
		sbSizer2911.Add( gbSizer91, 1, wx.EXPAND, 5 )
		bSizer321.Add( sbSizer2911, 1, wx.EXPAND, 5 )
		bSizer94.Add( bSizer321, 1, wx.EXPAND, 5 )
		self.m_panel30.SetSizer( bSizer94 )
		self.m_panel30.Layout()
		bSizer94.Fit( self.m_panel30 )
		self.output_to_box.AddPage( self.m_panel30, u"File", False )
		sbSizer52.Add( self.output_to_box, 0, wx.EXPAND |wx.ALL, 5 )
		sbSizer511.Add( sbSizer52, 0, wx.ALL|wx.EXPAND, 5 )
		bSizer93.Add( sbSizer511, 0, wx.ALL|wx.EXPAND, 5 )
		self.output_ctrl_panel.SetSizer( bSizer93 )
		self.output_ctrl_panel.Layout()
		bSizer93.Fit( self.output_ctrl_panel )
		
		self.Bind(wx.EVT_CHOICEBOOK_PAGE_CHANGED, self.OnSelectTo, id=self.output_to_box.GetId())
		self.Bind(wx.EVT_BUTTON, self.OnFileBrowse, id=self.output_file_select_button.GetId())

		self.Bind(wx.EVT_CHECKBOX, self.filter_what, id=self.output_filter_none_check.GetId())
		self.Bind(wx.EVT_CHECKBOX, self.filter_what, id=self.output_filter_0_check.GetId())
		self.Bind(wx.EVT_CHECKBOX, self.filter_what, id=self.output_filter_less_check.GetId())
		self.real_parent.Bind(wx.EVT_TEXT, self.filter_what, id=self.output_filter_less_text.GetId())
		self.Bind(wx.EVT_RADIOBUTTON, self.OnSep, id=self.output_sep_tab_radio.GetId())
		self.Bind(wx.EVT_RADIOBUTTON, self.OnSep, id=self.output_sep_space_radio.GetId())
		
		self._reset()
#






	def _reset(self):
		self.real_parent.params_output['to'] = WorkProcess.OUTPUT_TO_GUI
		self.param_filename = None
		self._update()
#






	def filter_what(self, event):
		if not 'filter' in self.real_parent.params_output:
			self.real_parent.params_output['filter'] = {}
		self.real_parent.params_output['filter']['None'] = self.output_filter_none_check.GetValue()
		self.real_parent.params_output['filter']['0'] = self.output_filter_0_check.GetValue()
		self.real_parent.params_output['filter']['less'] = self.output_filter_less_check.GetValue()
		try:
			if str(self.output_filter_less_text.GetValue()) == "":
				self.real_parent.params_output['filter']['value'] = None
			else:
				self.real_parent.params_output['filter']['value'] = float(self.output_filter_less_text.GetValue())
		except:
			self.output_filter_less_text.SetValue(str(""))
			self.filter_what(None)
#






	def _freeze(self):
		self.Disable()
#







	def OnSep(self,event):
		#print "OnSep"
		if self.output_sep_tab_radio.GetValue():
			self.real_parent.params_output['params']['sep'] = "\t"
		elif self.output_sep_space_radio.GetValue():
			self.real_parent.params_output['params']['sep'] = " "
#






	def _unfreeze(self):
		self.Enable()
#






	def OnFileBrowse(self, event):
		dialog = wx.FileDialog(None, style = wx.OPEN)
		if dialog.ShowModal() == wx.ID_OK:
			self.param_filename = dialog.GetPath()
			self._set_file_name()
#





	def _set_file_name(self):
		if self.param_filename == None:
			self.output_file_label.SetLabel(u"No file selected.")
		else:
			self.output_file_label.SetLabel(self.param_filename)
			self.real_parent.params_output['filename'] = self.param_filename
#






	def OnSelectTo(self, event):
		if DEBUG_LEVEL>0:
			print "OutputPanel: OnSelectTo()"
		self.real_parent.params_output['to'] = self.output_to_box.GetSelection()
#






	def _update(self):
		if DEBUG_LEVEL>0:
			print "OutputPanel: _update()"
		if self.real_parent.params_output['to'] == WorkProcess.OUTPUT_TO_GUI:
			self.output_to_box.SetSelection(0)
		elif self.real_parent.params_output['to'] == WorkProcess.OUTPUT_TO_FILE:
			self.output_to_box.SetSelection(1)
			
		if 'filename' in self.real_parent.params_output:
			self.param_filename = self.real_parent.params_output['filename']
			self._set_file_name()
			
		if 'filter' in self.real_parent.params_output:
			if 'None' in self.real_parent.params_output['filter']:
				self.output_filter_none_check.SetValue(bool(self.real_parent.params_output['filter']['None']))
			if '0' in self.real_parent.params_output['filter']:
				self.output_filter_0_check.SetValue(bool(self.real_parent.params_output['filter']['0']))
			if 'less' in self.real_parent.params_output['filter']:
				self.output_filter_less_check.SetValue(bool(self.real_parent.params_output['filter']['less']))
			if 'value' in self.real_parent.params_output['filter']:
				self.output_filter_less_text.SetValue(str(self.real_parent.params_output['filter']['value']))
#

#
			
			
		#if 'mixing_strategy' in self.real_parent.params_SS:
			#self.SS_mix_box.SetStringSelection(self.real_parent.params_SS['mixing_strategy'])
		#if 'measure' in self.real_parent.params_SS:
			#self.SS_measure_box.SetStringSelection(self.real_parent.params_SS['measure'])
		#if 'ontology' in self.real_parent.params_SS:
			#if self.real_parent.params_SS['ontology'] =="CC":
				#self.SS_GO_CC_radio.SetValue(1)
			#elif self.real_parent.params_SS['ontology'] == 'BP':
				#self.SS_GO_BP_radio.SetValue(1)
			#elif self.real_parent.params_SS['ontology'] == 'MF':
				#self.SS_GO_MF_radio.SetValue(1)
			#else:
				#self.SS_GO_MF_radio.SetValue(0)
				#self.SS_GO_BP_radio.SetValue(0)
				#self.SS_GO_CC_radio.SetValue(0)
#


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


class OutputWindow(wx.Frame):
	def __init__( self, parent):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "Output", pos = wx.DefaultPosition, size = wx.Size( 420,400 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL| wx.RESIZE_BORDER )
		self.real_parent = parent
		self.parent = parent
		
		self.SetSizeHintsSz( wx.Size( 330,-1 ), wx.DefaultSize )
		#self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		bSizer62 = wx.BoxSizer( wx.VERTICAL)
		
		sbSizer35 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, ""), wx.HORIZONTAL|wx.EXPAND )
		self.output_text = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_MULTILINE|wx.TE_READONLY)
		self.output_text.SetMinSize( wx.Size( 300, 400 ) )
		sbSizer35.Add( self.output_text, 0, wx.ALL|wx.EXPAND, 5 )
		bSizer62.Add( sbSizer35, 1, wx.ALL|wx.EXPAND, 5 )
		
		self.SetSizer( bSizer62 )
		self.Layout()
		bSizer62.Fit( self )
		#self.Show(True)
#
