import wx
# from fastSemSim.SemSim import SemSimMeasures
import fastsemsim.SemSim

DEBUG_LEVEL = 0

class SSPanel(wx.Panel):
	def __init__( self, real_parent, parent, id, pos, size, style):
		wx.Panel.__init__ ( self, parent, id, pos, size, style)
		self.real_parent = real_parent
		self.SS_panel = self # temporary workaround

		bSizer25 = wx.BoxSizer( wx.VERTICAL )
		bSizer27 = wx.BoxSizer( wx.HORIZONTAL )
		bSizer26 = wx.BoxSizer( wx.VERTICAL )
		sbSizer15 = wx.StaticBoxSizer( wx.StaticBox( self.SS_panel, wx.ID_ANY, u"Semantic Measure" ), wx.HORIZONTAL )
		self.SS_measure_box = wx.ComboBox( self.SS_panel, wx.ID_ANY, value="", pos=wx.DefaultPosition, size=wx.DefaultSize, style=wx.CB_READONLY)

		sbSizer15.Add( self.SS_measure_box, 0, wx.ALL, 5 )
		#self.SS_measure_settings_button = wx.BitmapButton( self.SS_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, wx.BU_AUTODRAW )
		#sbSizer15.Add( self.SS_measure_settings_button, 0, wx.ALL, 5 )
		bSizer26.Add( sbSizer15, 1, wx.EXPAND|wx.RIGHT, 5 )
		sbSizer16 = wx.StaticBoxSizer( wx.StaticBox( self.SS_panel, wx.ID_ANY, u"Mixing Strategy" ), wx.HORIZONTAL )

		self.SS_mix_box = wx.ComboBox( self.SS_panel, wx.ID_ANY, "", wx.DefaultPosition, wx.DefaultSize, style = wx.CB_READONLY)
		sbSizer16.Add( self.SS_mix_box, 0, wx.ALL, 5 )
		#self.SS_mix_settings_button = wx.BitmapButton( self.SS_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, wx.BU_AUTODRAW )
		#sbSizer16.Add( self.SS_mix_settings_button, 0, wx.ALL, 5 )
		bSizer26.Add( sbSizer16, 1, wx.EXPAND|wx.RIGHT|wx.TOP, 5 )
		bSizer27.Add( bSizer26, 0, 0, 5 )
		sbSizer17 = wx.StaticBoxSizer( wx.StaticBox( self.SS_panel, wx.ID_ANY, u"Ontology" ), wx.VERTICAL )
		self.SS_GO_MF_radio = wx.RadioButton( self.SS_panel, wx.ID_ANY, u"Molecular Function (MF)", wx.DefaultPosition, wx.DefaultSize, 0 )
		sbSizer17.Add( self.SS_GO_MF_radio, 0, wx.ALL, 5 )
		self.SS_GO_CC_radio = wx.RadioButton( self.SS_panel, wx.ID_ANY, u"Cellular Component (CC)", wx.DefaultPosition, wx.DefaultSize, 0 )
		sbSizer17.Add( self.SS_GO_CC_radio, 0, wx.ALL, 5 )
		self.SS_GO_BP_radio = wx.RadioButton( self.SS_panel, wx.ID_ANY, u"Biological Process (BP)", wx.DefaultPosition, wx.DefaultSize, 0 )
		sbSizer17.Add( self.SS_GO_BP_radio, 0, wx.ALL, 5 )
		bSizer27.Add( sbSizer17, 0, wx.LEFT, 5 )
		bSizer25.Add( bSizer27, 1, wx.ALL|wx.EXPAND, 5 )
		self.SS_panel.SetSizer( bSizer25 )
		self.SS_panel.Layout()
		bSizer25.Fit( self.SS_panel )
		
		self.real_parent.Bind(wx.EVT_COMBOBOX, self.OnSelectSS, id=self.SS_measure_box.GetId())
		self.real_parent.Bind(wx.EVT_COMBOBOX, self.OnSelectMS, id=self.SS_mix_box.GetId())
		self.real_parent.Bind(wx.EVT_RADIOBUTTON, self.OnSelectGO, id=self.SS_GO_MF_radio.GetId())
		self.real_parent.Bind(wx.EVT_RADIOBUTTON, self.OnSelectGO, id=self.SS_GO_BP_radio.GetId())
		self.real_parent.Bind(wx.EVT_RADIOBUTTON, self.OnSelectGO, id=self.SS_GO_CC_radio.GetId())
		self._reset()
#



	def _reset(self):
		if DEBUG_LEVEL>0:
			print "SSPanel: reset()"
		for i in fastsemsim.SemSim.term_SemSim_measures:
			self.SS_measure_box.Append(i)
		for i in fastsemsim.SemSim.mix_strategies:
			self.SS_mix_box.Append(i)
		self.SS_mix_box.SetSelection(0)
		self.SS_measure_box.SetSelection(0)
		self.SS_GO_MF_radio.SetValue(True)
		self.OnSelectGO(None)
		self.OnSelectSS(None)
		self.OnSelectMS(None)
		self.SS_to_send = True
		self._update()
#





	def _update(self):
		if DEBUG_LEVEL>0:
			print "SSPanel: _update()"
		if 'mixing_strategy' in self.real_parent.params_SS:
			self.SS_mix_box.SetStringSelection(self.real_parent.params_SS['mixing_strategy'])
			self.OnSelectMS(None)
		if 'measure' in self.real_parent.params_SS:
			self.SS_measure_box.SetStringSelection(self.real_parent.params_SS['measure'])
			self.OnSelectSS(None)
		if 'ontology' in self.real_parent.params_SS:
			if self.real_parent.params_SS['ontology'] =="CC":
				self.SS_GO_CC_radio.SetValue(1)
			elif self.real_parent.params_SS['ontology'] == 'BP':
				self.SS_GO_BP_radio.SetValue(1)
			elif self.real_parent.params_SS['ontology'] == 'MF':
				self.SS_GO_MF_radio.SetValue(1)
			else:
				self.SS_GO_MF_radio.SetValue(0)
				self.SS_GO_BP_radio.SetValue(0)
				self.SS_GO_CC_radio.SetValue(0)
#





	def _freeze(self):
		self.Disable()
#






	def _unfreeze(self):
		self.Enable()
#







	def OnSelectGO(self, event):
		if DEBUG_LEVEL>0:
			print "SSPanel: OnSelectGO()"
		if self.SS_GO_CC_radio.GetValue():
			self.real_parent.params_SS['ontology'] = "CC"
		elif self.SS_GO_BP_radio.GetValue():
			self.real_parent.params_SS['ontology'] = "BP"
		elif self.SS_GO_MF_radio.GetValue():
			self.real_parent.params_SS['ontology'] = "MF"
		else:
			self.real_parent.params_SS['ontology'] = None
		self.SS_to_send = True
#







	def OnSelectMS(self, event):
		if DEBUG_LEVEL>0:
			print "SSPanel: OnSelectMS()"
		c = self.SS_mix_box.GetCurrentSelection()
		#print c
		if c >= 0:
			self.real_parent.params_SS['mixing_strategy'] = self.SS_mix_box.GetString(self.SS_mix_box.GetCurrentSelection())
			#print self.real_parent.mixingstrategy
			#self.ok = True
			#if self.real_parent.ssmeasure is None:
				#self.real_parent.SetOperationOk(False)
			#else:
				#self.real_parent.SetOperationOk(True)
		else:
			self.real_parent.params_SS['mixing_strategy'] = None
		self.SS_to_send = True
#








	def OnSelectSS(self, event):
		if DEBUG_LEVEL>0:
			print "SSPanel: OnSelectSS()"
		c = self.SS_measure_box.GetCurrentSelection()
		if c >= 0:
			self.real_parent.params_SS['measure'] = self.SS_measure_box.GetString(self.SS_measure_box.GetCurrentSelection())
			#self.real_parent.ssmeasure = SemSimMeasures.SemSimMeasures[self.availablemeasures[self.ss.GetSelection()]]
			#self.real_parent.update_ssobject = True
			if fastsemsim.SemSim.term_SemSim_measures[self.real_parent.params_SS['measure']][1]:
				self.SS_mix_box.Enable()
				#if self.real_parent.SS_mixing_strategy == None:
					#self.real_parent.SetOperationOk(False)
					#pass
			else:
				self.SS_mix_box.Disable()
				#self.real_parent.SetOperationOk(True)
		else:
			self.real_parent.params_SS['measure'] = None
		self.SS_to_send = True
#
