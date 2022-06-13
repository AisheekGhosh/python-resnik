import wx

class ControlsPanel(wx.Panel):
	def __init__( self, real_parent, parent, id, pos, size, style):
		wx.Panel.__init__ ( self, parent, id, pos, size, style)
		self.parent = real_parent
		self.controls_panel = self # temporary workaround
		
		self.ok_bitmap = wx.Bitmap(self.parent.Ok_pic)
		self.warning_bitmap = wx.Bitmap(self.parent.Warning_pic)
		
		bSizer62 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer35 = wx.StaticBoxSizer( wx.StaticBox( self.controls_panel, wx.ID_ANY, u"Log" ), wx.VERTICAL )
		
		self.controls_log_text = wx.TextCtrl( self.controls_panel, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( -1,150 ), wx.TE_MULTILINE|wx.TE_READONLY )
		sbSizer35.Add( self.controls_log_text, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer62.Add( sbSizer35, 0, wx.ALL|wx.EXPAND, 5 )

		bSizer562 = wx.BoxSizer( wx.VERTICAL )
		bSizer57 = wx.BoxSizer( wx.HORIZONTAL )
		bSizer63 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.controls_start_button = wx.Button( self.controls_panel, wx.ID_ANY, u"Start", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer63.Add( self.controls_start_button, 0, wx.ALL, 5 )
		
		self.controls_stop_button = wx.Button( self.controls_panel, wx.ID_ANY, u"Stop", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer63.Add( self.controls_stop_button, 0, wx.ALL, 5 )
		
		bSizer562.Add( bSizer63, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		sbSizer321 = wx.StaticBoxSizer( wx.StaticBox( self.controls_panel, wx.ID_ANY, u"Tips" ), wx.HORIZONTAL )
		
		self.controls_tip_label = wx.StaticText( self.controls_panel, wx.ID_ANY, u"Suggested procedure:\n\n1) Load a Gene Ontology using the first panel.\n\n2)Load an Annotation Corpus using the first two panels.\n3)Select a Semantic Similarity measure using the Semantic Similarity panel.\n4) Specify a query using the Query panel.\n5) Set the output parameters in Output Control.\n6) Use this panel to start the computation.", wx.DefaultPosition, wx.Size( 250,-1 ), 0 )
		self.controls_tip_label.Wrap( -1 )

		sbSizer321.Add( self.controls_tip_label, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		sbSizer331 = wx.StaticBoxSizer( wx.StaticBox( self.controls_panel, wx.ID_ANY, u"Status" ), wx.VERTICAL )
		
		fgSizer4 = wx.FlexGridSizer( 5, 2, 2, 2 )
		fgSizer4.SetFlexibleDirection( wx.BOTH )
		fgSizer4.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.m_staticText481 = wx.StaticText( self.controls_panel, wx.ID_ANY, u"Gene Ontology", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText481.Wrap( -1 )
		fgSizer4.Add( self.m_staticText481, 0, wx.ALL, 5 )

		self.go_status_bitmap = wx.StaticBitmap( self.controls_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.go_status_bitmap, 0, wx.ALL, 5 )
		
		self.m_staticText49 = wx.StaticText( self.controls_panel, wx.ID_ANY, u"Annotation Corpus", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText49.Wrap( -1 )
		fgSizer4.Add( self.m_staticText49, 0, wx.ALL, 5 )
		
		self.ac_status_bitmap = wx.StaticBitmap( self.controls_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.ac_status_bitmap, 0, wx.ALL, 5 )
		
		self.m_staticText50 = wx.StaticText( self.controls_panel, wx.ID_ANY, u"SS Parameters", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText50.Wrap( -1 )
		fgSizer4.Add( self.m_staticText50, 0, wx.ALL, 5 )
		
		self.ss_status_bitmap = wx.StaticBitmap( self.controls_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.ss_status_bitmap, 0, wx.ALL, 5 )
		
		self.m_staticText512 = wx.StaticText( self.controls_panel, wx.ID_ANY, u"Query", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText512.Wrap( -1 )
		fgSizer4.Add( self.m_staticText512, 0, wx.ALL, 5 )
		
		self.query_status_bitmap = wx.StaticBitmap( self.controls_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.query_status_bitmap, 0, wx.ALL, 5 )
		
		self.m_staticText521 = wx.StaticText( self.controls_panel, wx.ID_ANY, u"Output Parameters", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText521.Wrap( -1 )
		fgSizer4.Add( self.m_staticText521, 0, wx.ALL, 5 )
		
		self.output_status_bitmap = wx.StaticBitmap( self.controls_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer4.Add( self.output_status_bitmap, 0, wx.ALL, 5 )

		sbSizer331.Add( fgSizer4, 1, wx.EXPAND, 5 )
		
		bSizer57.Add( sbSizer321, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		bSizer57.Add( sbSizer331, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer562.Add( bSizer57, 1, wx.EXPAND, 5 )
		
		bSizer62.Add( bSizer562, 1, wx.EXPAND, 5 )

		self.controls_panel.SetSizer( bSizer62 )
		self.controls_panel.Layout()
		bSizer62.Fit( self.controls_panel )

		self.parent.Bind(wx.EVT_BUTTON, self.OnStartCmd, id=self.controls_start_button.GetId())
		self.parent.Bind(wx.EVT_BUTTON, self.OnStopCmd, id=self.controls_stop_button.GetId())
		
		self.SetSSCtrlStatus(False)
		self.SetGOStatus(False)
		self.SetAcStatus(False)
		self.SetQueryStatus(False)
		self.SetOutputCtrlStatus(False)
		
		
#










	def OnStopCmd(self, event):
		self.parent.stop()
#










	def OnStartCmd(self, event):
		if self.parent.start():
			#self.controls_start_button.SetLabel(u"Pause")
			pass
#






	def setBitmap(self, handle, value):
		if value:
			handle.SetBitmap(self.ok_bitmap)
		else:
			handle.SetBitmap(self.warning_bitmap)
#

	def SetGOStatus(self, status):
		self.go_ok = status
		self.setBitmap(self.go_status_bitmap, status)
		
	def SetAcStatus(self, status):
		self.ac_ok = status
		self.setBitmap(self.ac_status_bitmap, status)
		
	def SetQueryStatus(self, status):
		self.query_ok = status
		self.setBitmap(self.query_status_bitmap, status)
		#self.activateGoCmd()
		
	def SetOutputCtrlStatus(self, status):
		self.outputctrl_ok = status
		#if self.show_pics:
		self.setBitmap(self.output_status_bitmap, status)
		#self.activateGoCmd()

	def SetSSCtrlStatus(self, status):
		#if self.show_pics:
		self.setBitmap(self.ss_status_bitmap, status)
		#self.activateGoCmd()
#
