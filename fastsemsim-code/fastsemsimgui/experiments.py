import sys
import os

import wx
import wx.lib.newevent


##################################################################################
#### start the main gui. Should be called from outside. Why? ####
##################################################################################

CustomEvent_Status, EVT_CUSTOM_STATUS = wx.lib.newevent.NewEvent()
CustomEvent_Run, EVT_CUSTOM_RUN = wx.lib.newevent.NewEvent()


class GOPanel(wx.Panel):
	def __init__( self, real_parent, parent, id, pos, size, style):
		wx.Panel.__init__ ( self, parent, id, pos, size, style)

		self.real_parent = real_parent
		GO_panel_sizer_1 = wx.BoxSizer( wx.VERTICAL )
		sbSizer3 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Source" ), wx.HORIZONTAL )
		self.GO_load_button = wx.Button( self, wx.ID_ANY, u"Load ...", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.GO_load_button.SetDefault() 
		sbSizer3.Add( self.GO_load_button, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		GO_panel_sizer_1.Add( sbSizer3, 0, wx.ALL|wx.EXPAND, 5 )
		
		self.Bind(wx.EVT_BUTTON, self.OnLoadButton, id=self.GO_load_button.GetId())

		self.real_parent.Bind(EVT_CUSTOM_RUN, self.func2)
		self.real_parent.Bind(EVT_CUSTOM_STATUS, self.func2)
#

	def OnLoadButton(self, event):
		wx.PostEvent(self.real_parent, CustomEvent_Run(data="prova"))
		wx.PostEvent(self.real_parent, CustomEvent_Status(data="prova"))
#

	def func2(self, event):
		print "GO: " + str(event.data)
		event.Skip()
#



class fastSemSimGui(wx.Frame):

	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "fastSemSimGui v.2 - Marco Mina", pos = wx.DefaultPosition, size = wx.Size( 652,471 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		

		fastSemSim_sizer_1 = wx.BoxSizer( wx.VERTICAL )
		self.fastSemSim_listbook = wx.Listbook( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LB_DEFAULT )
		fastSemSim_sizer_1.Add( self.fastSemSim_listbook, 1, wx.EXPAND |wx.ALL, 5 )
		self.SetSizer( fastSemSim_sizer_1 )
		self.Layout()
		self.Centre( wx.BOTH )

		self.GO_panel = GOPanel(self, self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.fastSemSim_listbook.AddPage( self.GO_panel, u"Gene Ontology", False )

		self.Bind(EVT_CUSTOM_STATUS, self.func)
#




	def func(self, event):
		print "main: " + str(event.data)
		event.Skip()
#



def _start():
	#multiprocessing.freeze_support()
	app = wx.App()
	window = fastSemSimGui(None)
	window.Centre()
	window.Show()  
	app.MainLoop()
#






if __name__=='__main__':
	_start()
#



