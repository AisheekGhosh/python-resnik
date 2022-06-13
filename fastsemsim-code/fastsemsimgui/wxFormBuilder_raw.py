# -*- coding: utf-8 -*- 

###########################################################################
## Python code generated with wxFormBuilder (version Sep  8 2010)
## http://www.wxformbuilder.org/
##
## PLEASE DO "NOT" EDIT THIS FILE!
###########################################################################

import wx

###########################################################################
## Class fastSemSim_gui
###########################################################################

class fastSemSim_gui ( wx.Frame ):
	
	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.DefaultPosition, size = wx.Size( 652,471 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		fastSemSim_sizer_1 = wx.BoxSizer( wx.VERTICAL )
		
		self.fastSemSim_listbook = wx.Listbook( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LB_DEFAULT )
		self.GO_panel = wx.Panel( self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
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
		
		sbSizer5.Add( gSizer2, 1, wx.ALIGN_CENTER|wx.SHAPED, 2 )
		
		bSizer5.Add( sbSizer5, 0, wx.ALL, 5 )
		
		sbSizer6 = wx.StaticBoxSizer( wx.StaticBox( self.GO_panel, wx.ID_ANY, u"Additional information" ), wx.VERTICAL )
		
		bSizer51 = wx.BoxSizer( wx.VERTICAL )
		
		self.m_staticText8 = wx.StaticText( self.GO_panel, wx.ID_ANY, u"what has been accepted?", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText8.Wrap( -1 )
		bSizer51.Add( self.m_staticText8, 0, wx.ALL, 5 )
		
		sbSizer6.Add( bSizer51, 1, wx.EXPAND, 5 )
		
		bSizer5.Add( sbSizer6, 0, wx.ALL, 5 )
		
		GO_panel_sizer_1.Add( bSizer5, 0, wx.EXPAND, 5 )
		
		bSizer4 = wx.BoxSizer( wx.HORIZONTAL )
		
		GO_panel_sizer_1.Add( bSizer4, 0, wx.ALIGN_RIGHT, 5 )
		
		self.GO_panel.SetSizer( GO_panel_sizer_1 )
		self.GO_panel.Layout()
		GO_panel_sizer_1.Fit( self.GO_panel )
		self.fastSemSim_listbook.AddPage( self.GO_panel, u"Gene Ontology", False )
		self.AC_panel = wx.Panel( self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
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
		
		self.m_staticText31 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText31.Wrap( -1 )
		gSizer21.Add( self.m_staticText31, 0, wx.ALL, 5 )
		
		self.m_staticText41 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"Edges", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText41.Wrap( -1 )
		gSizer21.Add( self.m_staticText41, 0, wx.ALL, 5 )
		
		self.m_staticText51 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText51.Wrap( -1 )
		gSizer21.Add( self.m_staticText51, 0, wx.ALL, 5 )
		
		self.m_staticText61 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"Categories", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText61.Wrap( -1 )
		gSizer21.Add( self.m_staticText61, 0, wx.ALL, 5 )
		
		self.m_staticText71 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText71.Wrap( -1 )
		gSizer21.Add( self.m_staticText71, 0, wx.ALL, 5 )
		
		sbSizer51.Add( gSizer21, 1, wx.ALIGN_CENTER|wx.SHAPED, 2 )
		
		bSizer52.Add( sbSizer51, 0, wx.ALL, 5 )
		
		sbSizer61 = wx.StaticBoxSizer( wx.StaticBox( self.AC_panel, wx.ID_ANY, u"Additional information" ), wx.VERTICAL )
		
		bSizer511 = wx.BoxSizer( wx.VERTICAL )
		
		gbSizer1 = wx.GridBagSizer( 0, 0 )
		gbSizer1.SetFlexibleDirection( wx.BOTH )
		gbSizer1.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.m_staticText18 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"File type", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText18.Wrap( -1 )
		gbSizer1.Add( self.m_staticText18, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.m_staticText19 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText19.Wrap( -1 )
		gbSizer1.Add( self.m_staticText19, wx.GBPosition( 0, 1 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.m_staticText20 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"MyLabel", wx.Point( -1,-1 ), wx.DefaultSize, 0 )
		self.m_staticText20.Wrap( -1 )
		gbSizer1.Add( self.m_staticText20, wx.GBPosition( 1, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.m_staticText211 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText211.Wrap( -1 )
		gbSizer1.Add( self.m_staticText211, wx.GBPosition( 1, 1 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.m_staticText22 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"MyLabel", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText22.Wrap( -1 )
		gbSizer1.Add( self.m_staticText22, wx.GBPosition( 2, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.m_staticText23 = wx.StaticText( self.AC_panel, wx.ID_ANY, u"-", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText23.Wrap( -1 )
		gbSizer1.Add( self.m_staticText23, wx.GBPosition( 2, 1 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		bSizer511.Add( gbSizer1, 1, wx.EXPAND, 5 )
		
		sbSizer61.Add( bSizer511, 1, wx.EXPAND, 5 )
		
		bSizer52.Add( sbSizer61, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer31.Add( bSizer52, 0, wx.EXPAND, 5 )
		
		self.AC_panel.SetSizer( bSizer31 )
		self.AC_panel.Layout()
		bSizer31.Fit( self.AC_panel )
		self.fastSemSim_listbook.AddPage( self.AC_panel, u"Annotation Corpus", False )
		self.SS_panel = wx.Panel( self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer25 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer27 = wx.BoxSizer( wx.HORIZONTAL )
		
		bSizer26 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer15 = wx.StaticBoxSizer( wx.StaticBox( self.SS_panel, wx.ID_ANY, u"Semantic Measure" ), wx.HORIZONTAL )
		
		SS_measure_boxChoices = []
		self.SS_measure_box = wx.ComboBox( self.SS_panel, wx.ID_ANY, u"Resnik", wx.DefaultPosition, wx.DefaultSize, SS_measure_boxChoices, 0 )
		sbSizer15.Add( self.SS_measure_box, 0, wx.ALL, 5 )
		
		self.SS_measure_settings_button = wx.BitmapButton( self.SS_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, wx.BU_AUTODRAW )
		sbSizer15.Add( self.SS_measure_settings_button, 0, wx.ALL, 5 )
		
		bSizer26.Add( sbSizer15, 1, wx.EXPAND|wx.RIGHT, 5 )
		
		sbSizer16 = wx.StaticBoxSizer( wx.StaticBox( self.SS_panel, wx.ID_ANY, u"Mixing Strategy" ), wx.HORIZONTAL )
		
		SS_mix_boxChoices = []
		self.SS_mix_box = wx.ComboBox( self.SS_panel, wx.ID_ANY, u"BMA", wx.DefaultPosition, wx.DefaultSize, SS_mix_boxChoices, 0 )
		sbSizer16.Add( self.SS_mix_box, 0, wx.ALL, 5 )
		
		self.SS_mix_settings_button = wx.BitmapButton( self.SS_panel, wx.ID_ANY, wx.NullBitmap, wx.DefaultPosition, wx.DefaultSize, wx.BU_AUTODRAW )
		sbSizer16.Add( self.SS_mix_settings_button, 0, wx.ALL, 5 )
		
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
		self.fastSemSim_listbook.AddPage( self.SS_panel, u"Semantic Similarity", False )
		self.query_panel = wx.Panel( self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer50 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer36 = wx.StaticBoxSizer( wx.StaticBox( self.query_panel, wx.ID_ANY, u"Input method" ), wx.VERTICAL )
		
		self.query_input_method_box = wx.Choicebook( self.query_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.CHB_DEFAULT )
		self.m_panel22 = wx.Panel( self.query_input_method_box, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer512 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer312 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel22, wx.ID_ANY, u"Input field" ), wx.VERTICAL )
		
		self.query_field_text = wx.TextCtrl( self.m_panel22, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 250,250 ), wx.TE_MULTILINE )
		self.query_field_text.SetMinSize( wx.Size( 250,250 ) )
		
		sbSizer312.Add( self.query_field_text, 0, wx.ALL, 5 )
		
		bSizer512.Add( sbSizer312, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer54 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer32 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel22, wx.ID_ANY, u"Input format" ), wx.VERTICAL )
		
		query_field_type_boxChoices = []
		self.query_field_type_box = wx.ComboBox( self.m_panel22, wx.ID_ANY, u"List", wx.DefaultPosition, wx.Size( 120,-1 ), query_field_type_boxChoices, wx.CB_READONLY )
		sbSizer32.Add( self.query_field_type_box, 0, wx.ALL, 5 )
		
		bSizer54.Add( sbSizer32, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		sbSizer33 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel22, wx.ID_ANY, u"Help" ), wx.VERTICAL )
		
		self.m_staticText53 = wx.StaticText( self.m_panel22, wx.ID_ANY, u"You can either specify a list of object or a set of pairs of objects. In the former case, the semantic similarity will be evaluated between each pair of objects in the list. In the latter, each pair will be considered", wx.DefaultPosition, wx.Size( 150,-1 ), 0 )
		self.m_staticText53.Wrap( -1 )
		sbSizer33.Add( self.m_staticText53, 0, wx.ALL, 5 )
		
		bSizer54.Add( sbSizer33, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL|wx.EXPAND, 5 )
		
		self.query_reset_button = wx.Button( self.m_panel22, wx.ID_ANY, u"Reset", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer54.Add( self.query_reset_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer512.Add( bSizer54, 1, wx.EXPAND, 5 )
		
		self.m_panel22.SetSizer( bSizer512 )
		self.m_panel22.Layout()
		bSizer512.Fit( self.m_panel22 )
		self.query_input_method_box.AddPage( self.m_panel22, u"Manual input", False )
		self.m_panel23 = wx.Panel( self.query_input_method_box, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer521 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer32 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer22 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer34 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel23, wx.ID_ANY, u"Source file" ), wx.HORIZONTAL )
		
		self.query_input_file_label = wx.StaticText( self.m_panel23, wx.ID_ANY, u"No file selected.", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.query_input_file_label.Wrap( -1 )
		sbSizer34.Add( self.query_input_file_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		bSizer22.Add( sbSizer34, 1, wx.LEFT|wx.RIGHT, 5 )
		
		bSizer20 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.query_load_button = wx.Button( self.m_panel23, wx.ID_ANY, u"Select file...", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer20.Add( self.query_load_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer22.Add( bSizer20, 0, wx.ALIGN_RIGHT|wx.RIGHT|wx.TOP, 5 )
		
		bSizer32.Add( bSizer22, 0, wx.ALIGN_CENTER|wx.EXPAND|wx.TOP, 5 )
		
		bSizer55 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer56 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer29 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel23, wx.ID_ANY, u"Options" ), wx.VERTICAL )
		
		fgSizer8 = wx.FlexGridSizer( 3, 3, 15, 5 )
		fgSizer8.SetFlexibleDirection( wx.BOTH )
		fgSizer8.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.m_staticText30 = wx.StaticText( self.m_panel23, wx.ID_ANY, u"Separator", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText30.Wrap( -1 )
		fgSizer8.Add( self.m_staticText30, 0, wx.ALL, 5 )
		
		bSizer46 = wx.BoxSizer( wx.VERTICAL )
		
		self.query_sep_tab_radio = wx.RadioButton( self.m_panel23, wx.ID_ANY, u"[tab]  '\\t'", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer46.Add( self.query_sep_tab_radio, 0, wx.LEFT|wx.RIGHT|wx.TOP, 5 )
		
		self.query_sep_space_radio = wx.RadioButton( self.m_panel23, wx.ID_ANY, u"[space]", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer46.Add( self.query_sep_space_radio, 0, wx.LEFT|wx.RIGHT, 5 )
		
		bSizer47 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.query_sep_custom_radio = wx.RadioButton( self.m_panel23, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer47.Add( self.query_sep_custom_radio, 0, wx.LEFT|wx.RIGHT, 5 )
		
		self.query_sep_custom_text = wx.TextCtrl( self.m_panel23, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 35,-1 ), 0 )
		bSizer47.Add( self.query_sep_custom_text, 0, wx.RIGHT, 5 )
		
		bSizer46.Add( bSizer47, 1, wx.EXPAND, 5 )
		
		fgSizer8.Add( bSizer46, 1, wx.EXPAND, 5 )
		
		self.m_staticText48 = wx.StaticText( self.m_panel23, wx.ID_ANY, u"Character used to separate fields within each row.", wx.DefaultPosition, wx.Size( 200,-1 ), 0 )
		self.m_staticText48.Wrap( -1 )
		fgSizer8.Add( self.m_staticText48, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		self.m_staticText511 = wx.StaticText( self.m_panel23, wx.ID_ANY, u"Input type", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText511.Wrap( -1 )
		fgSizer8.Add( self.m_staticText511, 0, wx.ALL, 5 )
		
		query_file_type_boxChoices = [ u"Pairs", u"List" ]
		self.query_file_type_box = wx.ComboBox( self.m_panel23, wx.ID_ANY, u"Pairs", wx.DefaultPosition, wx.Size( 100,-1 ), query_file_type_boxChoices, wx.CB_READONLY )
		fgSizer8.Add( self.query_file_type_box, 0, wx.ALL, 5 )
		
		self.m_staticText52 = wx.StaticText( self.m_panel23, wx.ID_ANY, u"Select whether the input is a list of objects or a set of pairs first field of each row is an object or a GO Term.", wx.DefaultPosition, wx.Size( 200,-1 ), 0 )
		self.m_staticText52.Wrap( -1 )
		fgSizer8.Add( self.m_staticText52, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		sbSizer29.Add( fgSizer8, 1, wx.EXPAND, 5 )
		
		bSizer56.Add( sbSizer29, 1, wx.EXPAND, 5 )
		
		bSizer55.Add( bSizer56, 1, wx.EXPAND, 5 )
		
		bSizer32.Add( bSizer55, 0, wx.ALL|wx.EXPAND, 5 )
		
		sbSizer18 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel23, wx.ID_ANY, u"Help" ), wx.VERTICAL )
		
		self.m_staticText311 = wx.StaticText( self.m_panel23, wx.ID_ANY, u"Select a file", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		self.m_staticText311.Wrap( -1 )
		self.m_staticText311.SetMaxSize( wx.Size( 500,-1 ) )
		
		sbSizer18.Add( self.m_staticText311, 0, wx.ALL, 5 )
		
		bSizer32.Add( sbSizer18, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer521.Add( bSizer32, 1, wx.EXPAND, 5 )
		
		self.m_panel23.SetSizer( bSizer521 )
		self.m_panel23.Layout()
		bSizer521.Fit( self.m_panel23 )
		self.query_input_method_box.AddPage( self.m_panel23, u"Load query from file", True )
		self.m_panel24 = wx.Panel( self.query_input_method_box, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer53 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer181 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel24, wx.ID_ANY, u"Help" ), wx.VERTICAL )
		
		self.m_staticText3111 = wx.StaticText( self.m_panel24, wx.ID_ANY, u"All the objects in the annotation corpus will be used as input list. Pairwise Semantic Similarity will be evaluated on such list", wx.DefaultPosition, wx.Size( 400,-1 ), 0 )
		self.m_staticText3111.Wrap( -1 )
		self.m_staticText3111.SetMaxSize( wx.Size( 500,-1 ) )
		
		sbSizer181.Add( self.m_staticText3111, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer53.Add( sbSizer181, 0, wx.EXPAND, 5 )
		
		self.m_panel24.SetSizer( bSizer53 )
		self.m_panel24.Layout()
		bSizer53.Fit( self.m_panel24 )
		self.query_input_method_box.AddPage( self.m_panel24, u"Annotation Corpus", False )
		sbSizer36.Add( self.query_input_method_box, 0, wx.EXPAND |wx.ALL, 5 )
		
		bSizer50.Add( sbSizer36, 0, wx.ALL|wx.EXPAND, 5 )
		
		self.query_panel.SetSizer( bSizer50 )
		self.query_panel.Layout()
		bSizer50.Fit( self.query_panel )
		self.fastSemSim_listbook.AddPage( self.query_panel, u"Query", False )
		self.output_ctrl_panel = wx.Panel( self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer93 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer511 = wx.StaticBoxSizer( wx.StaticBox( self.output_ctrl_panel, wx.ID_ANY, u"Output parameters" ), wx.VERTICAL )
		
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
		
		self.query_filter_less_check = wx.CheckBox( self.output_ctrl_panel, wx.ID_ANY, u"<", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer112.Add( self.query_filter_less_check, 0, 0, 5 )
		
		self.query_filter_less_text = wx.TextCtrl( self.output_ctrl_panel, wx.ID_ANY, u"0.0", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.query_filter_less_text.SetMaxSize( wx.Size( 70,-1 ) )
		
		bSizer112.Add( self.query_filter_less_text, 0, 0, 5 )
		
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
		self.output_to_box.AddPage( self.m_panel29, u"Output Form", False )
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
		
		self.output_sep_custom_radio = wx.RadioButton( self.m_panel30, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer4711.Add( self.output_sep_custom_radio, 0, wx.ALL|wx.LEFT|wx.RIGHT, 5 )
		
		self.output_sep_custom_text = wx.TextCtrl( self.m_panel30, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 35,-1 ), 0 )
		bSizer4711.Add( self.output_sep_custom_text, 0, wx.ALL|wx.RIGHT, 5 )
		
		bSizer4611.Add( bSizer4711, 1, wx.EXPAND, 5 )
		
		gbSizer91.Add( bSizer4611, wx.GBPosition( 0, 1 ), wx.GBSpan( 1, 1 ), wx.EXPAND, 5 )
		
		sbSizer2911.Add( gbSizer91, 1, wx.EXPAND, 5 )
		
		bSizer321.Add( sbSizer2911, 1, wx.EXPAND, 5 )
		
		bSizer94.Add( bSizer321, 1, wx.EXPAND, 5 )
		
		self.m_panel30.SetSizer( bSizer94 )
		self.m_panel30.Layout()
		bSizer94.Fit( self.m_panel30 )
		self.output_to_box.AddPage( self.m_panel30, u"File", True )
		sbSizer52.Add( self.output_to_box, 0, wx.EXPAND |wx.ALL, 5 )
		
		sbSizer511.Add( sbSizer52, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer93.Add( sbSizer511, 0, wx.ALL|wx.EXPAND, 5 )
		
		self.output_ctrl_panel.SetSizer( bSizer93 )
		self.output_ctrl_panel.Layout()
		bSizer93.Fit( self.output_ctrl_panel )
		self.fastSemSim_listbook.AddPage( self.output_ctrl_panel, u"Output Settings", False )
		self.controls_panel = wx.Panel( self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		bSizer62 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer35 = wx.StaticBoxSizer( wx.StaticBox( self.controls_panel, wx.ID_ANY, u"Log" ), wx.VERTICAL )
		
		self.controls_log_text = wx.TextCtrl( self.controls_panel, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( -1,200 ), wx.TE_MULTILINE|wx.TE_READONLY )
		sbSizer35.Add( self.controls_log_text, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer62.Add( sbSizer35, 0, wx.ALL|wx.EXPAND, 5 )
		
		bSizer63 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.controls_start_button = wx.Button( self.controls_panel, wx.ID_ANY, u"Start", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer63.Add( self.controls_start_button, 0, wx.ALL, 5 )
		
		self.controls_stop_button = wx.Button( self.controls_panel, wx.ID_ANY, u"Stop", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer63.Add( self.controls_stop_button, 0, wx.ALL, 5 )
		
		bSizer62.Add( bSizer63, 0, wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		sbSizer321 = wx.StaticBoxSizer( wx.StaticBox( self.controls_panel, wx.ID_ANY, u"Tips" ), wx.HORIZONTAL )
		
		self.controls_tip_label = wx.StaticText( self.controls_panel, wx.ID_ANY, u"FastSemSimGui version 2. \nYou should first load a Gene Ontology and then an Annotation Corpus (using the first two panels).\nSelect a Semantic Similarity measure using the Semantic Similarity panel. Enter a query using the Query panel, and", wx.DefaultPosition, wx.Size( 450,-1 ), 0 )
		self.controls_tip_label.Wrap( -1 )
		sbSizer321.Add( self.controls_tip_label, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		bSizer62.Add( sbSizer321, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		self.controls_panel.SetSizer( bSizer62 )
		self.controls_panel.Layout()
		bSizer62.Fit( self.controls_panel )
		self.fastSemSim_listbook.AddPage( self.controls_panel, u"Controls", True )
		self.output_panel = wx.Panel( self.fastSemSim_listbook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.fastSemSim_listbook.AddPage( self.output_panel, u"Output", False )
		
		fastSemSim_sizer_1.Add( self.fastSemSim_listbook, 1, wx.EXPAND |wx.ALL, 5 )
		
		self.SetSizer( fastSemSim_sizer_1 )
		self.Layout()
		self.fastSemSim_statusbar = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
		
		self.Centre( wx.BOTH )
	
	def __del__( self ):
		pass
	

###########################################################################
## Class GO_load_gui
###########################################################################

class GO_load_gui ( wx.Dialog ):
	
	def __init__( self, parent ):
		wx.Dialog.__init__ ( self, parent, id = wx.ID_ANY, title = u"Load Gene Ontology", pos = wx.DefaultPosition, size = wx.Size( 650,410 ), style = wx.DEFAULT_DIALOG_STYLE )
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		bSizer3 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer22 = wx.BoxSizer( wx.HORIZONTAL )
		
		sbSizer3 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Source file" ), wx.HORIZONTAL )
		
		self.GO_load_source_label = wx.StaticText( self, wx.ID_ANY, u"No file selected. The built-in version will be loaded.", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_RIGHT|wx.ST_NO_AUTORESIZE )
		self.GO_load_source_label.Wrap( -1 )
		sbSizer3.Add( self.GO_load_source_label, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL|wx.EXPAND, 5 )
		
		bSizer22.Add( sbSizer3, 1, wx.LEFT|wx.RIGHT, 5 )
		
		bSizer20 = wx.BoxSizer( wx.HORIZONTAL )
		
		self.GO_load_select_button = wx.Button( self, wx.ID_ANY, u"Select file...", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer20.Add( self.GO_load_select_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		self.GO_load_default_button = wx.Button( self, wx.ID_ANY, u"Set default", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer20.Add( self.GO_load_default_button, 0, wx.ALIGN_CENTER, 5 )
		
		bSizer22.Add( bSizer20, 0, wx.ALIGN_RIGHT|wx.RIGHT|wx.TOP, 5 )
		
		bSizer3.Add( bSizer22, 0, wx.ALIGN_CENTER|wx.EXPAND|wx.TOP, 5 )
		
		bSizer5 = wx.BoxSizer( wx.VERTICAL )
		
		bSizer23 = wx.BoxSizer( wx.VERTICAL )
		
		sbSizer5 = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Options" ), wx.HORIZONTAL )
		
		fgSizer1 = wx.FlexGridSizer( 2, 2, 0, 0 )
		fgSizer1.SetFlexibleDirection( wx.BOTH )
		fgSizer1.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.GO_load_ignore_haspart_check = wx.CheckBox( self, wx.ID_ANY, u"Ignore has_part", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer1.Add( self.GO_load_ignore_haspart_check, 1, wx.ALL, 5 )
		
		self.m_staticText59 = wx.StaticText( self, wx.ID_ANY, u"Check this to ignore 'has_part' relationships", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText59.Wrap( -1 )
		fgSizer1.Add( self.m_staticText59, 1, wx.ALIGN_RIGHT|wx.ALL, 5 )
		
		self.GO_load_ignore_regulates_button = wx.CheckBox( self, wx.ID_ANY, u"Ignore regulates", wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer1.Add( self.GO_load_ignore_regulates_button, 1, wx.ALL, 5 )
		
		self.m_staticText591 = wx.StaticText( self, wx.ID_ANY, u"Check this to ignore 'regulates' relationships", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText591.Wrap( -1 )
		fgSizer1.Add( self.m_staticText591, 1, wx.ALIGN_RIGHT|wx.ALL, 5 )
		
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
		
		self.m_staticText31 = wx.StaticText( self, wx.ID_ANY, u"The Gene Ontology file must be encoded in obo-xml format. It can also be gzipped. Please refer to www.geneontology.org to retrieve the most updated Gene Ontology. If not provided, the Gene Ontology provided with fastSemSim will be loaded.", wx.DefaultPosition, wx.Size( 500,-1 ), 0 )
		self.m_staticText31.Wrap( -1 )
		self.m_staticText31.SetMaxSize( wx.Size( 500,-1 ) )
		
		sbSizer18.Add( self.m_staticText31, 0, wx.ALL, 5 )
		
		bSizer3.Add( sbSizer18, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )
		
		self.SetSizer( bSizer3 )
		self.Layout()
		
		self.Centre( wx.BOTH )
	
	def __del__( self ):
		pass
	

###########################################################################
## Class AC_load_gui
###########################################################################

class AC_load_gui ( wx.Dialog ):
	
	def __init__( self, parent ):
		wx.Dialog.__init__ ( self, parent, id = wx.ID_ANY, title = u"Load Annotation Corpus", pos = wx.DefaultPosition, size = wx.Size( 685,442 ), style = wx.DEFAULT_DIALOG_STYLE )
		
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
		
		self.AC_load_tax_label = wx.TextCtrl( self.m_panel9, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		sbSizer23.Add( self.AC_load_tax_label, 0, wx.ALL, 5 )
		
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
		self.AC_load_type_box.AddPage( self.m_panel9, u"GAF-2", True )
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
		bSizer46.Add( self.AC_load_sep_tab_radio, 0, wx.LEFT|wx.RIGHT|wx.TOP, 5 )
		
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
		
		self.AC_loiad_multi_check = wx.CheckBox( self.m_panel10, wx.ID_ANY, u"Multiple associations", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer48.Add( self.AC_loiad_multi_check, 0, wx.LEFT|wx.RIGHT, 5 )
		
		fgSizer8.Add( bSizer48, 1, wx.EXPAND, 5 )
		
		self.m_staticText50 = wx.StaticText( self.m_panel10, wx.ID_ANY, u"Check this option if each row contains more than one association between objects and GO Terms.", wx.DefaultPosition, wx.Size( 325,-1 ), 0 )
		self.m_staticText50.Wrap( -1 )
		fgSizer8.Add( self.m_staticText50, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
		
		self.m_staticText511 = wx.StaticText( self.m_panel10, wx.ID_ANY, u"Data order", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText511.Wrap( -1 )
		fgSizer8.Add( self.m_staticText511, 0, wx.ALL, 5 )
		
		bSizer49 = wx.BoxSizer( wx.VERTICAL )
		
		self.AC_load_obj_first_radio = wx.RadioButton( self.m_panel10, wx.ID_ANY, u"Object -> GO Term(s)", wx.DefaultPosition, wx.DefaultSize, 0 )
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
		self.AC_load_type_box.AddPage( self.m_panel10, u"Plain File", False )
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
	
	def __del__( self ):
		pass
	

