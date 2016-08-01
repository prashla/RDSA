
/*-----------------------------------------------------------------------
 * Copyright (C) 2001 Green Light District Team, Utrecht University 
 *
 * This program (Green Light District) is free software.
 * You may redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by
 * the Free Software Foundation (version 2 or later).
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * See the documentation of Green Light District for further information.
 *------------------------------------------------------------------------*/

package gld;

import gld.edit.*;
import gld.sim.*;
import gld.xml.*;
import gld.tools.*;
import gld.infra.*;
import gld.utils.*;
import gld.config.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/**
 *
 * The main controller of the application . It starts the user interface and 
 * takes care of interaction with the user.
 *
 * @author Group GUI
 * @version 1.0
 */

public abstract class Controller extends Frame implements XMLSerializable
{
	/** Provides the online help system. */
	protected static HelpViewer helpViewer = null;
	/** The controller that started the app */
	protected static Controller mainController = null;

	/** The Model (needed to change its observers) */
	protected Model model;
	/** The MainView showing the model, controlled by this Controller */
	protected MainView view;
	/** ViewScroller that holds the view */
	protected ViewScroller viewScroller;
	/** The toolbar of this controller */
	protected GLDToolBar toolBar;
	/* The panel containing tool-specific control components. */
	protected Panel toolPanel;
	/** The statusbar (positioned at the bottom of the frame) */
	protected Label statusBar;
	/** The current settings */
	protected Settings settings;

	/** The currently selected tool */
	protected Tool currentTool;
	/** The name of the currently opened file, as used by 'Save' */
	protected String currentFilename;
	/** The current selection on the main view */
	protected Selection currentSelection;
	
	/** The configuration dialog */
	protected ConfigDialog configDialog;


	/**
	 * Constructs a new Controller. Initializes the whole frame and it components.
	 *
	 * @param m The <code>Model</code> to control.
	 */	
	public Controller(Model m, boolean splash)
	{
		SplashScreen ss = null;
		if (splash) {
			ss = new SplashScreen(this);
			ss.setBounds(300, 300, 396, 225);
			ss.setVisible(true);
		}
		
		model = m;

		setLayout(null);
		setBackground(Color.lightGray);
		addWindowListener(new WindowAdapter() { public void windowClosing(WindowEvent e) { quit(); } });
		addComponentListener(new ComponentAdapter() { public void componentResized(ComponentEvent e) { resizeComponents(); } });
		setMenuBar(createMenuBar());

		// creating our view
		view = new MainView(m.getInfrastructure());
		view.setViewPosition(new Point(-2500, -1500));
		view.setViewportPosition(new Point(-2500, -1500));
		viewScroller = new ViewScroller(view);
		add(viewScroller);

		// Create main selection
		Class[] filter = { Node.class, Drivelane.class };
		currentSelection = new Selection(view, filter, model.getInfrastructure());
		currentSelection.addObserver(view);
		view.addOverlay(currentSelection);

		// Create tool listener
		ToolListener toolListener = new ToolListener(this, view);
		view.addMouseMotionListener(toolListener);
    	view.addMouseListener(toolListener);

		// Create toolbar
		toolBar = createToolBar();
		add(toolBar);
		changeTool(new SelectTool(this));

		// Create status bar		
		statusBar = new Label();
		statusBar.setBackground(Color.lightGray);
		add(statusBar);
		
		
		// Misc
		setViewEnabled(true);
		setStatus("Ready.");
		setCurrentFilename(null);
		if (mainController == null) 
		{	mainController = this;
			settings=new GeneralSettings();
			try
			{	settings.doLoad();
			}
			catch (Exception e)
			{	showError ("Cannot load settings. Reverting to defaults. Cause : "+e);
			}
		}
		else
		{	settings=GeneralSettings.getCurrentSettings();
		}
		
		configDialog = new ConfigDialog(this);
		setSize(600, 400);
		
		
		// Give VM time to setup the window
		try { Thread.currentThread().sleep(splash ? 3000 : 1000); } catch (Exception e) { }

		// show the window
		setVisible(true);
		
		// Focus on the splash screen while VM lays out the main window
		if (splash) ss.requestFocus();

		viewScroller.center();
		
		if (splash) {
			ss.setVisible(false);
			ss.dispose();
		}
	}







	/*============================================*/
	/* GET and SET methods                        */
	/*============================================*/


	/** Returns the config dialog */
	public ConfigDialog getConfigDialog() { return configDialog; }
	/** Sets the config dialog */
	public void setConfigDialog(ConfigDialog cd) { configDialog = cd; }
	
	/** Returns the main controller */
	public static Controller getMainController() { return mainController; }
	/** Sets the main controller */
	public static void setMainController(Controller c) { mainController = c; }

  /** Returns the current tool panel */
  public Panel getToolPanel() { return toolPanel; }
  /** Sets the current tool panel */
  public void setToolPanel(Panel p) { toolPanel = p; }

	/** Returns the current selection */
	public Selection getCurrentSelection() { return currentSelection; }
	/** Sets thec current selection */
	public void setCurrentSelection(Selection s) { currentSelection = s; }

	/** Returns the current tool */
	public Tool getCurrentTool() { return currentTool; }
	
	/** Returns the current model */
	public Model getModel() { return model; }

  /** Sets the current tool */
  public void setCurrentTool(Tool t) { currentTool = t; }

	/** Returns the current <code>View</code> */
	public View getView() { return view; }
	
	/** Returns the container of the view */
	public ViewScroller getViewScroller() { return viewScroller; }

	/** Returns the text on the status bar */
	public String getStatus() { return statusBar.getText(); }
	/** Sets the text on the status bar */
	public void setStatus(String s) { statusBar.setText(s); }
	
	/** Sets visibility of View */
	public void setViewEnabled(boolean b) 
	{ 
		if(b) {
			setStatus("View enabled.");
			model.addObserver(view);
		} else {
			setStatus("View disabled.");
			model.deleteObserver(view);
		}
		view.setVisible(b); 
	}
	/** Returns true if the view is currently enabled */
	public boolean isViewEnabled() { return view.isVisible(); }

	/** Sets name to be the currently opened file. */
	public void setCurrentFilename(String name)
	{
		currentFilename = name;
		if (currentFilename == null)
			setTitle("Green Light District " + appName());
		else
			setTitle("Green Light District " + appName() + " - " + currentFilename);
	}










	/*============================================*/
	/* Some general methods                       */
	/*============================================*/

	/** User wants to quit */
	public void quit()
	{
		if (this == mainController)
		{
			try {
				settings.doSave();
			}
			catch (Exception e) {
				showError("Cannot save settings : "+e);
			}
			// TODO: ask if user wants to save work
			System.exit(0);
		}
		else {
			setVisible(false);
			dispose();
			configDialog.setVisible(false);
			configDialog.dispose();
		}
	}

	/** Resizes all <code>Components</code> in the frame to make them all fit. */
  protected void resizeComponents() 
  {
  	int width = this.getWidth();
  	int height = this.getHeight();
  	
  	toolBar.setBounds(7, 50, 800, 24);
  	viewScroller.setBounds(5, 80, width - 13, height - 113);
  	statusBar.setBounds(5, height - 27, width - 5, 25);
  }
  
	/**
	 * Shows an error dialog.
	 * @param msg The message to be shown.
	 */
	public void showError(String msg) {
		new ErrorDialog(this, msg);
	}
	
	/**
	 * Opens a help HTML-page.
	 *
	 * @param helpItem The specified item to be opened. (One of defined 'HelpViewer.HELP_xxx' constants.)
	 */
  public void showHelp(int helpItem)
	{
		Settings settings = GeneralSettings.getCurrentSettings();
		if(settings.getPropertyBooleanValue("use-jbrowser")) {
			if(helpViewer == null)
				helpViewer = new HelpViewer(this);
			helpViewer.showHelp(helpItem);
		} else {
			String cmd = settings.getPropertyStringValue("browser-command") + " \"";
			cmd += HelpViewer.getHelpItem(helpItem) + "\"";
			try { Runtime.getRuntime().exec(cmd); }
			catch(IOException e) { showError("Couldn't launch browser: "+cmd); }
		}
	}

	/** Shows the config dialog if it is hidden, hides it, if it is currently showing */
	public void switchConfigDialog() {
		configDialog.setVisible(!configDialog.isVisible());
	}
	/** Shows the config dialog */
	public void showConfigDialog() {
		configDialog.setVisible(true);
	}
	/** Shows the config dialog with the general panel. */
	public void showSettings() {
		configDialog.showGeneralPanel();
		configDialog.setVisible(true);
	}
	
	
	/** Reports unsuspected errors */
	public static void reportError(Throwable e) {
		System.out.println("Internal error:");
		System.out.println(e.getMessage());
		e.printStackTrace();
//		mainController.showError("An internal error ocurred");
	}
	
	/** Zooms the main view to the specified zoom factor index */
	public void zoomTo(int index) {
		viewScroller.zoomTo(index);
		zoomChanged();
	}
	
	/** Zooms the main view one step in */
	public void zoomIn(Point p) {
		viewScroller.zoomIn(p);
		zoomChanged();
	}
	
	/** Zooms the main view one step out 
	* @param p The Point to be centered at after zoom out
	*/
	public void zoomOut(Point p) {
		viewScroller.zoomOut(p);
		zoomChanged();
	}
	
	protected void zoomChanged() {
		toolBar.getZoom().select(view.getZoomIndex());
	}
	
	/** Runs the validator on the current infrastructure */
	public boolean validateInfra() {
		try {
			Validation v = new Validation(getModel().getInfrastructure());
			Vector errors = v.validate();
			if (!errors.isEmpty()) {
				String errorString ="";
				for (Enumeration e = errors.elements(); e.hasMoreElements();)
				{
					errorString+=e.nextElement()+"\n";
				}
				showError("Validation error(s):\n\n" + errorString);
				return false;
			}
		}
		catch (GLDException e) {
			reportError(e.fillInStackTrace());
		}
		setStatus("Infrastructure is valid.");
		view.redraw();
		return true;
	}
	
	/** Set a new infrastructure */
	public void newInfrastructure(Infrastructure infra) {
		//infra.reset();
		CustomFactory.reset();
		viewScroller.center();
		model.setInfrastructure(infra);
		view.setInfrastructure(infra);
		currentSelection.setSelectionStarter(infra);
		currentSelection.deselectAll();
		viewScroller.resizeInfra(infra.getSize());
		viewScroller.center();
		changeTool(new SelectTool(this));
	}
	
	/** Changes the current tool */
	public void changeTool(Tool t)
	{
  	if (toolPanel != null) toolBar.remComponent(toolPanel);
  	toolPanel = t.getPanel();
  	toolBar.addComponent(toolPanel);
		toolBar.doLayout();
		toolPanel.doLayout();
		toolBar.repaint();
		view.remOverlay(currentTool);
		view.addOverlay(t);
		currentTool = t;
	}





	
















	/*============================================*/
	/* Loading and saving                         */
	/*============================================*/

	public void newFile()
	{
		EditSizeDialog esd = new EditSizeDialog(this, 5000, 3000);
		esd.show();
				
		if (esd.ok()) {
			Dimension d = new Dimension(esd.getWidthI(), esd.getHeightI());
			setCurrentFilename(null);
			newInfrastructure(new Infrastructure(d));
			setStatus("New infrastructure.");
		}
	}

	/** User wants to load a file */
	public void openFile() {
		FileDialog diag = new FileDialog(this, "Open file...", FileDialog.LOAD);
		diag.show();
		String filename = diag.getFile();
		if (filename == null) return;
		filename = diag.getDirectory() + filename;
		tryLoad(filename);
	}
	
	/** User wants to save file as... */
	public void saveFileAs()
	{
		String filename;
		if ((filename = showSaveAsDialog()) == null) return;
		trySave(filename);
	}
	
	/** User wants to save */
	public void saveFile() {
		if (currentFilename == null) saveFileAs();
		else trySave(currentFilename);
	}
	
	/** Try to load a file, process exceptions
	  * @param filename The name of the file to load from
	 */
	protected void tryLoad(String filename)
	{	try {
			setStatus("Loading file " + filename);
			doLoad(filename);
			setCurrentFilename(filename);
			setStatus("Succesfully loaded "+filename);
		}
		catch (InvalidFilenameException e) {
			setStatus ("Cannot load "+filename);
			showError(e.getMessage());
		}
		catch (Exception e) {
			setStatus ("Cannot load "+filename);
			e.printStackTrace(System.out);
		}
		

	}
	
	/** Try to save file  , process exceptions 
	  * @param filename The name of the file to save to
	 */
	protected void trySave(String filename) {
		try {
			model.getInfrastructure().cachInboundLanes();
			doSave(filename);
			setCurrentFilename(filename);
		}
		catch (InvalidFilenameException e) {
			showError(e.getMessage());
		}
		catch (Exception e) {
			e.printStackTrace(System.out);
		}
	}

	/** Shows the 'Save as...' dialog and returns the filename chosen */
	protected String showSaveAsDialog()
	{
		FileDialog diag = new FileDialog(this, "Save as...", FileDialog.SAVE);
		diag.show();
		String filename;
		if((filename = diag.getFile()) == null) return null;
		return diag.getDirectory() + filename;
	}
	
	/** Load state from disk. Return true upon success. */
	protected abstract void doLoad(String filename) throws InvalidFilenameException, Exception;
	/** Save the current state to disk. Return true upon success. */
	protected abstract void doSave(String filename) throws InvalidFilenameException, Exception;
	
	/** Loads a model,the current view and the current controller from a XMLLoader 
	  * @param loader The loader to load everything from 
	  * @param model The model to load
	 */
	protected void loadAll (XMLLoader loader,Model model) throws IOException,XMLTreeException,XMLInvalidInputException
	{    	loader.load(null,model);
			loader.load(null,getView());
			loader.load(null,this);
	}
	
	/** Saves a model,the current view and the current controller to a XMLSaver 
	  * @param saver The XMLSaver to save everything to
	  * @param model The model to save
	 */
	protected void saveAll (XMLSaver saver,Model model) throws IOException,XMLTreeException,XMLCannotSaveException
	{	saver.saveObject(model); 
    	saver.saveObject(getView());
    	saver.saveObject(this); 
	}

	/*============================================*/
	/* XMLSerializable                            */
	/*============================================*/

 	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
 	{ 	int width=myElement.getAttribute("width").getIntValue();
  		int height=myElement.getAttribute("height").getIntValue();
  		setSize(width,height);
  		currentSelection.setSelectionStarter(model.getInfrastructure());
 	}

 	public XMLElement saveSelf () throws XMLCannotSaveException
 	{ 	XMLElement result=new XMLElement(getXMLName());
   		result.addAttribute(new XMLAttribute("width",super.getWidth()));
   		result.addAttribute(new XMLAttribute("height",super.getHeight()));
   		return result;
 	}
  
 	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
 	{ 	// A controller has no childs in the XML structure. The model and view
   		// are handled by the parser.
 	}

 	public String getXMLName ()
 	{ return "controller";
 	}
 
	public void setParentName (String parentName) throws XMLTreeException
	{	throw new XMLTreeException
	        ("Cannot set parent for XML root class controller");
	}



















	protected abstract String appName();
	public abstract void showFilePropertiesDialog();
	public abstract PopupMenu getPopupMenuFor(Selectable obj) throws PopupException;
	protected abstract GLDToolBar createToolBar();
	protected abstract MenuBar createMenuBar();
}