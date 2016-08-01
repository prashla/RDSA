
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

package gld.sim;

import gld.*;
import gld.infra.*;
import gld.edit.*;
import gld.sim.stats.*;
import gld.algo.dp.*;
import gld.algo.tlc.*;
import gld.tools.*;
import gld.xml.*;
import gld.utils.Arrayutils;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/**
 * The main controller for the simulation part of the application.
 *
 * @author Group GUI
 * @version 1.0
 */

public class SimController extends Controller implements Observer
{
	protected EditController editController = null;
	protected SimMenuBar menuBar;

	public static final String[] speedTexts = { "Low", "Medium", "High", "Maximum" };
	public static final int[] speedSettings = { 1000, 400, 50, 10 };

	protected Choice speedChoice;
	protected StatisticsOverlay statsOverlay;
	
	/**
	* Creates the main frame.
	*
	* @param m The <code>SimModel</code> to be controlled.
	*/
	public SimController(SimModel m, boolean splash) {
		super(m, splash);
		setSimModel(m);
		m.setSimController(this);
		
		speedChoice=new Choice();
		Enumeration e=Arrayutils.getEnumeration(speedTexts);
		while (e.hasMoreElements())
			speedChoice.add((String)(e.nextElement()));
		
		setSpeed((int)(speedTexts.length / 2));
		setCycleCounterEnabled(true);
		
		statsOverlay = new StatisticsOverlay(view,m.getInfrastructure());
		
		setTLC(0, 0);
		setDrivingPolicy(0);
	}







	/*============================================*/
	/* GET and SET methods                        */
	/*============================================*/

	/** Returns the current <code>SimModel</code> */
	public SimModel getSimModel() { return (SimModel)model; }
 
	/** Sets a new <code>SimModel</code> to be controlled */
	public void setSimModel(SimModel m) { model = m; }

	/** Enables or disables the cycle counter. */
	public void setCycleCounterEnabled(boolean b)
	{	if(b) 
			getSimModel().addObserver(this);
		else {
			setStatus("Cycle counter disabled at cycle " + getSimModel().getCurCycle() + ".");
			getSimModel().deleteObserver(this);
		}
	}















	/*============================================*/
	/* Load and save                              */
	/*============================================*/


	public void load(XMLElement myElement,XMLLoader loader)
		throws XMLTreeException, IOException, XMLInvalidInputException
	{	super.load(myElement,loader);
		// TODO restore menu options/choices in GUI
		statsOverlay = new StatisticsOverlay(view,getSimModel().getInfrastructure());
		if (XMLUtils.getLastName(statsOverlay).equals(loader.getNextElementName()))
		{	System.out.println("Loading stats");
			loader.load(this, statsOverlay);
		}
	}

	public XMLElement saveSelf() throws XMLCannotSaveException
	{	XMLElement result = super.saveSelf();
		/* This code is buggy
		result.addAttribute(new XMLAttribute("saved-by", "simulator"));
	 	result.addAttribute(new XMLAttribute("tlc-category",
				menuBar.getTLCMenu().getCategory()));
		result.addAttribute(new XMLAttribute("tlc-number",
				menuBar.getTLCMenu().getTLC()));
		result.addAttribute(new XMLAttribute("driving-policy",
				menuBar.getDPMenu().getSelectedIndex()));
		result.addAttribute(new XMLAttribute("speed",
				speedChoice.getSelectedIndex()));
		*/		
		return result;
	}
	
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{ 	saver.saveObject(statsOverlay);
	}

	public void doSave(String filename) throws InvalidFilenameException, Exception
	{	if(!filename.endsWith(".sim") )
			throw new InvalidFilenameException("Filename must have .sim extension.");
		setStatus("Saving simulation to " + filename);
		XMLSaver saver=new XMLSaver(new File(filename));
		saveAll(saver,getSimModel());
		saver.close();
		setStatus("Saved simulation to " + filename);
	}
	
	public void doLoad(String filename) throws InvalidFilenameException, Exception
	{	if(!filename.endsWith(".infra") && !filename.endsWith(".sim"))
			throw new InvalidFilenameException("You can only load .infra and .sim files.");
		stop();			
		TrackerFactory.purgeTrackers();
		XMLLoader loader=new XMLLoader(new File(filename));
		loadAll(loader,getSimModel());
		newInfrastructure(model.getInfrastructure());
		loader.close();
	}












	/*============================================*/
	/* Miscellanous                               */
	/*============================================*/

	/** Called by observable SimModel (if view enabled). */
	public void update(Observable o, Object arg) 
	{	int cycle=((SimModel)o).getCurCycle();
		if (cycle!=0)
			setStatus("Cycle: " + cycle);
	}

	/** Returns the name of this controller extension. */
	protected String appName() { return "simulator"; }
	
	protected MenuBar createMenuBar() {
		menuBar = new SimMenuBar(this, speedTexts);
		return menuBar;
	}
	
	protected GLDToolBar createToolBar() {
		return new SimToolBar(this);
	}










	/*============================================*/
	/* Invoked by Listeners                       */
	/*============================================*/

	/**
	* Opens the statistics viewer.
	*/
	public void showStatistics()
	{	new StatisticsController(getSimModel(), this);
	}

	/**
	* Shows the tracking window.
	*/	
	public void showTracker(int type)
	{	try	{	TrackerFactory.showTracker(getSimModel(), this, type);	}
		catch(GLDException e) {	reportError(e.fillInStackTrace());		}
	}

	/** Enables the statistics overlay */
	public void enableOverlay() {
		statsOverlay = new StatisticsOverlay(view, getSimModel().getInfrastructure());
    	getSimModel().addObserver(statsOverlay);
    	view.addOverlay(statsOverlay);
	}

	/** Enables the statistics overlay */
	public void disableOverlay() {
		getSimModel().deleteObserver(statsOverlay);
		view.remOverlay(statsOverlay);
	}
 
	public void setDrivingPolicy(int dp)
	{	try	{	getSimModel().setDrivingPolicy((new DPFactory(getSimModel(),getSimModel().getTLController())).getInstance(dp));	}
		catch (Exception e)	{	reportError(e);	}
	}

	public void setTLC(int cat, int nr)
	{	setColearningEnabled(cat == 1);
		try {
			SimModel sm = getSimModel();
			TLCFactory tlcf = new TLCFactory(sm.getInfrastructure(), sm.getRandom());
			TLController tlc = tlcf.genTLC(cat, nr);
			tlc.showSettings(this);
			sm.setTLController(tlc);
			setColearningEnabled((tlc instanceof Colearning));
		}
		catch (GLDException e) {
			reportError(e.fillInStackTrace());
		}
	}

	private void setColearningEnabled(boolean b) {
		if (!b && menuBar.getDPMenu().getSelectedIndex() == DPFactory.COLEARNING) {
			menuBar.getDPMenu().select(DPFactory.SHORTEST_PATH);
			setDrivingPolicy(DPFactory.SHORTEST_PATH);
		}
		((CheckboxMenuItem)menuBar.getDPMenu().getItem(DPFactory.COLEARNING)).setEnabled(b);
	}

	/** Shows the file properties dialog */
	public void showFilePropertiesDialog()
	{
		String simName = getSimModel().getSimName();
		Infrastructure infra = getSimModel().getInfrastructure();
		String comments = infra.getComments();
		
		SimPropDialog propDialog = new SimPropDialog(this, simName, comments);
		
		propDialog.show();
		if(propDialog.ok())	{
			getSimModel().setSimName(propDialog.getSimName());
			infra.setComments(propDialog.getComments());
		}
		this.setStatus("Simulation \"" + getSimModel().getSimName() + "\".");
	}

	/** Creates a right-click popup-menu for the given object */
	public PopupMenu getPopupMenuFor(Selectable obj) throws PopupException {
		SimPopupMenuFactory pmf = new SimPopupMenuFactory(this);
		return pmf.getPopupMenuFor(obj);
	}
	
	/** Returns the filename of the currently loaded file */
	public String getCurrentFilename() {
		return currentFilename;
	}
	
	/** Sets the speed of the simulation */
	public void setSpeed(int speed) {
		((SimToolBar)toolBar).getSpeed().select(speed);
		menuBar.getSpeedMenu().select(speed);
		getSimModel().setSpeed(speedSettings[speed]);
	}
	
	/** Makes model do one step */
	public void doStep() { getSimModel().doStep(); }
	
	/** Paues the simulation */
	public void pause() {
		setStatus("Paused at cycle " + getSimModel().getCurCycle() + ".");
		getSimModel().pause();
	}
	
	/** Resumes or starts the simulation */
	public void unpause() {
		setStatus("Simulation running.");
		getSimModel().unpause();
	}
	
	/** Stops the simulation and resets the infrastructure */
	public void stop() {
		int cycle=getSimModel().getCurCycle() ;
		if (cycle!=0)
			setStatus("Stopped at cycle " + ".");
		try {
			getSimModel().pause();
			getSimModel().reset();
		}
		catch (SimulationRunningException ex) {
			reportError(ex.fillInStackTrace());
			getSimModel().unpause();
		}
	}
	
	/** Starts a series of 10 simulations */
	public void runSeries() {
		setStatus("Running a Series of simulations.");
		TrackerFactory.purgeTrackers();
		try	{
			TrackerFactory.showTracker(getSimModel(), this, TrackerFactory.TOTAL_QUEUE);
			TrackerFactory.showTracker(getSimModel(), this, TrackerFactory.TOTAL_WAIT);
			TrackerFactory.showTracker(getSimModel(), this, TrackerFactory.TOTAL_ROADUSERS);
			TrackerFactory.showTracker(getSimModel(), this, TrackerFactory.TOTAL_JUNCTION);
			TrackerFactory.disableTrackerViews();
		}
		catch(GLDException e) { 
			reportError(e.fillInStackTrace());
		}
		menuBar.setViewEnabled(false);
		menuBar.setCycleCounterEnabled(false);
		this.setViewEnabled(false);
		this.setCycleCounterEnabled(false);
		
		this.setSpeed((speedSettings.length-1));
		
		getSimModel().runSeries();
	}
	
	public void nextSeries() {
		getSimModel().pause();
		int curSeries = getSimModel().getCurSeries();
		// If we have data, save it
		if(curSeries>0) {
			String simName = getSimModel().getSimName();
			String tlcName = TLCFactory.getDescription(TLCFactory.getNumberByXMLTagName(XMLUtils.getLastName(getSimModel().getTLController().getXMLName())));
			String dpName  = DPFactory.getDescription(DPFactory.getNumberByXMLTagName(XMLUtils.getLastName(getSimModel().getDrivingPolicy().getXMLName())));
			TrackingController[] tca = TrackerFactory.getTrackingControllers();
			for(int i=0;i<tca.length;i++) {
				TrackingView tv = tca[i].getTrackingView();
				String filename = simName + " - " + tlcName + " - " + dpName + " - " + tv.getDescription() + " - run " + curSeries + ".dat";
				try { tv.saveData(filename,getSimModel()); }
				catch(IOException exc) { showError("Couldn't save statistical data from series!");}
			}			
		}
		// If we have more runs to run, do so.
		if (curSeries < getSimModel().getNumSeries()) {
			setStatus("Running a series of simulations, currently at: "+curSeries);
			try { getSimModel().reset(); }
			catch(SimulationRunningException e) {}
			getSimModel().nextCurSeries();
			getSimModel().unpause();
		}
		else {
			setStatus("Done running Series of simulations.");
			getSimModel().stopSeries();
			TrackerFactory.purgeTrackers();
		}
	}
	
	/** Opens the editor */
	public void openEditor() {
		if (editController == null) editController = new EditController(new EditModel(), false);
		editController.show();
		editController.requestFocus();
	}

	/** Set temp debug infra */
	protected void setInfra(int nr)
	{	Infrastructure infra;
		switch (nr)
		{	case 1 : infra=new SimpleInfra(); break;
			case 2 : infra=new LessSimpleInfra(); break;
			case 3 : infra=new NetTunnelTest1(); break;
			case 4 : infra=new NetTunnelTest2(); break;
			default: infra=new Infrastructure(); break;
		}
		try
		{	Vector errors=(new Validation(infra)).validate();
			if(!errors.isEmpty())
				showError(errors.toString());
		}
		catch (InfraException e)
		{	reportError(e);
		}
		getSimModel().setCurCycle(0);
		newInfrastructure(infra);
	}
}