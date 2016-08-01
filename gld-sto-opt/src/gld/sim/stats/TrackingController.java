
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

package gld.sim.stats;

import gld.ErrorDialog;
import gld.sim.SimModel;
import gld.sim.SimController;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
*
* The controller for the tracking window, it controls a TrackingView.
*
* @author Group GUI
* @version 1.0
*/

public class TrackingController extends Frame
{
	private SimModel model;
	private SimController controller;
	private TrackingView view;

	/**
	* Creates a <code>TrackingController</code>.
	*
	* @param _model The <code>SimModel</code> statistics should be read from.
	* @param _controller The parent <code>SimController</code>.
	* @param _view The <code>TrackingView</code> to be shown.
	*/
	public TrackingController(SimModel _model, SimController _controller, TrackingView _view)
	{
		model = _model;
		controller = _controller;
		view = _view;
		
		addWindowListener(new WindowAdapter() { public void windowClosing(WindowEvent e) { closeWindow(); } });
		setBounds(200, 200, 400, 300);
		setBackground(Color.lightGray);
		setTitle("Tracking - " + view.getDescription());
		setVisible(true);

 		add(view);
 		model.addObserver(view);
		setViewEnabled(true);

		setMenuBar(makeMenuBar());
	}
	
	/**
	* Closes the <code>TrackingController</code>.
	*/
	public void closeWindow()
	{
		setVisible(false);
		model.deleteObserver(view);
		dispose();
	}
	
	/** Enables or disables the view. */
	public void setViewEnabled(boolean enable) { view.setVisible(enable); view.redraw(); }
	/** Returns true if the view is enabled. */
	public boolean isViewEnabled() { return view.isVisible(); }

	/**
	* Creates the <code>MenuBar</code> to be used.
	*/
	public MenuBar makeMenuBar()
	{
		MenuBar bar = new MenuBar();
		Menu menu; MenuItem item;
		
		menu = new Menu("File");
  	bar.add(menu);
  	MenuListener ml = new MenuListener(this);

  	item = new MenuItem("Export...", new MenuShortcut(KeyEvent.VK_S));
  	menu.add(item);
  	item.addActionListener(ml);
  	
  	menu.add(new MenuItem("-"));

  	item = new MenuItem("Close", new MenuShortcut(KeyEvent.VK_W));
  	menu.add(item);
  	item.addActionListener(ml);

		/*  Options */

    menu = new Menu("Options");
  	bar.add(menu);

  	CheckboxMenuItem citem = new CheckboxMenuItem("Toggle view", true);
  	menu.add(citem);
   	citem.addItemListener(ml);

		addToOptionsMenu(menu);

 		return bar;
 	}
 	
 	protected void addToOptionsMenu(Menu menu) {}

	/**
	 * Shows an error dialog.
	 * @param msg The message to be shown.
	 */
	public void showError(String msg) {
		new ErrorDialog(this, msg);
	}

	/** Resets the view. */
	public void reset() {
		view.reset();
	}

	/**
	* Listens to the menus.
	*/
	protected class MenuListener implements ActionListener, ItemListener
	{
		TrackingController controller;
		
		public MenuListener(TrackingController _controller)
		{
			controller = _controller;
		}
		
		/**
		* Handles the <code>ActionEvent</code> action.
		* @param e The ActionEvent that has occured.
		*/
		public void actionPerformed( ActionEvent e ) 
		{
			String sel = ((MenuItem) e.getSource()).getLabel();
			if (sel.equals("Close")) { closeWindow(); }
			else if(sel.equals("Export...")) 
			{ 
				FileDialog diag = new FileDialog(controller, "Export...", FileDialog.SAVE);
				diag.setFile("export - " + view.getDescription() + ".dat");
				diag.show();
				String filename;
				if((filename = diag.getFile()) == null) return;
				filename = diag.getDirectory() + filename;
				try { view.saveData(filename, model); }
				catch(IOException exc) { controller.showError("Couldn't export data to \"" + filename + "\"!"); }
			}
		}
		
		public void itemStateChanged(ItemEvent e) {
			CheckboxMenuItem item = (CheckboxMenuItem)e.getItemSelectable();
			setViewEnabled(item.getState());
		}
	}
	
	public TrackingView getTrackingView() {
		return view;
	}
}