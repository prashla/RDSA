
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

import gld.infra.*;
import gld.sim.SimModel;
import gld.sim.SimController;
import gld.utils.CheckMenu;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.IOException;

/**
*
* The controller for the statistics viewer, it controlls the StatisticsView.
*
* @author Group GUI
* @version 1.0
*/

public class StatisticsController extends Frame
{
	/** Currently available views. */
	protected final static String[] viewDescs = { "Summary", "Table" }; //, "Graphical" };
	
	protected SimController parent;
	protected StatisticsView view = null;
	protected StatisticsModel stats;
	
	protected CheckMenu viewCM, modeCM;
	protected Scrollbar sbHorizontal, sbVertical;

	/**
	* Creates a <code>StatisticsController</code>.
	*
	* @param m The <code>SimModel</code> statistics should be read from.
	* @param controller The parent <code>SimController</code>.
	*/
	public StatisticsController(SimModel model, SimController _parent)
	{
		stats = new StatisticsModel(model);
		parent = _parent;
		addWindowListener(new WindowAdapter() { public void windowClosing(WindowEvent e) { closeWindow(); } });
		setBounds(200, 200, 400, 300);
		setBackground(Color.lightGray);

		add(sbVertical = new Scrollbar(Scrollbar.VERTICAL, 0, 1, 0, 1), BorderLayout.EAST);
		sbVertical.addAdjustmentListener(new ScrollListener());
		add(sbHorizontal = new Scrollbar(Scrollbar.HORIZONTAL, 0, 1, 0, 1), BorderLayout.SOUTH);
		sbHorizontal.addAdjustmentListener(new ScrollListener());

		setMenuBar(makeMenuBar());
		setView("Summary");
		refresh();
		setVisible(true);
		view.requestFocus();
	}



	public void setScrollMax(int hor, int ver)
	{
		sbHorizontal.setMaximum(hor);
		sbVertical.setMaximum(ver);
	}




	/*============================================*/
	/* Invoked by listeners                       */
	/*============================================*/
	
	/**
	* Closes the <code>StatisticsController</code>.
	*/
	private void closeWindow()
	{
		setVisible(false);
		dispose();
	}
	
	/** Refreshes the statistical data shown. */
	private void refresh()
	{
		setTitle("Statistics (at cycle " + parent.getSimModel().getCurCycle() + ")");
		stats.refresh();
	}

	/** 
	* Sets the current view mode.	
	* @param view One of the constants in <code>viewDescs[]</code>.
	*/
	protected void setView(String desc)
	{
		modeCM.setEnabled(true);
		if(view != null) {
			remove(view);
			stats.deleteObserver(view);
		}
		
		if(desc.equals("Summary")) {
			view = new StatsSummaryView(this, stats);
			modeCM.setEnabled(false);
		}
		else if(desc.equals("Table")) view = new StatsTableView(this, stats);
		else if(desc.equals("Graphical")) view = new StatsBarView(this, stats);

		add(view);
		stats.addObserver(view);
		doLayout();
		sbHorizontal.setValue(0);
		sbVertical.setValue(0);
		view.requestFocus();
	}




	/**
	 * Exports data to a CSV file. 
	 */
	protected void exportData()
	{
		FileDialog diag = new FileDialog(parent, "Export...", FileDialog.SAVE);
		diag.setFile("export - statistics " + stats.getSimName() + ".dat");
		diag.show();
		String filename;
		if((filename = diag.getFile()) == null) return;
		filename = diag.getDirectory() + filename;
		try { stats.saveData(filename);	}
		catch(IOException exc) { parent.showError("Couldn't export data to \"" + filename + "\"!"); }
	}	






	/*============================================*/
	/* Menubar                                    */
	/*============================================*/

	/**
	* Creates the <code>MenuBar</code> to be used.
	*/
	public MenuBar makeMenuBar()
	{
		MenuBar bar = new MenuBar();
		Menu menu; MenuItem item;
		
		menu = new Menu("File");
  	bar.add(menu);
  	MenuListener ml = new MenuListener();

 		item = new MenuItem("Export...");
 		menu.add(item);
  	item.addActionListener(ml);
  	
  	menu.add(new MenuItem("-"));
  	
  	item = new MenuItem("Close", new MenuShortcut(KeyEvent.VK_W));
  	menu.add(item);
  	item.addActionListener(ml);

 		menu = new Menu("Options");
 		bar.add(menu);

 		item = new MenuItem("Refresh", new MenuShortcut(KeyEvent.VK_R));
 		menu.add(item);
 		item.addActionListener(ml);

		menu.add(new MenuItem("-"));

		viewCM = new CheckMenu("View", viewDescs, false);
		viewCM.addItemListener(ml);
		viewCM.select(0);
		menu.add(viewCM);

		String[] modes = {"all roadusers", "last " + Node.STAT_NUM_DATA + " roadusers"};
		modeCM = new CheckMenu("Show average of", modes, false);
		modeCM.addItemListener(ml);
		modeCM.select(0);
		menu.add(modeCM);

 		return bar;
 	}







	/*============================================*/
	/* Listeners                                  */
	/*============================================*/

	/**
	* Listens to the menus.
	*/
	public class MenuListener implements ActionListener, ItemListener
	{
		/**
		* Handles the <code>ActionEvent</code> action.
		* @param e The ActionEvent that has occured.
		*/
		public void actionPerformed( ActionEvent e ) 
		{
			String sel = ((MenuItem) e.getSource()).getLabel();

			if (sel.equals("Export...")) exportData();
			else if (sel.equals("Close")) closeWindow();
			else if(sel.equals("Refresh")) refresh();
		}
		
		public void itemStateChanged(ItemEvent e)
		{
			CheckMenu menu = (CheckMenu)e.getItemSelectable();
			if(menu == modeCM)
				stats.setAllTimeAvg(menu.getSelectedIndex() == 0);
			else
				setView(menu.getSelectedItem().getLabel());
		}
	}
	
	protected class ScrollListener implements AdjustmentListener
	{
		public void adjustmentValueChanged(AdjustmentEvent e)
		{
			if(e.getSource() == sbHorizontal)
				view.setHorScroll(e.getValue());
			else
				view.setVerScroll(e.getValue());
		}	
	}
}