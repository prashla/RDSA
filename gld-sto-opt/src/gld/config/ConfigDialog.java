
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

package gld.config;

import gld.*;
import gld.infra.*;
import gld.tools.*;
import gld.utils.*;
import gld.sim.SimController;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

/**
 *
 * @author Group GUI
 * @version 1.0
 */

public class ConfigDialog extends Frame implements Observer, WindowListener, ActionListener
{
	protected static final int WIDTH = 400;
	protected static final int HEIGHT = 350;

	protected static final int TITLE_XPOS = 5;
	protected static final int TITLE_YPOS = 25;
	protected static final int TITLE_WIDTH = WIDTH - TITLE_XPOS * 2;
	protected static final int TITLE_HEIGHT = 25;

	protected static final int CLOSE_WIDTH = 50;
	protected static final int CLOSE_HEIGHT = 24;
	protected static final int CLOSE_XPOS = (int)(WIDTH / 2 - CLOSE_WIDTH / 2) - 1;
	protected static final int CLOSE_YPOS = HEIGHT - CLOSE_HEIGHT - 10;

	protected static final int PANEL_XPOS = 10;
	protected static final int PANEL_YPOS = TITLE_YPOS + TITLE_HEIGHT + 10;
	protected static final int PANEL_WIDTH = WIDTH - PANEL_XPOS * 2;
	protected static final int PANEL_HEIGHT = CLOSE_YPOS - PANEL_YPOS - 10;
	
	
	public static boolean AlwaysOnTop = true;


	Controller controller;
	SubPanel subPanel;
	Button close;
	PanelFactory cpf;
	Label title;

	public ConfigDialog(Controller con)
	{
		super("Configuration dialog");

		controller = con;
		con.getCurrentSelection().addObserver(this);
		con.getModel().addObserver(this);

		cpf = new PanelFactory(this,
			con instanceof SimController ? PanelFactory.TYPE_SIM : PanelFactory.TYPE_EDIT);

		setBounds(100, 100, WIDTH, HEIGHT);
		setLayout(null);
		addWindowListener(this);
		setResizable(false);
		setBackground(SystemColor.control);

		title = new Label("", Label.CENTER);
		title.setFont(new Font(null, 0, 16));
		title.setBounds(TITLE_XPOS, TITLE_YPOS, TITLE_WIDTH, TITLE_HEIGHT);
		add(title);

		subPanel = new SubPanel();
		subPanel.setBounds(PANEL_XPOS, PANEL_YPOS, PANEL_WIDTH, PANEL_HEIGHT);
		add(subPanel);

		close = new Button("OK");
		close.setBounds(CLOSE_XPOS, CLOSE_YPOS, CLOSE_WIDTH, CLOSE_HEIGHT);
		close.addActionListener(this);
		add(close);
		
		subPanel.setConfigPanel(new GeneralPanel(this));
	}
	
	
	public void update(Observable o, Object arg) {
		if (o instanceof Model) {
			ConfigPanel cp = subPanel.getConfigPanel();
			if (cp != null) cp.reset();
			return;
		}
		Selection s = (Selection)o;
		try {
			ConfigPanel cp = subPanel.getConfigPanel();
			cp.ok();
			subPanel.setConfigPanel(cpf.createPanel(s));
			if (AlwaysOnTop) toFront();
		}
		catch (ConfigException e) {
			Controller.reportError(e);
		}
	}

	
	/** Returns the controller that created this config dialog */
	public Controller getController() { return controller; }
	/** Sets the controller that created this config dialog */
	public void setController(Controller con) { controller = con; }
	
	/** Returns the current configuration panel */
	public ConfigPanel getConfigPanel() { return subPanel.getConfigPanel(); }
	/** Sets the current configuration panel */
	public void setConfigPanel(ConfigPanel cp) { subPanel.setConfigPanel(cp); }
	
	/** Returns the title of this dialog */
	public String getTitle() { return title.getText(); }
	/** Sets the title of this dialog */
	public void setTitle(String newtitle) {
		title.setText(newtitle);
		super.setTitle(newtitle);
	}

	/** Shows a message dialog with an OK button */
	public void showError(String msg) {
		new ErrorDialog(this, msg);
	}
	
	/** Changes the current selection to select the given object */
	public void selectObject(Selectable s) {
		controller.getCurrentSelection().newSelection(s);
	}
	
	/** Show the general panel. */
	public void showGeneralPanel() {
		ConfigPanel cp = subPanel.getConfigPanel();
		cp.ok();
		subPanel.setConfigPanel(new GeneralPanel(this));
	}
	
	/** Shows the road user panel. */
	public void showRoaduser(Roaduser ru) {
		ConfigPanel cp = subPanel.getConfigPanel();
		cp.ok();
		subPanel.setConfigPanel(new RoaduserPanel(this, ru));
	}





	public void setVisible(boolean b) {
		if (b == false) subPanel.getConfigPanel().ok();
		super.setVisible(b);
	}

	public void actionPerformed(ActionEvent e) { setVisible(false); }
	public void windowClosing(WindowEvent e) { hide(); }
	public void windowActivated(WindowEvent e) { }
	public void windowClosed(WindowEvent e) { }
	public void windowDeactivated(WindowEvent e) { }
	public void windowIconified(WindowEvent e) { }
	public void windowDeiconified(WindowEvent e) { }
	public void windowOpened(WindowEvent e) { }

	
	private class SubPanel extends Panel
	{
		ConfigPanel current;
		public SubPanel() {
			setLayout(null);
		}
		public SubPanel(ConfigPanel cp) {
			this();
			setConfigPanel(cp);
		} 
		public void setConfigPanel(ConfigPanel cp)
		{
			if (cp != current && cp != null) {
				if (current != null) remove(current);
				add(cp);
				cp.setBounds(0, 0, getWidth(), getHeight());
				current = cp;
				doLayout();
			}
		}
		public ConfigPanel getConfigPanel() {
			return current;
		}
	}
}