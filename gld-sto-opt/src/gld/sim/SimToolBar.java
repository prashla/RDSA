
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

import gld.GLDToolBar;
import gld.utils.*;

import java.awt.*;
import java.awt.event.*;


/**
 *
 * The ToolBar for the simulator
 *
 * @author Group GUI
 * @version 1.0
 */

public class SimToolBar extends GLDToolBar
{
	protected static final int STEP = APPBUTTON;
	protected static final int RUN = APPBUTTON + 1;
	protected static final int PAUSE = APPBUTTON + 2;
	protected static final int STOP = APPBUTTON + 3;
	
	protected Choice speed;

	public SimToolBar(SimController sc) {
		super(sc, false);

		addSeparator();
	}
	
	protected void addTools() {
  	addButton("gld/images/edgenode.gif", this, EDGENODE);
  	
  	addSeparator();
  	
  	addButton("gld/images/step.gif", this, STEP);
  	addButton("gld/images/run.gif", this, RUN);
  	addButton("gld/images/pause.gif", this, PAUSE);
  	addButton("gld/images/stop.gif", this, STOP);
  	
  	addSeparator();
		
		speed = new Choice();
		for (int i=0; i < SimController.speedTexts.length; i++)
			speed.add(SimController.speedTexts[i]);
		speed.addItemListener(this);
		speed.setSize(100, 20);
		addComponent(speed);
	}
	
	public Choice getSpeed() { return speed; }

	public void actionPerformed(ActionEvent e) {
		super.actionPerformed(e);
		SimController sc = (SimController)controller;
		int Id = ((IconButton)e.getSource()).getId();
		switch (Id) {
			case STEP  : { sc.doStep(); break; }
			case RUN   : { sc.unpause(); break; }
			case PAUSE : { sc.pause(); break; }
			case STOP  : { sc.stop(); break; }
		}
	}
	
	public void itemStateChanged(ItemEvent e) {
		SimController sc = (SimController)controller;
		Choice s = (Choice)e.getItemSelectable();
		if (s == speed) sc.setSpeed(speed.getSelectedIndex());
		else super.itemStateChanged(e);
	}
}