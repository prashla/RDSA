
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

import java.awt.*;
import java.awt.event.*;

import gld.*;
import gld.infra.*;
import gld.utils.*;

/**
 *
 * @author Group GUI
 * @version 1.0
 */

public class EditJunctionPanel extends ConfigPanel implements ActionListener
{
	Junction junction;
	
	Hyperlink[] roadLinks;
	Label nrsigns, width;
	
	public EditJunctionPanel(ConfigDialog cd, Junction j) {
		super(cd);
		
		String[] dirs = { "north", "east", "south", "west" };
		
		roadLinks = new Hyperlink[4];
		
		for (int i=0; i < 4; i++) {
			Label lab = new Label("Road " + dirs[i] + ": ");
			lab.setBounds(0, i * 20, 100, 20);
			add(lab);
			
			roadLinks[i] = new Hyperlink();
			roadLinks[i].addActionListener(this);
			roadLinks[i].setBounds(100, i * 20, 100, 20);
			add(roadLinks[i]);
		}
		
		nrsigns = new Label();
		nrsigns.setBounds(0, 100, 200, 20);
		add(nrsigns);
		
		width = new Label();
		width.setBounds(0, 120, 200, 20);
		add(width);


		setJunction(j);
	}
	
	public void reset() {
		Road[] roads = junction.getAllRoads();
		
		for (int i=0; i < 4; i++) {
			if (roads[i] != null) {
				roadLinks[i].setText(roads[i].getName());
				roadLinks[i].setEnabled(true);
			}
			else {
				roadLinks[i].setText("null");
				roadLinks[i].setEnabled(false);
			}
		}
		
		nrsigns.setText("Junction has " + junction.getNumRealSigns() + " trafficlights");
		width.setText("Junction is " + junction.getWidth() + " units wide");
	}

	public void setJunction(Junction j) {
		junction = j;
		confd.setTitle(junction.getName());
		reset();
	}
	
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		for (int i=0; i < 4; i++)
			if (source == roadLinks[i]) confd.selectObject(junction.getAllRoads()[i]);
	}
}