
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

import java.awt.*;
import java.awt.event.*;

/**
 *
 * The HelpMenu. Used by both editor and simulator
 *
 * @author Group GUI
 * @version 1.0
 */

public class HelpMenu extends Menu implements ActionListener
{
	Controller controller;
	
	public HelpMenu(Controller c) {
		super("Help");
		controller = c;
		
		MenuItem item;
		
  	item = new MenuItem("Help", new MenuShortcut(KeyEvent.VK_H));
  	add(item);
  	item.addActionListener(this);

		addSeparator();

  	item = new MenuItem("Specifications");
  	add(item);
  	item.addActionListener(this);

    item = new MenuItem("License");
    add(item);
    item.addActionListener(this);

    item = new MenuItem("Website");
    add(item);
    item.addActionListener(this);

		addSeparator();
    
    item = new MenuItem("About");
    add(item);
    item.addActionListener(this);
	}
	
	public void actionPerformed(ActionEvent e) {
		String sel = e.getActionCommand();
		if (     sel.equals("Help"))           controller.showHelp(HelpViewer.HELP_INDEX);
		else if (sel.equals("Specifications")) controller.showHelp(HelpViewer.HELP_SPECS);
		else if (sel.equals("Website"))        controller.showHelp(HelpViewer.HELP_WEBSITE);
		else if (sel.equals("License"))        controller.showHelp(HelpViewer.HELP_LICENSE);
		else if (sel.equals("About"))          controller.showHelp(HelpViewer.HELP_ABOUT);
	}
}