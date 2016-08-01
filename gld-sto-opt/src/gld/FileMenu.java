
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
 * The FileMenu. Used by both editor and simulator
 *
 * @author Group GUI
 * @version 1.0
 */

public class FileMenu extends Menu implements ActionListener
{
	Controller controller;
	
	public FileMenu(Controller c, boolean newmenu) {
		super("File");
		controller = c;
		MenuItem item;
		
		if (newmenu) {
	  	item = new MenuItem("New", new MenuShortcut(KeyEvent.VK_N));
			add(item);
			item.addActionListener(this);
		}

  	item = new MenuItem("Open...", new MenuShortcut(KeyEvent.VK_O));
  	add(item);
  	item.addActionListener(this);

  	item = new MenuItem("Save", new MenuShortcut(KeyEvent.VK_S));
  	add(item);
  	item.addActionListener(this);

    item = new MenuItem("Save as...");
    add(item);
    item.addActionListener(this);

		addSeparator();

    item = new MenuItem("Properties...");
    add(item);
    item.addActionListener(this);

		addSeparator();
    
    item = new MenuItem("Quit", new MenuShortcut(KeyEvent.VK_Q));
    add(item);
    item.addActionListener(this);
	}
	
	public void actionPerformed(ActionEvent e) {
		String s = e.getActionCommand();
			
		if (s.equals("New")) controller.newFile();
		else if (s.equals("Open...")) controller.openFile();
		else if (s.equals("Save")) controller.saveFile();
		else if (s.equals("Save as...")) controller.saveFileAs();
		else if (s.equals("Properties...")) controller.showFilePropertiesDialog();
		else if (s.equals("Quit")) controller.quit();
	}
}