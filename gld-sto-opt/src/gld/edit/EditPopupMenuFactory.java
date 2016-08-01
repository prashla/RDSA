
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

package gld.edit;

import gld.*;
import gld.infra.*;

import java.awt.MenuItem;
import java.awt.MenuShortcut;
import java.awt.PopupMenu;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.lang.Exception;


/**
 *
 * Factory for creating popup menus for editor
 *
 * @author Group GUI
 * @version 1.0
 */

public class EditPopupMenuFactory
{
	protected EditController controller;
	
	public EditPopupMenuFactory(EditController con) {
		controller = con;
	}
	
	/**
	 * Creates a right-click PopupMenu for the given object.
	 * A listener is added to the menu as well.
	 */
	public PopupMenu getPopupMenuFor(Selectable obj) throws PopupException
	{
		if (obj instanceof Node) return getNodeMenu((Node)obj);
		if (obj instanceof Road) return getRoadMenu((Road)obj);
		if (obj instanceof Drivelane) return getDrivelaneMenu((Drivelane)obj);
		throw new PopupException("Unknown object type");
	}



	protected PopupMenu getNodeMenu(Node n) throws PopupException
	{
		return getGenericMenu(new DefaultGUIObjectListener());
	}
	
	protected PopupMenu getRoadMenu(Road r) throws PopupException
	{
		return getGenericMenu(new DefaultGUIObjectListener());
	}
	
	protected PopupMenu getDrivelaneMenu(Drivelane l) throws PopupException
	{
		return getGenericMenu(new DefaultGUIObjectListener());
	}
	
	
	protected PopupMenu getGenericMenu(PopupMenuListener pml)
	{
		PopupMenu menu = new PopupMenu();
  	MenuItem item = new MenuItem("Delete", new MenuShortcut(KeyEvent.VK_DELETE));
  	item.addActionListener(pml);
  	menu.add(item);
  	
  	menu.add(new MenuItem("-"));
  	
  	item = new MenuItem("Properties...", new MenuShortcut(KeyEvent.VK_ENTER));
  	item.addActionListener(pml);
  	menu.add(item);

  	return menu;
	}
	

	protected class DefaultGUIObjectListener implements PopupMenuListener
	{
		public void actionPerformed(ActionEvent e) {
			String s = e.getActionCommand();
			if (s.equals("Delete")) controller.deleteSelection();
			else if (s.equals("Properties..."))
				controller.showConfigDialog();
		}
	}







	protected static interface PopupMenuListener extends ActionListener { }
}