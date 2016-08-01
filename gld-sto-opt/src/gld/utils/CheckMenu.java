
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

package gld.utils;

import java.awt.*;
import java.awt.event.*;
import java.util.Vector;

/**
 * This class is used to construct menus containing only CheckboxMenuItems
 *
 * @author Joep Moritz
 * @version 1.0
 */

public class CheckMenu extends Menu implements ItemSelectable
{
	Vector listeners;
	Vector selectedItems;
	boolean allowMultipleSelections;
	int selectedIndex;
	
	/**
	 * Creates a new CheckMenu. The new menu does not allow multiple selections.
	 *
	 * @param name The name of the menu
	 * @param items The CheckboxMenuItems that are added to this menu
	 */
	public CheckMenu(String name, String[] items) {
		this(name, items, false);
	}
	
	/**
	 * Creates a new CheckMenu.
	 *
	 * @param name The name of the new menu
	 * @param items The CheckboxMenuItems that are added to this menu
	 * @param allowmul Allows multiple selections if true
	 */
	public CheckMenu(String name, String[] items, boolean allowmul) {
		super(name);
		
		listeners = new Vector(1);
		selectedItems = new Vector(1);
		allowMultipleSelections = allowmul;
		
		Listener lis = new Listener(this);
		
		CheckboxMenuItem citem;
		for (int i=0; i < items.length; i++) {
			citem = new CheckboxMenuItem(items[i]);
			add(citem);
			citem.addItemListener(lis);
		}
	}
	
	/** Adds given item listener */
	public void addItemListener(ItemListener il) { listeners.add(il); }
	/** Removes given item listener */
	public void removeItemListener(ItemListener il) { listeners.remove(il); }
	/** Returns an array of the objects currently selected */
	public Object[] getSelectedObjects() { return selectedItems.toArray(new CheckboxMenuItem[1]); }
	/** Returns an array of all CheckboxMenuItems in this Checkmenu */
	public CheckboxMenuItem[] getItems() {
		CheckboxMenuItem[] citems = new CheckboxMenuItem[getItemCount()];
		for(int i=0; i<citems.length; i++)
			citems[i] = (CheckboxMenuItem)getItem(i);
		return citems;
	}
	
	/** Deselects all items. This will send an itemStateChanged message to all listeners */
	public void deselectAll() {
		for (int i=0; i < selectedItems.size(); i++) {
			((CheckboxMenuItem)selectedItems.get(i)).setState(false);
		}
		selectedItems.clear();
	}
	
	/** Deselects given item. This will send an itemStateChanged message to all listeners */
	public void deselect(CheckboxMenuItem cmi) {
		if (allowMultipleSelections) {
			cmi.setState(false);
			selectedItems.remove(cmi);
		}
	}
	
	/** Deselects the item with given index. This will send an itemStateChanged message to all listeners */
	public void deselect(int i) {
		deselect((CheckboxMenuItem)getItem(i));
	}
	
	/** Selects the given item. This will send an itemStateChanged message to all listeners */
	public void select(CheckboxMenuItem cmi) {
		select(getIndex(cmi));
	}
	
	/** Selects the item with given index. This will send an itemStateChanged message to all listeners */
	public void select(int i) {
		if (!allowMultipleSelections) {
			deselectAll();
			selectedIndex = i;
		}
		CheckboxMenuItem cmi = ((CheckboxMenuItem)getItem(i));
		cmi.setState(true);
		selectedItems.add(cmi);
	}
	
	/**
	 * Returns the index of the given item.
	 * Returns -1 if the item is not part of this menu.
	 */
	public int getIndex(CheckboxMenuItem cmi) {
		int nritems = getItemCount();
		for (int i=0; i < nritems; i++) if (getItem(i) == cmi) return i;
		return -1;
	}

	/**
	 * Returns the index of the currently selected item.
	 * Returns -1 if no item is selected, or if this is a multiple selections menu
	 */
	public int getSelectedIndex() { return selectedIndex; }
	
	/**
	 * Returns the currently selected item.
	 * Returns null if no item is selected, or if this is a multiple selections menu
	 */
	public CheckboxMenuItem getSelectedItem() {
		if (!allowMultipleSelections && selectedItems.size() > 0) {
			return (CheckboxMenuItem)selectedItems.get(0);
		}
		return null;
	}


	private class Listener implements ItemListener
	{
		private CheckMenu menu;
		public Listener(CheckMenu cm) { menu = cm; }
		public void itemStateChanged(ItemEvent e) {
			boolean sendmsg = true;
			CheckboxMenuItem cmi = (CheckboxMenuItem)e.getItemSelectable();
		
			if (e.getStateChange() == ItemEvent.SELECTED) {
				if (!allowMultipleSelections) {
					deselectAll();
					selectedIndex = getIndex(cmi);
				}
				if (!selectedItems.contains(cmi)) selectedItems.add(cmi);
			}
			else {
				if (!allowMultipleSelections) {
					sendmsg = false;
					if (selectedItems.contains(cmi)) select(cmi);
				}
				else {
					selectedItems.remove(cmi);
				}
			}

			if (sendmsg) {
				ItemEvent ie = new ItemEvent(menu, e.getID(), e.getItem(), e.getStateChange());
				for (int i=0; i < listeners.size(); i++)
					((ItemListener)listeners.get(i)).itemStateChanged(ie);
			}
		}
	}
}