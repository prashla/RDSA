
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
import java.util.*;

/**
 * Selections can be made in the simulator and editor
 *
 * @author Group GUI
 * @version 1.0
 */

public class Selection extends Observable implements Overlay
{
	LinkedList selectedObjects;
	View view;
	Class[] selectionFilter = { };
	SelectionStarter starter;
	
	public Selection(View v) {
		this(v, null, null);
	}
	public Selection(View v, Class[] sf) {
		this(v, sf, null);
	}
	public Selection(View v, Class[] sf, SelectionStarter ss) {
		view = v;
		selectionFilter = sf;
		starter = ss;
		selectedObjects = new LinkedList();
	}
	
	/** Creates a new Selection object that is a clone of the given Selection object */
	public Selection(Selection copy) {
		view = copy.getView();
		selectedObjects = (LinkedList)copy.getSelectedObjects().clone();
		selectionFilter = copy.getSelectionFilter();
		starter = copy.getSelectionStarter();
	}
	
	public int overlayType() { return 2; }
	
	public void paint(Graphics gr) throws GLDException {
		Graphics2D g = (Graphics2D)gr;
		Selectable obj = null;
		//g.setXORMode(Color.darkGray);
		//g.setColor(Color.lightGray);
		g.setPaintMode();
		g.setColor(new Color(0, 0.8f, 0, 0.4f));
		for (ListIterator it = selectedObjects.listIterator(); it.hasNext(); )
			g.fill(((Selectable)it.next()).getComplexBounds());
	}
	
	/** Returns the view associated with this selection */
	public View getView() { return view; }

	/** Returns the current selection starter */
	public SelectionStarter getSelectionStarter() { return starter; }
	/** Sets the current selection starter */
	public void setSelectionStarter(SelectionStarter ss) { starter = ss; }
	
	/** Returns the current selection filter */
	public Class[] getSelectionFilter() { return selectionFilter; }
	/** Sets the current selection filter */
	public void setSelectionFilter(Class[] sf) { selectionFilter = sf; }
	
	/** Returns the list of all selected objects */
	public LinkedList getSelectedObjects() { return selectedObjects; }
	/** Sets the list of all selected objects */
	public void setSelectedObjects(LinkedList l) {
		selectedObjects = l;
		setChanged();
		notifyObservers();
	}
	/** Return the number of currently selected object */
	public int getNumSelectedObjects() { return selectedObjects.size(); }
	/** Checks if the selection is empty */
	public boolean isEmpty() { return selectedObjects.isEmpty(); }
	
	/** Sets this selection to the given selection */
	public void setSelection(Selection s) { selectedObjects = (LinkedList)s.getSelectedObjects().clone(); }
	
	/**
	 * Removes all objects from the list of selected objects.
	 */
	public void deselectAll() {
		selectedObjects = new LinkedList();
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Removes a Selectable from the list of selected objects.
	 *
	 * @param o The Selectable to remove from the list
	 */
	public void deselect(Selectable o) {
		selectedObjects.remove(o);
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Adds all objects to the list of selected objects.
	 *
	 * @param o The Selectable that contains a recursive list of all objects to add
	 */
	public void selectAll() {
		selectedObjects = new LinkedList();
		selectAllObjects(starter);
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Adds all children (recursively) of a SelectionStarter to the list of selected objects.
	 *
	 * @param o The Selectable of which to add the children to the list
	 */
	private void selectAllObjects(SelectionStarter ss) {
		Selectable child = null;
		for (Enumeration e = ss.getChildren(); e.hasMoreElements(); ) {
			child = (Selectable)e.nextElement();
			if (child.isSelectable()) selectedObjects.add(child);
			selectAllObjects(child);
		}
	}
	
	/**
	 * Creates a new selection containing the given Selectable
	 *
	 * @param o The Selectable the new selection should contain
	 */
	public void newSelection(Selectable o) {
		selectedObjects = new LinkedList();
		selectedObjects.add(o);
		setChanged();
		notifyObservers();
	}

	/**
	 * Creates a new selection holding the Selectable located at
	 * given point. If no such Selectable exists, the selection is emptied.
	 *
	 * @param p The point of the Selectable to add to the new selection
	 */
	public void newSelection(Point p) {
		selectedObjects = new LinkedList();
		Selectable obj = selectBest(selectItems(starter, p), p);
		if (obj != null) selectedObjects.add(obj);
		setChanged();
		notifyObservers();
	}

	/**
	 * Adds the Selectable located at given point to the selection
	 *
	 * @param p The location of the Selectable to add
	 */
	public void addToSelection(Point p) {
		Selectable obj = selectBest(selectItems(starter, p), p);
		if (obj != null) {
			selectedObjects.add(obj);
			setChanged();
			notifyObservers();
		}
	}
	
	/**
	 * Removes the Selectable located at given point from the selection.
	 * If no such Selectable is selected, adds the Selectable located at
	 * given point to the selection.
	 *
	 * @param p The point the Selectable is located at
	 */
	public void invertWithSelection(Point p) {
		LinkedList items = selectItems(starter, p);
		Selectable obj = null;
		for (ListIterator it = items.listIterator(); it.hasNext(); ) {
			obj = (Selectable)it.next();
			if (selectedObjects.contains(obj)) {
				selectedObjects.remove(obj);
				setChanged();
				notifyObservers();
				return;
			}
		}
		obj = selectBest(items, p);
		if (obj != null) {
			selectedObjects.add(obj);
			setChanged();
			notifyObservers();
		}
	}
	
	
	/**
	 * Adds all Selectables inside given rectangle to selection.
	 *
	 * @param r The rectangle enclosing Selectables to add
	 */
	public void addToSelection(Rectangle r) {
		LinkedList items = selectItems(starter, r);
		Selectable obj = null;
		for (ListIterator it = items.listIterator(); it.hasNext(); ) {
			obj = (Selectable)it.next();
			if (!selectedObjects.contains(obj))
				selectedObjects.add(obj);
		}
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Creates a new selection holding all Selectables inside given rectangle
	 *
	 * @param r The rectangle enclosing Selectables
	 */
	public void newSelection(Rectangle r) {
		selectedObjects = selectItems(starter, r);
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Returns the closest object selected by given point using given SelectionStarter
	 *
	 * @param sf The selection filter to aply
	 * @param ss The selection starter
	 * @param p The point
	 */
	public static Object selectObject(Class[] sf, SelectionStarter ss, Point p)
	{
		Selection s = new Selection(null, sf, ss);
		s.newSelection(p);
		if (s.isEmpty()) return null;
		return s.getSelectedObjects().getFirst();
	}
	
	/**
	 * Selects the Selectable from given list that is most likely
	 * to be targeted by given point.
	 *
	 * @param items The list of Selectables to choose from
	 * @param p The targeted point.
	 * @return The Selectable from list closest to given point
	 */
	private Selectable selectBest(LinkedList items, Point p) {
		Selectable obj = null, best = null;
		int dist = Integer.MAX_VALUE;
		
		for (ListIterator it = items.listIterator(); it.hasNext(); ) {
			obj = (Selectable)it.next();
			if (!selectedObjects.contains(obj) && obj.getDistance(p) < dist) {
				dist = obj.getDistance(p);
				best = obj;
			}
		}
		return best;
	}
	
	/**
	 * Finds out if given Selectable is located inside given rectangle
	 *
	 * @param o Selectable to check
	 * @param r Rectangle to check
	 * @return True if o is located inside r. False otherwise.
	 */
	private boolean selects(Rectangle r, Selectable o) {
		return o.isSelectable() && r.contains(o.getSelectionPoint()) && inFilter(o);
	}
	
	/**
	 * Finds out if given point is located inside bounding box of given Selectable
	 *
	 * @param p Point to check
	 * @param o Selectable to check
	 * @return True if p is located inside bounds of o. False otherwise.
	 */
	private boolean selects(Point p, Selectable o) {
		return o.isSelectable() && o.getComplexBounds().contains(p) && inFilter(o);
	}
	
	/**
	 * Finds out if given object is an instance of any of the classes in the current selection filter
	 *
	 * @param o The Object to check
	 * @return True if given object is an instance of at least 1 of the classes in the current filter
	 */
	private boolean inFilter(Object o) {
		for (int i=0; i < selectionFilter.length; i++)
			if (selectionFilter[i].isInstance(o)) return true;
		return false;
	}
	
	
	/**
	 * Creates a list of all Selectables selected by given point.
	 * Uses the @ref selects(Point p, Selectable o) function, and the
	 * getChildren() functionality of SelectionStarter.
	 *
	 * @param s SelectionStarter to start search at.
	 * @param p Point that selects objects.
	 * @return a list of all Selectables that are selected by given point
	 * @see Selectable#getChildren()
	 * @see Selection#selects(Point p, Selectable o)
	 */
	private LinkedList selectItems(SelectionStarter ss, Point p)
	{
		LinkedList list = new LinkedList();
		Selectable child;
		if (ss.hasChildren()) {
			for (Enumeration e = ss.getChildren(); e.hasMoreElements(); ) {
				child = (Selectable)e.nextElement();
				if (selects(p, child)) list.add(child);
				list.addAll(selectItems(child, p));
			}
		}
		return list;
	}

	/**
	 * Creates a list of all Selectables selected by given rectangle.
	 * Uses the @ref selects(Rectangle r, Selectable o) function, and
	 * the getChildren() functionality of Selectable.
	 *
	 * @param o SelectionStarter to start search at.
	 * @param p Point that selects objects.
	 * @return a list of all Selectables that are selected by given point
	 * @see Selectable#getChildren()
	 * @see Selection#selects(Rectangle r, Selectable o)
	 */
	private LinkedList selectItems(SelectionStarter ss, Rectangle r)
	{
		LinkedList list = new LinkedList();
		Selectable child;
		if (ss.hasChildren()) {
			for (Enumeration e = ss.getChildren(); e.hasMoreElements(); ) {
				child = (Selectable)e.nextElement();
				if (selects(r, child)) list.add(child);
				list.addAll(selectItems(child, r));
			}
		}
		return list;
	}
}