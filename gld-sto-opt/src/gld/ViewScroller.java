
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

import gld.infra.*;

import java.awt.*;
import java.util.*;
import java.awt.event.*;

/**
 *
 * @author  Group GUI
 * @version 1.0
 *
 * A <code>ViewScroller</code> is a component that lets the user scroll over
 * a <code>View</code> using scrollbars.
 *
 */

public class ViewScroller extends ScrollPane implements AdjustmentListener, ComponentListener
{
	protected View view;

	/**
	 * Create a new <code>ViewScroller</code> for a given <code>Container</code> and <code>View</code>.
	 *
	 * @param parent the parent <code>Container</code> for this <code>ViewScroller</code>
	 * @param view  the <code>View</code> for this <code>ViewScroller</code>
	 */
	public ViewScroller(View view)
	{
		super(ScrollPane.SCROLLBARS_ALWAYS);
		setup(view);
	}


	/**
	 * Create a new <code>ViewScroller</code> for a given <code>Container</code> and <code>View</code>.
	 *
	 * @param parent   the parent <code>Container</code> for this <code>ViewScroller</code>
	 * @param view    the <code>View</code> for this <code>ViewScroller</code>
	 */
	public ViewScroller(View view, boolean scrollbars)
	{
		super(scrollbars ? ScrollPane.SCROLLBARS_ALWAYS : ScrollPane.SCROLLBARS_NEVER);
		setup(view);
	}

	private void setup(View v)
	{
		view = v;
		getHAdjustable().addAdjustmentListener(this);
		getVAdjustable().addAdjustmentListener(this);
		add(view);
		addComponentListener(this);
	}








	/** Called when infrastructure is resized */
	public void resizeInfra(Dimension infras) {
		Point p = view.toInfra(view.getViewportPosition());
		view.resizeInfra(infras);
		scrollTo(p);
		doLayout();
		scrollTo(p);
	}

	/** Centers the view */
	public void center() {
		setScrollPosition(new Point(
			(int)(view.getWidth() / 2 - view.getViewportSize().width / 2),
			(int)(view.getHeight() / 2 - view.getViewportSize().height / 2)));
	}
	
	/** Scrolls the view to the specified infrastructure point */
	public void scrollTo(Point p) {
		Point vp = view.toView(p);
		setScrollPosition(vp);
	}
	
	/** Centers the view on the specified infrastructure point */
	public void center(Point p) {
		p = view.toView(p);
		p.x -= (int)(view.getViewportSize().width / 2);
		p.y -= (int)(view.getViewportSize().height / 2);
		setScrollPosition(p);
	}
	
	/** Returns the current center point in infrastructure coordinate space */
	public Point getCurrentCenter() {
		Point p = getScrollPosition();
		p.x += (int)(getViewportSize().width / 2);
		p.y += (int)(getViewportSize().height / 2);
		return view.toInfra(p);
	}
	
	/** Zooms the view to the specified index */
	public void zoomTo(int index) {
		Point cc = getCurrentCenter();
		view.zoomTo(index);
		center(cc);
		doLayout();
		center(cc);
	}
	
	/** Zooms the view in, centering on the given point */
	public void zoomIn(Point p) {
		view.zoomIn();
		center(p);
		doLayout();
		center(p);
	}
	
	/** Zooms the view out, centering on the given point */
	public void zoomOut(Point p) {
		view.zoomOut();
		center(p);
		doLayout();
		center(p);
	}		









	public void adjustmentValueChanged(AdjustmentEvent e) {
		view.scrollViewport(getScrollPosition());
	}
	
	public void componentHidden(ComponentEvent e) { }
	public void componentMoved(ComponentEvent e) { }
	public void componentResized(ComponentEvent e) {
		view.resizeViewport(getViewportSize());
		doLayout();
	}
	public void componentShown(ComponentEvent e) { }
}