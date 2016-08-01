
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
import java.util.*;
import java.awt.event.*;

/**
 *
 * Basic ToolBar class.
 * Designed on Windows, so it may not look as good on other platforms...
 *
 * @author Group Joep Moritz (No, I'm really not schizo.)
 * @version 1.0
 */

public class ToolBar extends Panel
{
//	Vector tools;
	protected int separatorWidth = 8;
	protected int buttonWidth = 24;
	protected int buttonHeight = 24;
	protected int totalWidth = 0;
	
	public ToolBar() {
		super();
//		tools = new Vector(3);

		setLayout(null);
		setBackground(Color.lightGray);
	}
	
	public int getSeparatorWidth() { return separatorWidth; }
	public void setSeparatorWidth(int w) { separatorWidth = w; }
	public int getButtonWidth() { return buttonWidth; }
	public void setButtonWidth(int w) { buttonWidth = w; }
	public int getButtonHeight() { return buttonHeight; }
	public void setButtonHeight(int h) { buttonHeight = h; }
	
	
	/**
	 * Adds an IconButton to this toolbar.
	 *
	 * @param imgurl The url of the image for the new IconButton
	 * @param al The ActionListener to add to the new IconButton
	 * @param id The id for the new button
	 */
	public void addButton(String imgurl, ActionListener al, int id) {
  	Toolkit tk = Toolkit.getDefaultToolkit();
  	IconButton b = new IconButton(tk.getImage(imgurl).
  		getScaledInstance(buttonWidth - 6, buttonHeight - 6, Image.SCALE_SMOOTH), id);
  	b.setBounds(totalWidth, 0, buttonWidth, buttonHeight);
  	b.addActionListener(al);
  	add(b);
  	
  	totalWidth += buttonWidth;
	}
	
	/** Adds a separator to this toolbar */
	public void addSeparator() {
		Panel p = new Panel(null);
		p.setSize(separatorWidth, buttonHeight);
		add(p);
		
		totalWidth += separatorWidth;
	}
	
	/** Adds a component to this toolbar */
	public void addComponent(Component c) {
		add(c);
		c.setLocation(totalWidth, 0);
		totalWidth += c.getWidth();
	}
	
	public void remComponent(Component c) {
		remove(c);
		
		totalWidth -= c.getWidth();
	}
}