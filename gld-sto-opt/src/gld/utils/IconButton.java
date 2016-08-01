
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

import java.awt.Button;
import java.awt.Graphics;
import java.awt.Image;

/**
 *
 */
public class IconButton extends Button
{
	/** The image on this IconButton */
	protected Image img = null;
	protected int id = 0;
	
	/** Create a new IconButton without an Image */
	public IconButton() {
		super();
	}
	
	/** Create a new IconButton with a given Image */
	public IconButton(Image i) {
		super();
		img = i;
	}
	
	/** Create a new IconButton with a given Image and Id */
	public IconButton(Image i, int idn) {
  	super();
		img = i;
		id = idn;
	}




	/**
	 * Draw this IconButton
	 * @param g The Graphics object to draw onto
	 */
	public void paint(Graphics g)
	{
		super.paint(g);
		if (img != null)
		{
			g.drawImage(img,(int)(this.getWidth() / 2) - (int)(img.getWidth(this) / 2),
											(int)(this.getHeight() / 2) - (int)(img.getHeight(this) / 2),
											this);
		}
	}
	
	/** Get the image for this IconButton */
	public Image getImage() { return img; }
	/** Set the image for this IconButton */
	public void setImage(Image i) { img = i; }
	
	/** Returns the id of this button */
	public int getId() { return id; }
	/** Sets the Id of this button */
	public void setId(int i) { id = i; }
}
		
		