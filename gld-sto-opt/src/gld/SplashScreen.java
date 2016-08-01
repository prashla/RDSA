
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
/*
 *The Splash screen for the program, which gets the focus on starting up the application
 */

public class SplashScreen extends Window
{	private static final int WIDTH=550,HEIGHT=424;
	Image image;
	
	public SplashScreen(Controller parent) {
		super(parent);
  		Toolkit tk = Toolkit.getDefaultToolkit();
  		image = tk.createImage("gld/images/splash.jpg");
  		requestFocus();
		MediaTracker mt=new MediaTracker(this);
		mt.addImage(image,0);
		try {mt.waitForID(0);} catch (InterruptedException e) {}
	}
	
	public void paint(Graphics g) {
		setSize(WIDTH,HEIGHT);
		g.setColor(Color.white);
		g.fillRect(0, 0, WIDTH, HEIGHT);
		g.drawImage(image, 0, 0, this);
		toFront();
	}
}
