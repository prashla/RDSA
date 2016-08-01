
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
import java.net.*;
import java.io.File;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.text.html.*;

/**
 *
 * Class used to view and browse through online help.
 *
 * @author Group GUI
 * @version 1.0
 */

public class HelpViewer extends JFrame
{
	/** Constant for use with <code>showHelp(int index)</code>, indicating to show the help index. */
	public static final int HELP_INDEX = 0;
	/** Constant for use with <code>showHelp(int index)</code>, indicating to show the specifications. */
	public static final int HELP_SPECS = 1;
	/** Constant for use with <code>showHelp(int index)</code>, indicating to show the license for this application. */
	public static final int HELP_LICENSE = 2;
	/** Constant for use with <code>showHelp(int index)</code>, indicating to show the website of this application. */
	public static final int HELP_WEBSITE = 3;
	/** Constant for use with <code>showHelp(int index)</code>, indicating to show the about page. */
	public static final int HELP_ABOUT = 4;
	
	/** The Controller that created this HelpViewer. */
	protected Controller controller;

	/**  The pane viewing the help pages. */
	protected JEditorPane viewPane;
	
	/** The base URL for all online help files. */
	protected static String base = determineBase();
	
	
	/** Creates a HelpViewer showing nothing. 
	* @param _controller The Controller to be used
	*/
	public HelpViewer(Controller _controller)
	{
		super("GLD Online Help");
		controller = _controller;
		this.setSize(820,625);
		
		// add view
		Container c = this.getContentPane();
		viewPane = new JEditorPane();
		viewPane.setEditable(false);
		viewPane.addHyperlinkListener(new HyperListener());
		JScrollPane scroller = new JScrollPane(viewPane);
		c.add(scroller);
	}
	
	/** Determines the base URL for all local help files. */
	public static String determineBase()
	{
		// determine base URL
		File file = new File("");
		String base = "file:/" + file.getAbsolutePath();		
		if(!base.endsWith("/gld"))
			base += "/gld";
		base += "/docs/";
		return base;
	}

	/** 
	* Returns the URL for the helpItem specified.
	* @param helpItem One of constants HELP_INDEX, HELP_SPECS, HELP_LICENSE, HELP_WEBSITE
	*                 and HELP_ABOUT
	*/
	public static String getHelpItem(int helpItem)
	{
		switch(helpItem)
		{
			case HELP_INDEX: return base + "index.html";
			case HELP_SPECS: return base + "specs/index.html";
			case HELP_LICENSE: return base + "license.html";
			case HELP_WEBSITE: return "http://www.students.cs.uu.nl/swp/2001/isg";
			case HELP_ABOUT: return base + "about.html";
		}
		return "";
	}

	/** 
	* Shows the help page belonging to the specified HELP_xxx constant.
	* @param helpItem One of constants HELP_INDEX, HELP_SPECS, HELP_LICENSE, HELP_WEBSITE
	*                 and HELP_ABOUT
	*/
	public void showHelp(int helpItem)
	{
		toURL(getHelpItem(helpItem));
		setVisible(true);
		requestFocus();
	}
	
	/**
	* Sets the view pane to a given URL.
	* @param s String representing the URL to show.
	*/
	protected void toURL(String s)
	{
		try
		{
			URL url = new URL(s);
			toURL(url);
		}
		catch(Exception exc) 
		{ 
			controller.showError("Couldn't view help item \"" + s + "\" : " + exc);
		}
	}

	/**
	* Sets the view pane to a given URL.
	* @param url URL to view.
	*/
	protected void toURL(URL url)
	{
		try
		{
			viewPane.setPage(url);
		}
		catch(Exception exc) { controller.showError("Couldn't view help item \"" + url + "\" : " + exc); }
	}
	
	/**
	* Simple inner class which listens to hyperlink events.
	*/
	public class HyperListener implements HyperlinkListener
	{
		/**
		* Refer to a link's destination when activated
		* @param e The HyperlinkEvent to be checked
		
		*/
		public void hyperlinkUpdate(HyperlinkEvent e)
		{
			if(e.getEventType() == HyperlinkEvent.EventType.ACTIVATED)
				toURL(e.getURL());
		}
	}
}