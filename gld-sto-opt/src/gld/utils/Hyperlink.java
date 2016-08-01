
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
 *
 */
public class Hyperlink extends Component
{
	/** Indicates that the label should be left justified. */
	public static final int LEFT = 1;
	/** Indicates that the label should be centered. */
	public static final int CENTER = 2;
	/** Indicates that the label should be right justified. */
	public static final int RIGHT = 3;

	protected String label;
	protected boolean underline;
	protected int alignment;

	protected Vector listeners;
	protected Rectangle textBounds;

	/**
	 * Constructs an empty hyperlink.
	 */
	public Hyperlink() {
		this("unlabeled", LEFT);
	}
	
	/**
	 * Constructs a new hyperlink
	 * with the specified string of text, left justified.
	 *
	 * @param text The lavel of the new hyperlink.
	 */
	public Hyperlink(String text) {
		this(text, LEFT);
	}
	
	/**
	 * Constructs a new hyperlink that presents
	 * the specified string of text with the specified alignment.
	 *
	 * @param text The label of the new hyperlink.
	 * @param al The alignment of the label.
	 */
	public Hyperlink(String text, int al) {
		super();
		setForeground(Color.blue);
		textBounds = null;
		listeners = new Vector(1);

		Listener listener = new Listener();
		addFocusListener(listener);
		addKeyListener(listener);
		addMouseListener(listener);
		addMouseMotionListener(listener);

		label = text;
		underline = true;
		alignment = al;
	}

	/** Returns the text of this hyperlink. */
	public String getText() { return label; }
	/** Sets the text of this hyperlink to the specified text. */
	public void  setText(String text) { label = text; repaint(); }

	/** Checks if the text of this hyperlink is underlined. */
	public boolean getUnderline() { return underline; }
	/** Underlines the text of this hyperlink if the specified value is true. */
	public void setUnderline(boolean ul) { underline = ul; repaint(); }

	/** Returns the alignment of this hyperlink. */
	public int getAlignment() { return alignment; }
	/** Sets the alignment of this hyperlink to the specified alignment. */
	public void setAlignment(int al) { alignment = al; repaint(); }

	/**
	 * Adds the specified action listener to receive action events from this hyperlink.
	 * If l is null, no exception is thrown and no action is performed.
	 *
	 * @param l The action listener to add.
	 */
	public void addActionListener(ActionListener l) {
		if (l != null) {
			listeners.add(l);
			enableEvents(AWTEvent.MOUSE_EVENT_MASK | AWTEvent.KEY_EVENT_MASK);
		}
	}
	
	/**
	 * Removes the specified action listener so that it no longer receives action events from this hyperlink.
	 * Action events occur when a user releases the left mouse button
	 * when the mouse cursor is over this hyperlink.
	 * If l is null, no exception is thrown and no action is performed.
	 *
	 * @param l The action listener to remove.
	 */
	public void removeActionListener(ActionListener l) {
		if (l != null) {
			listeners.remove(l);
			if (listeners.isEmpty())
				disableEvents(AWTEvent.MOUSE_EVENT_MASK | AWTEvent.KEY_EVENT_MASK);
		}
	}
	
	/**
	 * Removes all action listeners.
	 */
	public void removeAllActionListeners() {
		listeners.clear();
	}
	
	/** Returns an array of the action listeners. */
	public ActionListener[] getActionListeners() {
		return (ActionListener[])listeners.toArray(new ActionListener[1]);
	}
	

	/** Returns true to allow hyperlinks to receive focus. */
	public boolean isFocusTraversable() { return false; }
	
	
	
	private class Listener implements FocusListener, MouseListener, KeyListener, MouseMotionListener
	{
		/**
		 * Dispatches ActionEvents to the registered listeners
		 * and requests focus if the user clicked this hyperlink.
		 * 
		 * @param e The mouse event.
		 */
		public void mouseClicked(MouseEvent e) {
			if ((e.getModifiers() & MouseEvent.BUTTON1_MASK) == MouseEvent.BUTTON1_MASK && isEnabled()) {
				requestFocus();
				ActionEvent ae = new ActionEvent(e.getSource(), e.getID(), label, e.getModifiers());
				for (int i=0; i < listeners.size(); i++)
					((ActionListener)(listeners.get(i))).actionPerformed(ae);
			}
		}
	
		/**
		 * Dispatches ActionEvents to the registered listeners
		 * if the user pressed the return or space key.
		 * 
		 * @param e The key event.
		 */
		public void keyReleased(KeyEvent e) {
			int key = e.getKeyCode();
			if ((key == KeyEvent.VK_SPACE || key == KeyEvent.VK_ENTER) && isEnabled()) {
				ActionEvent ae = new ActionEvent(e.getSource(), e.getID(), label, e.getModifiers());
				for (int i=0; i < listeners.size(); i++)
					((ActionListener)(listeners.get(i))).actionPerformed(ae);
			}
		}
		
		/**
		 * Changes the cursor into a hand if the cursor position lies within the bounds of the text
		 */
		public void mouseMoved(MouseEvent e) {
			if (textBounds != null && textBounds.contains(e.getPoint()) && isEnabled())
				setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
			else
				setCursor(null);
		}

		/** Calls repaint to show this hyperlink received focus. */
		public void focusGained(FocusEvent e) { repaint(); }
		/** Calls repaint to show this hyperlink lost focus. */
		public void focusLost(FocusEvent e) { repaint(); }
		
		public void mousePressed(MouseEvent e) { }
		public void mouseReleased(MouseEvent e) { }
		public void mouseEntered(MouseEvent e) { }
		public void mouseExited(MouseEvent e) { }
		public void mouseDragged(MouseEvent e) { }

		public void keyPressed(KeyEvent e) { }
		public void keyTyped(KeyEvent e) { }
	}
		


	/**
	 * Paints this hyperlink on the given Graphics object.
	 * It uses the current font and color of the Graphics object
	 * to draw the label and possible line under the label.
	 * You should use Component.setForeground() and Component.setFont()
	 * to change the font and/or color that is used to paint the label.
	 * See the AWT docs about the Component class for more information.
	 *
	 * @param g The Graphics object to paint this hyperlink on.
	 */
	public void paint(Graphics gr)
	{
		Graphics2D g = (Graphics2D)gr;
		
		FontMetrics fm = g.getFontMetrics(g.getFont());
		int ascent = fm.getAscent();
		int descent = fm.getDescent();
		int height = ascent + descent;
		int width = fm.stringWidth(label);
		int x = 0, y = (int)((getHeight() - height) / 2);

		if (alignment == CENTER) x = (int)((getWidth() - width) / 2);
		if (alignment == RIGHT)  x = getWidth() - width - 1;
		
		g.drawString(label, x, y + ascent - 1);
		if (underline) g.drawLine(x, y + height - 1, x + width, y + height - 1);
		
		if (hasFocus()) {
			float[] pattern = { 1.0f, 1.0f };
			g.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, pattern, 0.0f));
			g.setColor(Color.black);
			g.drawRect(x, y, width, height);
		}

		textBounds = new Rectangle(x, y, width, height);
	}
	
	public void paintAll(Graphics g) {
		paint(g);
	}
}