
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
import gld.tools.*;
import gld.xml.*;

import java.awt.*;
import java.awt.event.*;

/**
 *
 * The dialog used to change the size of the infrastructure
 *
 * @author Group GUI
 * @version 1.0
 */

public class EditSizeDialog extends Dialog implements ActionListener, WindowListener
{
	SizePanel sizePanel;
	int width, height;
	boolean ok;
	
	/** Creates an <code>EditPropDialog</code>. */
	public EditSizeDialog(Controller c, int width, int height)
	{
		super(c, "Change size", true);
		this.setResizable(false);
		this.setSize(400, 200);
		this.addWindowListener(this);
		this.setLayout(new BorderLayout());

		sizePanel = new SizePanel(this, width, height);
		this.add(sizePanel, BorderLayout.CENTER);
		this.add(new OkCancelPanel(this), BorderLayout.SOUTH);
	}






	
	
	/*============================================*/
	/* GET, SET, ok() and show()                  */
	/*============================================*/
	
	// I & S erachter anders clash met Component.getWidth() en zo
	
	/** Returns the width entered in the textfield */
	public int getWidthI() { return width; }
	/** Sets the text of the width textfield */
	public void setWidthI(int w) { width = w; sizePanel.setWidthS(w); }
	
	/** Returns the height entered in the textfield */
	public int getHeightI() { return height; }
	/** Sets the text of the height textfield */
	public void setHeightI(int h) { height = h; sizePanel.setHeightS(h); }
	
	/** Returns true if the user clicked 'Ok' to close this dialog */
	public boolean ok() { return ok; }
	
	/** Shows the dialog */
	public void show() {
		ok = false;
		super.show();
	}











	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Ok") || e.getSource() instanceof TextField) {
 			try {
 				width = Integer.parseInt(sizePanel.getWidthS());
 				height = Integer.parseInt(sizePanel.getHeightS());
 			}
 			catch (NumberFormatException exp) {
 				new ErrorDialog(this, "You must enter an integer");
 				ok = false;
 				return;
 			}
 			ok = true;
 		}
 		else ok = false;
		hide();
 	}
 	
	public void windowClosing(WindowEvent e) { hide(); }
	public void windowActivated(WindowEvent e) { }
	public void windowClosed(WindowEvent e) { }
	public void windowDeactivated(WindowEvent e) { }
	public void windowIconified(WindowEvent e) { }
	public void windowDeiconified(WindowEvent e) { }
	public void windowOpened(WindowEvent e) { }











	/*============================================*/
	/* Panels                                     */
	/*============================================*/

  /** Panel containing the necessary components to set the infrastructure properties. */
  private static class SizePanel extends Panel
  {
  	TextField tfWidth, tfHeight;
  	EditSizeDialog esd;

		public SizePanel(ActionListener al, int width, int height)
		{
			this.esd = esd;
			GridBagLayout gridbag = new GridBagLayout();
			this.setLayout(gridbag);
			
			tfWidth = makeRow(gridbag, "Width", "" + width);
			tfHeight = makeRow(gridbag, "Height", "" + height);
			
			tfWidth.addActionListener(al);
			tfHeight.addActionListener(al);
		}
		
		private TextField makeRow(GridBagLayout gridbag, String label, String text)
  	{
			GridBagConstraints c = new GridBagConstraints();
			Label lbl;
			TextField tf;
			
			c.fill = GridBagConstraints.BOTH;
			c.weightx = 1.0;
			lbl = new Label(label);
			gridbag.setConstraints(lbl, c);
			this.add(lbl);
			c.gridwidth = GridBagConstraints.REMAINDER;
			c.weightx = 2.0;
			tf = new TextField(text, 40);
			gridbag.setConstraints(tf, c);
			this.add(tf);
			return tf;
 		}
 		
 		public String getWidthS() { return tfWidth.getText().trim(); }
 		public void setWidthS(int w) { tfWidth.setText("" + w); }
 		
 		public String getHeightS() { return tfHeight.getText().trim(); }
 		public void setHeightS(int h) { tfHeight.setText("" + h); }
	}



  /** Panel containing buttons "Ok" and "Cancel". */
  public class OkCancelPanel extends Panel
  {
		public OkCancelPanel(ActionListener al)
		{  
			this.setLayout(new FlowLayout(FlowLayout.CENTER));
			String[] labels = {"Ok", "Cancel"};
 			Button b;
 			for(int i=0; i<labels.length; i++)
 			{
 				b = new Button(labels[i]);
 				b.addActionListener(al);
 				this.add(b);
 			}
 		}
	}
}