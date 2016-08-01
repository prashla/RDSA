
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
 * The dialog used to set <code>Infrastructure</code> properties.
 *
 * @author Group GUI
 * @version 1.0
 */

public class EditPropDialog extends Dialog
{
	TextField infranameText, authorText, commentsText;
	boolean ok;
	
	/** Creates an <code>EditPropDialog</code>. */
	public EditPropDialog(Controller c, String infraname, String author, String comments)
	{
		super(c, "Properties...", true);
		this.setResizable(false);
		this.setSize(400, 200);
		this.addWindowListener(new EPWindowListener());
		this.setLayout(new BorderLayout());

		ActionListener al = new EPActionListener();
		this.add(new PropPanel(infraname, author, comments), BorderLayout.CENTER);
		this.add(new OkCancelPanel(al), BorderLayout.SOUTH);
	}



	
	
	/*============================================*/
	/* GET, SET, ok() and show()                  */
	/*============================================*/
	
	/** Returns the infrastructure name the user entered. */
	public String getInfraname() { return infranameText.getText().trim(); }
	/** Sets the infrastructure name as shown in the dialog. */
	public void setInfraname(String s) { infranameText.setText(s); }
	/** Returns the author name the user entered. */
	public String getAuthor() { return authorText.getText().trim(); }
	/** Sets the author name as shown in the dialog. */
	public void setAuthor(String s) { authorText.setText(s); }
	/** Returns the comments the user entered. */
	public String getComments() { return commentsText.getText().trim(); }
	/** Sets the comments as shown in the dialog. */
	public void setComments(String s) { commentsText.setText(s); }
	
	/** Returns true if the user clicked 'Ok' to close this dialog. */
	public boolean ok() { return ok; }
	
	/** Shows the dialog. */
	public void show()
	{
		ok = false;
		super.show();
	}
	






	/*============================================*/
	/* Listeners                                  */
	/*============================================*/
	
	/** Listens to the buttons of the dialog. */
	public class EPActionListener implements ActionListener
	{
		public void actionPerformed(ActionEvent e)
		{
			String sel = ((Button)e.getSource()).getLabel();
			if(sel.equals("Ok"))
				ok = true;
			hide();
		}
	}
	
	
 	/** Listens to the <code>EditPropDialog</code> window. */
 	public class EPWindowListener extends WindowAdapter
 	{
	  public void windowClosing(WindowEvent e) 
  	{
   		hide();
  	}
  }






	/*============================================*/
	/* Panels                                     */
	/*============================================*/

  /** Panel containing the necessary components to set the infrastructure properties. */
  public class PropPanel extends Panel
  {
		public PropPanel(String infraname, String author, String comments)
		{ 
			GridBagLayout gridbag = new GridBagLayout();
			this.setLayout(gridbag);
			
			infranameText = makeRow(gridbag, "Title", infranameText, infraname);
			authorText = makeRow(gridbag, "Author", authorText, author);
			commentsText = makeRow(gridbag, "Comments", commentsText, comments);
		}
		
		private TextField makeRow(GridBagLayout gridbag, String label, TextField textField, String text)
  	{
			GridBagConstraints c = new GridBagConstraints();
			Label lbl;
			
			c.fill = GridBagConstraints.BOTH;
			c.weightx = 1.0;
			lbl = new Label(label);
			gridbag.setConstraints(lbl, c);
			this.add(lbl);
			c.gridwidth = GridBagConstraints.REMAINDER;
			c.weightx = 2.0;
			textField = new TextField(text, 40);
			gridbag.setConstraints(textField, c);
			this.add(textField);
			return textField;
 		}
	}

  /** Panel containing buttons "Ok" and "Cancel". */
  public class OkCancelPanel extends Panel
  {
		public OkCancelPanel(ActionListener action)
		{  
			this.setLayout(new FlowLayout(FlowLayout.CENTER));
			String[] labels = {"Ok", "Cancel"};
 			Button b;
 			for(int i=0; i<labels.length; i++)
 			{
 				b = new Button(labels[i]);
 				b.addActionListener(action);
 				this.add(b);
 			}
 		}
	}

}