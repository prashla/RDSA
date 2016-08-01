
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

import gld.edit.EditController;
import gld.edit.EditModel;
import gld.sim.SimController;
import gld.sim.SimModel;

/**
 *
 * The general class which is used to start up the editor or simulator
 *
 * @author Group Model
 * @version 1.0
 */
public class GLDStarter
{	protected boolean splashScreen=false,loadFile=false,noFurtherOptions=false;
	protected String filename="";
	protected String[] params;
	protected static int type=-1;
	protected final static int EDITOR=0;
    public static final int SIMULATOR=1; 
										
	/** Make a new GLDStarter
	  * @param params The command line parameters with which the real main
	  *        class was loaded.
	  * @param type The type of the controller
	 */
	public GLDStarter (String[] params,int startType)
	{	type=startType;
		this.params=params;
	}									
	
	/** Indicates what kind of program is running in this VM (EDITOR or
	 * SIMULATOR). The method returns -1 if the program is not yet running.
	 */ 
	public static int getProgramModus ()
	{	return type;
	}
	
	/** Process the command line parameters that where specified when
	  * the program was started. 
	 */
 	public void processParams () 
	{	if (params.length>0)
		{	for (int i=0;i<params.length;i++)
			{	if (params[i].startsWith("-") && ! noFurtherOptions)
				{	processOption (params[i]);
				}
				else if (! loadFile)
				{	filename=params[i];
					loadFile=true;
				}
				else
				{	illegalParametersError();
				}
			}
		}
   }
	
	/** Each option parameter (parameter that begins with a dash) is processed
	  * by this method
	  * @param The option, including dash sign
	 */
	public void processOption (String option)
	{	if ("-splash".equals(option))
		{	System.out.println("Splash!");
			splashScreen=true;
		}
		else if ("-nosplash".equals(option))
		{	System.out.println("No Splash!");
			splashScreen=false;
		}
		else if ("--".equals(option))
			noFurtherOptions=true;
		else
			illegalParametersError();
	}
	
	/** This method is called when a command line parameter is encountered that
	  * we do not understand. It prints an error message.
	 */
	public void illegalParametersError ()
	{	System.out.println("Illegal parameters");
		System.out.println("Usage : java "+getStarterName()+" [--] [filename] [-splash]");
		System.exit (1);
	}
	
	public void start()
	{	processParams();
		Controller controller=getController();
		if (loadFile)
			controller.tryLoad(filename);
	}
	
	/** Method which gives the right controller for starting */
	public Controller getController () 
	{	switch (type)
		{ 	case EDITOR :  return new EditController (new EditModel(),splashScreen);
			case SIMULATOR:return new SimController (new SimModel(),splashScreen);
		}
		System.out.println
			("GLDStarter was called with unknown controller type "+type+
			 ". Cannot get new instance of Controller.");
		System.exit(1);	 
		return null; // Dummy to satisfy the compiler
	}
	
	/** Method which gives the right main class name for error messages */
	public String getStarterName () 
	{	switch (type)
		{ 	case EDITOR   : return "GLDEdit";
			case SIMULATOR: return "GLDSim";
		}
		System.out.println
			("GLDStarter was called with unknown controller type "+type+
			 ". Cannot get name of caller.");
		System.exit(1);	 
		return null; // Dummy to satisfy the compiler
	}

}
