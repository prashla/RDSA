
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

package gld.infra;

import java.awt.*;
import java.util.*;
import java.io.*;
import java.applet.*;

import gld.*;
import gld.xml.*;

/**
 * 
 * @author Group Datastructures
 * @version 1.0
 *
 * This class implements customizable roadusers.
 * It can represent cars, busses and bicycles.
 */

public class CustomRoaduser extends Roaduser
{
  protected int vehicle;
  protected int driver;
  protected int[] passengers = { };


	public CustomRoaduser(Node start, Node dest, int pos) {
		super(start, dest, pos);
	}
	
	/** Empty constructor for loading */
	public CustomRoaduser() { }

	/** Returns the ID of the vehicle. */	
	public int getVehicle() { return vehicle; }
	/** Returns the name of thevehicle. */
	public String getVehicleName() { return CustomFactory.getVehicleName(vehicle); }
	/** Sets the ID of the vehicle. */
	public void setVehicle(int v) { vehicle = v; }

	/** Returns the ID of the driver. */
	public int getDriver() { return driver; }
	/** Returns the name of the driver. */
	public String getDriverName() { return CustomFactory.getPersonName(driver); }
	/** Sets the ID of the driver. */
	public void setDriver(int d) { driver = d; }

	/** Returns an array of passenger IDs. */
	public int[] getPassengers() { return passengers; }
	public int getNumPassengers() { return (passengers.length < 1 ? 1 : passengers.length); }
	/** Sets the list of passengers. */
	public void setPassengers(int[] p) { passengers = p; }

	/** Returns the speed of this custom. */
	public int getSpeed() { return CustomFactory.getSpeed(this); }
	/** Returns the length of this custom. */
	public int getLength() { return CustomFactory.getLength(this); }
	/** Returns the roaduser type of this custom. */
	public int getType() { return CustomFactory.getType(this); }
	/** Returns the name of this custom. */
	public String getName() { return CustomFactory.getName(this); }
	/** Returns the description of this custom. */
	public String getDescription() { return CustomFactory.getDescription(this); }
	/** Returns the picture of this custom. */
	public String getPicture() { return CustomFactory.getPicture(this); }
	/** Returns the sound of this custom. */
	public String getSound() { return CustomFactory.getSound(this); }

	
	public void paint(Graphics g, int x, int y, float zf) {
		paint(g,x,y,zf,(double)0.0);
	}

	public void paint(Graphics g, int x, int y, float zf, double dlangle)
	{
    	g.setColor(CustomFactory.getColor(this));
    	double angle = dlangle - Math.toRadians(45.0);
    	int width = CustomFactory.getGraphicalWidth(this);
    	int length = CustomFactory.getGraphicalLength(this);
    
    	double corr = 1 - (width > length ? (double)length / (double)width : (double)width / (double)length);
    	
    	int[] cx = new int[4];
    	cx[0] = (int)(Math.round((double)x + Math.sin(angle) * width));
    	cx[1] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(90.0)) * width));
    	cx[2] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(180.0) + Math.atan(corr)) * length));
    	cx[3] = (int)(Math.round((double)x + Math.sin(angle + Math.toRadians(270.0) - Math.atan(corr)) * length));

    	int[] cy = new int[4];
    	cy[0] = (int)(Math.round((double)y + Math.cos(angle) * width));
    	cy[1] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(90.0)) * width));
    	cy[2] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(180.0) + Math.atan(corr)) * length));
    	cy[3] = (int)(Math.round((double)y + Math.cos(angle + Math.toRadians(270.0) - Math.atan(corr)) * length));

    g.fillPolygon(cx,cy,4);
  }

	
	// Specific XMLSerializable implementation 
 	public String getXMLName() {
 		return parentName+".roaduser-custom";
 	}
 	public void load (XMLElement myElement, XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
 	{
 		super.load(myElement, loader);
 		vehicle = myElement.getAttribute("vehicle").getIntValue();
 		driver = myElement.getAttribute("driver").getIntValue();
		passengers = (int[])XMLArray.loadArray(this, loader);
	}

 	public XMLElement saveSelf () throws XMLCannotSaveException
 	{
 		XMLElement result = super.saveSelf();
 		result.addAttribute(new XMLAttribute("vehicle", vehicle));
 		result.addAttribute(new XMLAttribute("driver", driver));
 		return result;
 	}

 	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
 	{
 		super.saveChilds(saver);
 		XMLArray.saveArray(passengers, this, saver,"passengers");
	}
}