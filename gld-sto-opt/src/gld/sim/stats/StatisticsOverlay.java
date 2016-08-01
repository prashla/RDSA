
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

package gld.sim.stats;

import gld.*;
import gld.sim.SimModel;
import gld.infra.*;
import gld.xml.*;

import java.awt.*;
import java.io.IOException;
import java.util.Observer;
import java.util.Observable;

/**
*
* Overlay for <code>gld.View</code>. Shows waiting queue lengths and relative 
* average junction waiting times.
*
* @author Group GUI
* @version 1.0
*/

public class StatisticsOverlay implements Overlay, Observer, XMLSerializable
{
	/** Current infrastructure. */
	Infrastructure infra;
	
	SpecialNode[] specialNodes;
	Junction[] junctions;
	float[] specialData, junctionData;
	float specialMax, junctionMax;
	int specialNum, junctionNum;
	String parentName="controller";
	
	/** Creates a default <code>StatisticsOverlay</code>. */
	public StatisticsOverlay(View _view, Infrastructure _infra)
	{	setInfrastructure(_infra);
	}
	
	/** Sets a new infrastructure as the current one and rereads all data. */
	public void setInfrastructure(Infrastructure _infra)
	{
		infra = _infra;
		specialNodes = infra.getSpecialNodes();
		specialNum = specialNodes.length;
		specialData = new float[specialNum];
		junctions = infra.getJunctions();
		junctionNum = junctions.length;
		junctionData = new float[junctionNum];
		
		specialMax = junctionMax = 0;
		
		for(int i=0; i<specialNum; i++)
			specialData[i] = 0;
		for(int i=0; i<junctionNum; i++)
			junctionData[i] = 0;
	}
	
	/** Updates the current data. */
	public void update(Observable obs, Object obj)
	{
		Infrastructure _infra = ((SimModel)obs).getInfrastructure();
		if(infra != _infra)
			setInfrastructure(_infra);
		refreshData();
	}
	
	/** Rereads the statistical data from the model. */
	public void refreshData()
	{
		specialMax = 0;
		junctionMax = 0.01f;
		
		for(int i=0; i<specialNum; i++)
			if((specialData[i] = specialNodes[i].getWaitingQueueLength()) > specialMax) 
				specialMax = specialData[i];
			
		for(int i=0; i<junctionNum; i++)
			if((junctionData[i] = junctions[i].getStatistics(0).getAvgWaitingTime(true)) > junctionMax) 
				junctionMax = junctionData[i];
	}
	
	public int overlayType() { return 2; }
	
	public void paint(Graphics g) throws GLDException
	{
		g.setPaintMode();
		
		// Waiting queues
		g.setColor(Color.blue);
		for(int i=0; i<specialNum; i++)
		{
			Rectangle r = specialNodes[i].getBounds();
			int dy = (int)(specialData[i] / 2);
			if(specialNodes[i].getRoadPos() == 3)
				r.x += r.width + 5;
			else
				r.x -= 10;
			r.y += r.height - dy + 1;
			r.width = 5;
			r.height = dy;
			g.fillRect(r.x, r.y, r.width, r.height);
		}
	
		// Average junction waiting times
		for(int i=0; i<junctionNum; i++)
		{
			Rectangle r = junctions[i].getBounds();
			g.setColor(new Color(0,0,0.7f,0.6f*junctionData[i]/junctionMax));
			g.fillRect(r.x, r.y, r.width, r.height);
		}
	}
	
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{ 	specialMax=myElement.getAttribute("special-max").getFloatValue();
		junctionMax=myElement.getAttribute("junction-max").getFloatValue();
		specialData=(float[])XMLArray.loadArray(this,loader);
		junctionData=(float[])XMLArray.loadArray(this,loader);	
	}

	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	XMLElement result=new XMLElement("overlaystats");
		result.addAttribute(new XMLAttribute("special-max",specialMax));
		result.addAttribute(new XMLAttribute("junction-max",junctionMax));
	  	return result;
	}
  
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{ 	XMLArray.saveArray(specialData,this,saver,"special-data");
		XMLArray.saveArray(junctionData,this,saver,"junction-data");
	}

	public String getXMLName ()
	{ 	return parentName+".overlaystats";
	}
	
	public void setParentName (String parentName)
	{	this.parentName=parentName; 
	}
}