
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

import gld.*;
import gld.utils.*;
import gld.xml.*;

import java.awt.*;
import java.awt.geom.*;
import java.io.IOException;
import java.util.*;

/**
 *
 * Basic Road.
 *
 * @author Group Datastructures
 * @version 1.0
 */

public class Road implements Selectable, XMLSerializable, TwoStageLoader
{
	/** Id of this Road. */
	protected int roadId;
	/** Length of this Road in blocks. */
	protected int length;
	/** Alpha Node of this Road. */
	protected Node alphaNode = null;
	/** Beta Node of this Road. */
	protected Node betaNode = null;
	/* For both lane[] goes that they fill from the right side of the road.
	   That is, the rightmost lane seen with the road is [0]. */
	/** Lanes where Roadusers move towards alphaNode. */
	protected Drivelane[] alphaLanes = { };
	/** Lanes where Roadusers move towards betaNode. */
	protected Drivelane[] betaLanes = { };
	/** Turns in the Road. */
	protected Turn[] turns = { };
	/** Container object to 1tranfer information from the first stage
	  * loader to the second stage loader.
	 */
	protected TwoStageLoaderData loadData=new TwoStageLoaderData();
	protected String parentName="model.infrastructure.node";

	/**
	 * Constructs a new Road
	 *
	 * @param _alpha The alpha node of this road
	 * @param _beta The beta node of this road
	 * @param _length The length of this road
	 */
	public Road(Node _alpha, Node _beta, int _length) {
		roadId = -1;
		alphaNode = _alpha;
		betaNode = _beta;
		length = _length;
	}

	/** Empty constructor for loading */
	public Road () { }








	/*============================================*/
	/* LOAD and SAVE                              */
	/*============================================*/

 	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
 	{ 	loadData.alphaNodeId=myElement.getAttribute("alpha-node").getIntValue();
   		loadData.betaNodeId=myElement.getAttribute("beta-node").getIntValue();
   		//System.out.println("Road "+roadId+" Loaded alphaNodeId:"+loadData.alphaNodeId+" and betaNodeId:"+loadData.betaNodeId);
   		roadId=myElement.getAttribute("id").getIntValue();
   		length=myElement.getAttribute("length").getIntValue();
		turns=(Turn[])XMLArray.loadArray(this,loader);
		loadData.alphaLaneIds=(int[])XMLArray.loadArray(this,loader);
		loadData.betaLaneIds=(int[])XMLArray.loadArray(this,loader);
 	}

 	public XMLElement saveSelf () throws XMLCannotSaveException
 	{ 	XMLElement result=new XMLElement("road");
   		result.addAttribute(new XMLAttribute("alpha-node",alphaNode.getId()));
   		result.addAttribute(new XMLAttribute("beta-node",betaNode.getId()));
   		result.addAttribute(new XMLAttribute("id",roadId));
   		result.addAttribute(new XMLAttribute("length",length));
   		result.addAttribute(new XMLAttribute("num-turns",turns.length));
   		return result;
 	}

 	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
 	{ 	XMLArray.saveArray (turns,this,saver,"turns");
		XMLArray.saveArray (getDrivelaneIdList(alphaLanes),this,saver,
					"alpha-lanes");
		XMLArray.saveArray (getDrivelaneIdList(betaLanes),this,saver,
		             		"beta-lanes");
	}


 	public String getXMLName ()
 	{ 	return parentName+".road";
 	}

	public void setParentName (String parentName)
	{	this.parentName=parentName;
	}

 	protected int[] getDrivelaneIdList ( Drivelane[] list)
 	{ 	int[] result=new int[list.length];
   		for (int t=0;t<list.length;t++)
   		{ 	if (list[t]==null)
     				result[t]=-1;
     			else
     				result[t]=list[t].getId();
   		}
   		return result;
 	}

 	class TwoStageLoaderData
 	{ 	int alphaNodeId,betaNodeId;
   		int[] alphaLaneIds,betaLaneIds;
 	}

 	public void loadSecondStage (Dictionary dictionaries) throws XMLInvalidInputException,XMLTreeException
 	{ 	// Load nodes
   		Dictionary nodeDictionary=(Dictionary)(dictionaries.get("node"));
   		alphaNode=(Node)(nodeDictionary.get(new Integer(loadData.alphaNodeId)));
   		betaNode=(Node)(nodeDictionary.get(new Integer(loadData.betaNodeId)));
/*   		System.out.println("Road:"+roadId+" alpaNodeId");
   		System.out.println(loadData.alphaNodeId+" ==> "+alphaNode);
   		System.out.println(alphaNode+" == alphaNode");
   		System.out.println("Road:"+roadId+" betaNodeId");
   		System.out.println(loadData.betaNodeId+" ==> "+betaNode);
   		System.out.println(betaNode+" == betaNode");
*/
   		// Load lanes
   		Dictionary laneDictionary=(Dictionary)(dictionaries.get("lane"));
   		alphaLanes=new Drivelane[loadData.alphaLaneIds.length];
   		betaLanes=new Drivelane[loadData.betaLaneIds.length];
   		for (int t=0;t<alphaLanes.length;t++)
      			alphaLanes[t]=(Drivelane)(laneDictionary.get
                    	(new Integer(loadData.alphaLaneIds[t])));
   		for (int t=0;t<betaLanes.length;t++)
      		betaLanes[t]=(Drivelane)(laneDictionary.get
                    	(new Integer(loadData.betaLaneIds[t])));
 	}






















	/*============================================*/
	/* Basic GET and SET methods                  */
	/*============================================*/


	/** Returns the Id of this Road */
	public int getId() { return roadId; }
	/** Sets the Id of this Road */
	public void setId(int id) { roadId = id; }

	/** Returns the alpha Node of this road */
	public Node getAlphaNode() { return alphaNode; }
	/** Sets the alpha Node of this road */
	public void setAlphaNode(Node n) { alphaNode = n; }

	/** Returns the beta Node of this road */
	public Node getBetaNode() { return betaNode; }
	/** Sets the beta Node of this road */
	public void setBetaNode(Node n) { betaNode = n; }

	/** Returns the length of this Road in blocks */
	public int getLength() { return length; }
	/** Sets the length of this Road in blocks */
	public void setLength(int l) { length = l; }

	/** Returns the Turns that this Road makes in the representation */
	public Turn[] getTurns() { return turns; }
	/** Sets the Turns that this Road makes in the representation */
	public void setTurns(Turn[] t) { turns = t; }

	/** Returns the Drivelanes where Roadusers move towards the alphaNode */
	public Drivelane[] getAlphaLanes() { return alphaLanes; }
	/** Sets the Drivelanes where Roadusers move towards the alphaNode */
	public void setAlphaLanes(Drivelane[] l) { alphaLanes = l; }

	/** Returns the Drivelanes where Roadusers move towards the betaNode */
	public Drivelane[] getBetaLanes() { return betaLanes; }
	/** Sets the Drivelanes where Roadusers move towards the betaNode */
	public void setBetaLanes(Drivelane[] l) { betaLanes = l; }

	/** Returns the name of this road. The name of a road is unique in an infrastructure. */
	public String getName() { return "Road " + roadId; }














	/*============================================*/
	/* Selectable                                 */
	/*============================================*/


	public Rectangle getBounds() { return getComplexBounds().getBounds(); }
	public Shape getComplexBounds() {
		Area a = new Area();
		for (int i=0; i < alphaLanes.length; i++)
			a.add(new Area(alphaLanes[i].getComplexBounds()));
		for (int i=0; i < betaLanes.length; i++)
			a.add(new Area(betaLanes[i].getComplexBounds()));
		return a;
	}

	public int getDistance(Point p) { return (int)getCenterPoint().distance(p); }
	public Point getSelectionPoint() { return getCenterPoint(); }
	public Point[] getSelectionPoints() { return null; }
	public Point getCenterPoint() {
		Rectangle r = getBounds();
		return new Point(r.x + r.width / 2, r.y + r.height / 2);
	}

	public boolean hasChildren() { return true; }
	public boolean isSelectable() { return true; }
	public Enumeration getChildren() { return new ArrayEnumeration(alphaLanes, betaLanes); }


















	/*============================================*/
	/* MODIFYING DATA                             */
	/*============================================*/


	/**
	 * Resets the road.
	 * This will reset all Drivelanes.
	 * @see Drivelane#reset()
	 */
	public void reset() {
		for (int i=0; i < alphaLanes.length; i++) {
			alphaLanes[i].reset();
		}
		for (int i=0; i < betaLanes.length; i++) {
			betaLanes[i].reset();
		}
	}

	/**
	 * Adds a Drivelane leading to the alpha node
	 *
	 * @param l The drivelane to add
	 * @throw InfraException if l is null
	 */
	public void addAlphaLane(Drivelane l) throws InfraException {
		if (l == null) throw new InfraException("Parameter l is null");
		alphaLanes = (Drivelane[])Arrayutils.addElement(alphaLanes, l);
	}

	/**
	 * Removes a Drivelane leading to the alpha node
	 *
	 * @param l The drivelane to remove
	 * @throw InfraException if l is null
	 * @throw InfraException if the lane is not an alpha lane on this road
	 */
	public void remAlphaLane(Drivelane l) throws InfraException {
		if (l == null) throw new InfraException("Parameter l is null");
		int i = Arrayutils.findElement(alphaLanes, l);
		if (i == -1) throw new InfraException("Lane is not an alpha lane on this road");
		alphaLanes = (Drivelane[])Arrayutils.remElement(alphaLanes, i);
	}

	/**
	 * Adds a Drivelane leading to the beta node
	 *
	 * @param l The drivelane to add
	 * @throw InfraException if l is null
	 */
	public void addBetaLane(Drivelane l) throws InfraException {
		if (l == null) throw new InfraException("Parameter l is null");
		betaLanes = (Drivelane[])Arrayutils.addElement(betaLanes, l);
	}


	/**
	 * Removes a Drivelane leading to the beta node
	 *
	 * @param l The drivelane to remove
	 * @throw InfraException if l is null
	 * @throw InfraException if the lane is not a beta lane on this road
	 */
	public void remBetaLane(Drivelane l) throws InfraException {
		if (l == null) throw new InfraException("Parameter l is null");
		int i = Arrayutils.findElement(betaLanes, l);
		if (i == -1) throw new InfraException("Lane is not an beta lane on this road");
		alphaLanes = (Drivelane[])Arrayutils.remElement(betaLanes, i);
	}


	/**
	 * Adds a Turn to this road
	 *
	 * @param t The turn to add
	 * @throw InfraException if t is null
	 */
	public void addTurn(Turn t) throws InfraException {
		if (t == null) throw new InfraException("Parameter t is null");
		turns = (Turn[])Arrayutils.addElement(turns, t);
	}

	/**
	 * Removes a Turn from this road
	 *
	 * @param t The turn to remove
	 * @throw InfraException if t is null
	 * @throw InfraException if the turn is not part of this road
	 */
	public void remTurn(Turn t) throws InfraException {
		if (t == null) throw new InfraException("Parameter t is null");
		int i = Arrayutils.findElement(turns, t);
		if (i == -1) throw new InfraException("Turn is not part of this road");
		turns = (Turn[])Arrayutils.remElement(turns, i);
	}


	/**
	 * Adds a Drivelane to this road
	 *
	 * @param lane The drivelane to add
	 * @param to The Node this Drivelane leads to
	 * @throw InfraException if lane or to is null
	 * @throw InfraException if the road is not connected to the given node
	 */
	public void addLane(Drivelane lane, Node to) throws InfraException {
		if (lane == null) throw new InfraException("Parameter lane is null");
		if (to == null) throw new InfraException("Parameter to is null");
		if (to == alphaNode) addAlphaLane(lane);
		else if (to == betaNode) addBetaLane(lane);
		else throw new InfraException("Road is not connected to this node");
	}


	/**
	 * Removes a Drivelane from this road
	 *
	 * @param lane The drivelane to remove
	 * @throw InfraException if lane is null
	 * @throw InfraException if lane is neither an alpha, nor a beta lane on this road
	 */
	public void remLane(Drivelane lane) throws InfraException {
		if (lane == null) throw new InfraException("Parameter lane is null");
		int i = Arrayutils.findElement(alphaLanes, lane);
		if (i != -1) {
			alphaLanes = (Drivelane[])Arrayutils.remElement(alphaLanes, i);
			return;
		}
		i = Arrayutils.findElement(betaLanes, lane);
		if (i != -1) {
			betaLanes = (Drivelane[])Arrayutils.remElement(betaLanes, i);
			return;
		}
		throw new InfraException("Lane is neither an alpha, nor a beta lane on this road");
	}















	/*============================================*/
	/* COMPLEX GET                                */
	/*============================================*/


	/** Returns the other node this road is connected to, seen from the given node */
	public Node getOtherNode(Node seenfrom) { return seenfrom == alphaNode ? betaNode : alphaNode; }


	/** Returns the width of this road in number of lanes */
	public int getWidth() { return getNumAllLanes(); }

	/** Returns the number of total Drivelanes in this road */
	public int getNumAllLanes() { return alphaLanes.length + betaLanes.length; }

	public int getNumAlphaLanes() { return alphaLanes.length; }
	public int getNumBetaLanes() { return betaLanes.length; }

	/**
	 * Returns the number of Drivelanes that are outbound from the given node
	 *
	 * @throw InfraException if n is null
	 * @throw InfraException if this road is not connected to the given node
	 */
	public int getNumOutboundLanes(Node n) throws InfraException {
		return getOutboundLanes(n).length;
	}

	/**
	 * Returns the number of Drivelanes that are inbound on the given node
	 *
	 * @throw InfraException if n is null
	 * @throw InfraException if this road is not connected to the given node
	 */
	public int getNumInboundLanes(Node n) throws InfraException {
		return getInboundLanes(n).length;
	}

	/**
	 * Returns an array of Drivelanes which are outbound from the given node.
	 *
	 * @throw InfraException if n is null
	 * @throw InfraException if this road is not connected to the given node
	 */
	public Drivelane[] getOutboundLanes(Node n) throws InfraException
	{
		if (n == null) throw new InfraException("Parameter n is null");
		else if (n == alphaNode) return betaLanes;
		else if (n == betaNode) return alphaLanes;
		throw new InfraException("Node is not connected to this road");
	}

	/**
	 * Returns an array of Drivelanes which are inbound on the node given.
	 *
	 * @throw InfraException if n is null
	 * @throw InfraException if this road is not connected to the given node
	 */
	public Drivelane[] getInboundLanes(Node n) throws InfraException
	{
		//System.out.println("Requested Node.getId(): "+n.getId());
		//System.out.println("Alpha Node.getId(): "+((alphaNode==null)?"null!":""+alphaNode.getId()));
		//System.out.println("Beta Node.getId(): "+betaNode.getId());
		if (n == null) throw new InfraException("Parameter n is null");
		else if (n == alphaNode) return alphaLanes;
		else if (n == betaNode) return betaLanes;
		throw new InfraException("Node is not connected to this road");
	}

	/** Returns an array of Drivelanes which are part of this road */
	public Drivelane[] getAllLanes() {
		return (Drivelane[])Arrayutils.addArray(alphaLanes, betaLanes);
	}


 /**
   * Draws this Road, with all it's DriveLanes and Roadusers onto a Graphics object, mostly
   * a View
   * @param g The Graphics object of the View
   * @param x The horizontal scroll number of the View
   * @param y The vertical scroll number of the View
   * @param zf The zoom factor of the View
   */
  public void paint(Graphics g, int x, int y, float zf) throws GLDException
  {	this.paint(g,x,y,zf,0.0);
  }

  public void paint(Graphics g) throws GLDException
  {   paint(g, 0, 0, 1.0f, 0.0);
  }


  /**
   * Paints an arbitrary segment of a drivelane and returns the outlines of this segment
   * as a Shape
   */
  public Shape paintDrivelaneSegment(Graphics g, Point p1_1, Point p1_2, Point p2_1, Point p2_2, double a1, double a2, Drivelane dl)
  {
  	switch(dl.getType())
    {
        case 1: g.setColor(new Color(230,230,255)); break;      // Cars only
        case 2: g.setColor(new Color(255,230,230)); break;      // Busses only
        case 3: g.setColor(new Color(255,230,255)); break;      // Cars and Busses only (default, geloof ik)
        case 4: g.setColor(new Color(230,255,230)); break;      // Bicycles only
        case 5: g.setColor(new Color(230,255,255)); break;      // Cars and Bicycles only
        case 6: g.setColor(new Color(255,255,230)); break;      // Busses and Bicycles only
        default: g.setColor(Color.white); break;                 // All Roadusers
    }
    try
    {
        TurnCurve arc1 = CurveUtils.createCurve(p1_1,p1_2,a1,a2);
        TurnCurve arc2 = CurveUtils.createCurve(p2_1,p2_2,a1,a2);
        GeneralPath path = CurveUtils.createPath(arc1,arc2);
        ((Graphics2D)g).fill(path);
        g.setColor(Color.black);
        ((Graphics2D)g).draw(arc1);
        ((Graphics2D)g).draw(arc2);
		return path;
    }
    catch(CurveException exc)
    {
		System.out.println("CurveException: " + exc.getMessage());
    }
	return null;
  }


	/**
	 * Draws a Roaduser on a (segment of a) Drivelane
	 * @param g The Graphics object to draw onto
	 * @param r The Roaduser to draw
	 * @param p1 The first point of the segment
	 * @param p2 The last point of the segment
	 * @param a1 The angle of the first point of the segment
	 * @param a2 The angle of the last point of the segment
	 * @param pos The position of the Roaduser relative to the start of the segment
	 */
 	public void paintRoaduser(Graphics g, Roaduser r, Point p1, Point p2, double a1, double a2, int pos)
  	{
  		try
  		{
  			TurnCurve arc = CurveUtils.createCurve(p1,p2,a1,a2);
//			Object[] o = CurveUtils.getPointAndAngle(arc,pos,10);
//			Point  p = (Point)o[0];
//			double a = ((Double)o[1]).doubleValue();

			Point p = CurveUtils.getPoint(arc,pos,10);
			double a = CurveUtils.getAngle(arc,pos,10);

			r.paint(g,p.x,p.y,(float)1.0,a);
 		}
  		catch(Exception x)
  		{
  			System.out.println("Exception: " + x.getMessage());
  		}
  	}

	/**
	 * Draws a Roaduser on a Node
	 * @param g The Graphics object to draw onto
	 * @param ru The Roaduser to draw
	 * @param node The Node this Roaduser is on
	 * @param p The first point of the new Drivelane the Roaduser will be on
	 * @param pos The position of the Roaduser on the Node
	 * @param d The new Drivelane
	 */
	public void paintRoaduserOnNode(Graphics g, Roaduser ru, Node node, Point p, int pos, Drivelane d)
	{
		double na = 0.0;
		try
		{
			int id = ru.getPrevSign();
			Drivelane[] dls = node.getLanesLeadingTo(d,0);
			for(int i = 0; i < dls.length; i++)
				if(dls[i].getId() == id)
				{
					Road road = dls[i].getRoad();
					int c = node.isConnectedAt(road);
					switch(c)
					{
						case(0): na = Math.PI * 0.5;	break;
						case(1): na = 0.0;		break;
						case(2): na = -Math.PI * 0.5;	break;
						case(3): na = Math.PI;		break;
					}
					dls = road.getInboundLanes(node);
					for(int j = 0; j < dls.length; j++)
						if(dls[j].getId() == id)
						{
							Point q = new Point((int)((double)node.getCoord().x + Math.cos(na) * 5.0 * node.getWidth() - Math.sin(na) * 10.0 * (double)j - Math.sin(na) * 5.0),
                                                (int)((double)node.getCoord().y - Math.sin(na) * 5.0 * node.getWidth() - Math.cos(na) * 10.0 * (double)j - Math.cos(na) * 5.0));

                            double ra;
                            if(q.x == p.x)
                                ra = 0.0;
                            else if(q.y == p.y)
                                ra = Math.PI * 0.5;
                            else
                                ra = Math.atan((double)(q.x - p.x) / (double)(q.y - p.y));

                            if(pos < 2 && pos >= -2)
                            	q = avg3(p,p,q);
                            else
                            	q = avg3(p,q,q);

                            ru.paint(g,q.x,q.y,1.0f,ra);
                            return;
                        }
                }
      	}
  		catch(Exception x)
		{}
	}


	public void paintDrivelaneAttributes(Graphics g, Point p, double a, Drivelane d)
	{
// draw sign
		Point q = new Point((int)(p.x - 7 * Math.cos(a)),(int)(p.y + 7 * Math.sin(a)));
		Sign s = d.getSign();
		if(s != null && s.getType() == Sign.TRAFFICLIGHT)
		{
			g.setColor(s.mayDrive() ? Color.green : Color.red);
			g.fillOval(q.x-2,q.y-2,6,6);
		}

		// return if lane leads to SpecialNodes
		if(d.getNodeLeadsTo() instanceof SpecialNode) return;

// draw directions
		Point q1 = new Point((int)(p.x + 5 * Math.cos(a)),(int)(p.y - 5 * Math.sin(a)));
		Point q2 = new Point((int)(p.x + 8 * Math.cos(a)),(int)(p.y - 8 * Math.sin(a)));
		Point q3 = new Point((int)(p.x + 3 * Math.cos(a)),(int)(p.y - 3 * Math.sin(a)));
		Point q4 = new Point((int)(p.x + Math.cos(a)),(int)(p.y - Math.sin(a)));
		g.setColor(Color.black);
		g.drawLine(q1.x,q1.y,q2.x,q2.y);
		
		boolean[] tars = d.getTargets();
		if(tars[0])	// left
		{
			switch((int)Math.toDegrees(CurveUtils.normalize(a)))
			{
				case(90):  g.drawLine(q1.x,q1.y,q1.x + 3,q3.y); break;
				case(0):   g.drawLine(q1.x,q1.y,q3.x,q1.y + 3); break;
				case(360): g.drawLine(q1.x,q1.y,q3.x,q1.y + 3); break;
				case(270): g.drawLine(q1.x,q1.y,q1.x - 3,q3.y); break;
				case(180): g.drawLine(q1.x,q1.y,q3.x,q1.y - 3); break;
			}
		}
		if(tars[1]) // straight
			g.drawLine(q1.x,q1.y,q4.x,q4.y);
		if(tars[2]) // right
		{
			switch((int)Math.toDegrees(CurveUtils.normalize(a)))
			{
				case(90):  g.drawLine(q1.x,q1.y,q1.x - 3,q3.y); break;
				case(0):   g.drawLine(q1.x,q1.y,q3.x,q1.y - 3); break;
				case(360): g.drawLine(q1.x,q1.y,q3.x,q1.y - 3); break;
				case(270): g.drawLine(q1.x,q1.y,q1.x + 3,q3.y); break;
				case(180): g.drawLine(q1.x,q1.y,q3.x,q1.y + 3); break;
			}
		}
	}		


  private Point avg(Point a, Point b) {
  	return new Point(a.x + ((b.x - a.x) / 2), a.y + ((b.y - a.y) / 2));
  }
  
  private Point avg3(Point a, Point b, Point c) {
  	Point p = avg(a,b);
  	Point q = avg(b,c);
  	return avg(p,q);
  }


  private Point getCurvePoint(Point p, double a, int d, int w, int x, int y)
  {
    return new Point((int)(x + p.x + d * Math.sin(a) * w * 10),
                     (int)(y + p.y + d * Math.cos(a) * w * 10));
  }

  /**
   * Paint this road with it's roadusers
   */
  public void paint(Graphics g, int x, int y, float zf, double bogus) throws GLDException
  {
      int sx1 = 0, sx2 = 0, ex1 = 0, ex2 = 0;
      int sy1 = 0, sy2 = 0, ey1 = 0, ey2 = 0;
      int alphaMaxWidth = alphaNode.getWidth();
      int betaMaxWidth = betaNode.getWidth();
      int width = getWidth();

      int cax, cay, cbx, cby;
      double alphaAngle = 0.0, betaAngle = 0.0;

      int acp, bcp;
      acp = alphaNode.isConnectedAt(this);
      bcp = betaNode.isConnectedAt(this);

// calculate the starting positions of the road
      if(acp == 0)    // alpha top
      {
          cax = alphaNode.coord.x;
          cay = alphaNode.coord.y - alphaMaxWidth * 5 - 1;
          alphaAngle = Math.toRadians(90); // Math.PI * 0.5
      }
      else if(acp == 1) // alpha right
      {
          cax = alphaNode.coord.x + alphaMaxWidth * 5 + 1;
          cay = alphaNode.coord.y;
          alphaAngle = 0.0;
      }
      else if(acp == 2) // alpha bottom
      {
          cax = alphaNode.coord.x;
          cay = alphaNode.coord.y + alphaMaxWidth * 5 + 1;
          alphaAngle = Math.toRadians(270); // Math.PI * 1.5
      }
      else    // alpha left
      {
          cax = alphaNode.coord.x - alphaMaxWidth * 5 - 1;
          cay = alphaNode.coord.y;
          alphaAngle = Math.toRadians(180); // Math.PI
      }

      if(bcp == 0)    // beta top
      {
          cbx = betaNode.coord.x;
          cby = betaNode.coord.y - betaMaxWidth * 5 - 1;
          betaAngle = Math.toRadians(270);
      }
      else if(bcp == 1) // beta right
      {
          cbx = betaNode.coord.x + betaMaxWidth * 5 + 1;
          cby = betaNode.coord.y;
          betaAngle = Math.toRadians(180);
      }
      else if(bcp == 2) // beta bottom
      {
          cbx = betaNode.coord.x;
          cby = betaNode.coord.y + betaMaxWidth * 5 + 1;
          betaAngle = Math.toRadians(90);
      }
      else    // beta left
      {
          cbx = betaNode.coord.x - betaMaxWidth * 5 - 1;
          cby = betaNode.coord.y;
          betaAngle = 0.0;
      }

      Point alphaCoord = new Point(cax,cay);
      Point betaCoord = new Point(cbx,cby);

	  boolean vertical = false;
	  if(Math.abs(alphaCoord.x - betaCoord.x) < 5)
	  {
	  	if(turns.length > 0)
	  	{
	  		if(Math.abs(alphaCoord.x - turns[0].getCoord().x) < 5)
	  			vertical = true;
	  	}
	  	else
	  		vertical = true;
	  }


      if(turns.length > 0 && !vertical)
      {
          double ta = CurveUtils.normalize(turns[0].getAngle() + Math.PI);
          for(int n = 0; n < alphaLanes.length; n++)
          {
  			  Point p1_1 = getCurvePoint(alphaCoord,alphaAngle,-1,n,x,y);
			  Point p1_2 = getCurvePoint(turns[0].getCoord(),ta,1,n,x,y);
			  Point p2_1 = getCurvePoint(alphaCoord,alphaAngle,-1,n+1,x,y);
			  Point p2_2 = getCurvePoint(turns[0].getCoord(),ta,1,n+1,x,y);
			  Shape bounds = paintDrivelaneSegment(g,p1_1,p1_2,p2_1,p2_2,alphaAngle,ta,alphaLanes[n]);
			  alphaLanes[n].clearCurveBounds();
			  alphaLanes[n].addCurveBounds(bounds);
			  paintDrivelaneAttributes(g,avg(p1_1,p2_1),alphaAngle,alphaLanes[n]);
          }
          for(int n = 0; n < betaLanes.length; n++)
          {
              Point p1_1 = getCurvePoint(alphaCoord,alphaAngle,1,n,x,y);
              Point p1_2 = getCurvePoint(turns[0].getCoord(),ta,-1,n,x,y);
              Point p2_1 = getCurvePoint(alphaCoord,alphaAngle,1,n+1,x,y);
              Point p2_2 = getCurvePoint(turns[0].getCoord(),ta,-1,n+1,x,y);
			  Shape bounds = paintDrivelaneSegment(g,p1_1,p1_2,p2_1,p2_2,alphaAngle,ta,betaLanes[n]);
			  betaLanes[n].clearCurveBounds();
			  betaLanes[n].addCurveBounds(bounds);
		  }
      }

      for(int t = 0; t < turns.length - 1; t++)
      {
          double ta1 = CurveUtils.normalize(turns[t].getAngle());
          double ta2 = CurveUtils.normalize(turns[t+1].getAngle() + Math.PI);
          Point tp1 = turns[t].getCoord();
          Point tp2 = turns[t+1].getCoord();
          for(int n = 0; n < alphaLanes.length; n++)
          {
			  Point p1_1 = getCurvePoint(tp1,ta1,-1,n,x,y);
			  Point p1_2 = getCurvePoint(tp2,ta2,1,n,x,y);
			  Point p2_1 = getCurvePoint(tp1,ta1,-1,n+1,x,y);
			  Point p2_2 = getCurvePoint(tp2,ta2,1,n+1,x,y);
			  Shape bounds = paintDrivelaneSegment(g,p1_1,p1_2,p2_1,p2_2,ta1,ta2,alphaLanes[n]);
			  alphaLanes[n].addCurveBounds(bounds);
          }
          for(int n = 0; n < betaLanes.length; n++)
          {
			  Point p1_1 = getCurvePoint(tp1,ta1,1,n,x,y);
			  Point p1_2 = getCurvePoint(tp2,ta2,-1,n,x,y);
			  Point p2_1 = getCurvePoint(tp1,ta1,1,n+1,x,y);
			  Point p2_2 = getCurvePoint(tp2,ta2,-1,n+1,x,y);
              Shape bounds = paintDrivelaneSegment(g,p1_1,p1_2,p2_1,p2_2,ta1,ta2,betaLanes[n]);
              betaLanes[n].addCurveBounds(bounds);
          }
      }

      if(turns.length > 0 && !vertical)
      {
          int last = turns.length - 1;
          double ta = CurveUtils.normalize(turns[last].getAngle());
          Point tp = turns[last].getCoord();
          for(int n = 0; n < betaLanes.length; n++)
          {
			  Point p1_1 = getCurvePoint(tp,ta,1,n,x,y);
			  Point p1_2 = getCurvePoint(betaCoord,betaAngle,1,n,x,y);
			  Point p2_1 = getCurvePoint(tp,ta,1,n+1,x,y);
			  Point p2_2 = getCurvePoint(betaCoord,betaAngle,1,n+1,x,y);
              Shape bounds = paintDrivelaneSegment(g,p1_1,p1_2,p2_1,p2_2,ta,betaAngle,betaLanes[n]);
			  paintDrivelaneAttributes(g,avg(p1_2,p2_2),betaAngle + Math.PI,betaLanes[n]);
			  betaLanes[n].addCurveBounds(bounds);
          }
          for(int n = 0; n < alphaLanes.length; n++)
          {
			  Point p1_1 = getCurvePoint(tp,ta,-1,n,x,y);
			  Point p1_2 = getCurvePoint(betaCoord,betaAngle,-1,n,x,y);
			  Point p2_1 = getCurvePoint(tp,ta,-1,n+1,x,y);
			  Point p2_2 = getCurvePoint(betaCoord,betaAngle,-1,n+1,x,y);
              Shape bounds = paintDrivelaneSegment(g,p1_1,p1_2,p2_1,p2_2,ta,betaAngle,alphaLanes[n]);
              alphaLanes[n].addCurveBounds(bounds);
          }




// draw roadusers:
		for(int n = 0; n < alphaLanes.length; n++)
		{
			Drivelane dl = alphaLanes[n];
			LinkedList queue = dl.getQueue();
			if(queue.size() > 0)
			{
				ListIterator li = queue.listIterator();
				while(li.hasNext())
				{
					Roaduser r = (Roaduser)li.next();
					int pos = r.getPosition();
					
					alphaturnloop:
					{
						for(int t = 0; t <= turns.length; t++)
						{
							if(t < turns.length && pos < length)
							{
								if(turns[t].getPosition() >= pos)
								{
									Point tp1, tp2;
									double ta1, ta2;
									int cpos = pos;
									if(t==0)	// alpha->turns[0]
									{
										tp1 = new Point((int)(alphaCoord.x + 5 * Math.cos(alphaAngle)),
										                (int)(alphaCoord.y - 5 * Math.sin(alphaAngle)));
										tp2 = turns[t].getCoord();
										ta1 = alphaAngle;
										ta2 = turns[t].getAngle() + Math.PI;
									}
									else
									{
										tp1 = turns[t-1].getCoord();
										tp2 = turns[t].getCoord();
										ta1 = turns[t-1].getAngle();
										ta2 = turns[t].getAngle() + Math.PI;
										cpos = pos - (int)(turns[t-1].getPosition());
									}
									Point p1 = avg(getCurvePoint(tp1,ta1,-1,n,x,y),
								    	           getCurvePoint(tp1,ta1,-1,(n+1),x,y));
									Point p2 = avg(getCurvePoint(tp2,ta2,1,n,x,y),
								    	           getCurvePoint(tp2,ta2,1,(n+1),x,y));
									paintRoaduser(g,r,p1,p2,ta1,ta2,cpos);
									break alphaturnloop;
								}
							}
							else if(pos < length)
							{
								Point tp1 = turns[t-1].getCoord();
								Point tp2 = new Point((int)(betaCoord.x - 5 * Math.cos(betaAngle)),
								                      (int)(betaCoord.y + 5 * Math.sin(betaAngle)));
								double ta1 = turns[t-1].getAngle();
								double ta2 = betaAngle;
								int cpos = pos - (int)(turns[t-1].getPosition()) + 1;
								Point p1 = avg(getCurvePoint(tp1,ta1,-1,n,x,y),
											   getCurvePoint(tp1,ta1,-1,n+1,x,y));
								Point p2 = avg(getCurvePoint(tp2,ta2,-1,n,x,y),
								               getCurvePoint(tp2,ta2,-1,n+1,x,y));
								paintRoaduser(g,r,p1,p2,ta1,ta2,cpos);
							}
							else
							{
								Point tp1 = new Point((int)(betaCoord.x - 5 * Math.cos(betaAngle)),
								                      (int)(betaCoord.y + 5 * Math.sin(betaAngle)));
								Point p = avg(getCurvePoint(tp1,betaAngle,-1,n,x,y),
								              getCurvePoint(tp1,betaAngle,-1,n+1,x,y));
                                paintRoaduserOnNode(g,r,betaNode,p,pos-length,dl);
							}
								
						}
					} // alphaturnloop
				}		
			}
		}
		
		for(int n = 0; n < betaLanes.length; n++)
		{
			Drivelane dl = betaLanes[n];
			LinkedList queue = dl.getQueue();
			if(queue.size() > 0)
			{
				ListIterator li = queue.listIterator();
				while(li.hasNext())
				{
					Roaduser r = (Roaduser)li.next();
					int pos = length - r.getPosition() - 2;
					betaturnloop:
					{
						for(int t = 0; t <= turns.length; t++)
						{
                            if(t < turns.length && pos >= 0)
							{
								if(turns[t].getPosition() >= pos)
								{
									Point tp1, tp2;
									double ta1, ta2;
									int cpos = pos;
									if(t==0)	// alpha->turns[0]
									{
										tp1 = new Point((int)(alphaCoord.x + 5 * Math.cos(alphaAngle)),
										                (int)(alphaCoord.y - 5 * Math.sin(alphaAngle)));
										tp2 = turns[t].getCoord();
										ta1 = alphaAngle;
										ta2 = turns[t].getAngle() + Math.PI;
									}
									else
									{
										tp1 = turns[t-1].getCoord();
										tp2 = turns[t].getCoord();
										ta1 = turns[t-1].getAngle();
										ta2 = turns[t].getAngle() + Math.PI;
										cpos = pos - (int)(turns[t-1].getPosition());
									}
									Point p1 = avg(getCurvePoint(tp1,ta1,1,n,x,y),
								    	           getCurvePoint(tp1,ta1,1,(n+1),x,y));
									Point p2 = avg(getCurvePoint(tp2,ta2,-1,n,x,y),
								    	           getCurvePoint(tp2,ta2,-1,(n+1),x,y));
									paintRoaduser(g,r,p1,p2,ta1,ta2,cpos);
									break betaturnloop;
								}
							}
                            else if(pos >= 0)
							{
								Point tp1 = turns[t-1].getCoord();
								Point tp2 = new Point((int)(betaCoord.x - 5 * Math.cos(betaAngle)),
								                      (int)(betaCoord.y + 5 * Math.sin(betaAngle)));
								double ta1 = turns[t-1].getAngle();
								double ta2 = betaAngle;
								int cpos = pos - (int)(turns[t-1].getPosition()) + 1;
								Point p1 = avg(getCurvePoint(tp1,ta1,1,n,x,y),
											   getCurvePoint(tp1,ta1,1,n+1,x,y));
								Point p2 = avg(getCurvePoint(tp2,ta2,1,n,x,y),
								               getCurvePoint(tp2,ta2,1,n+1,x,y));
								paintRoaduser(g,r,p1,p2,ta1,ta2,cpos);
							}
							else
							{
                              	Point tp1 = new Point((int)(alphaCoord.x - 5 * Math.cos(alphaAngle)),
					            					  (int)(alphaCoord.y + 5 * Math.sin(alphaAngle)));
								Point p = avg(getCurvePoint(tp1,alphaAngle,1,n,x,y),
								              getCurvePoint(tp1,alphaAngle,1,n+1,x,y));
                                paintRoaduserOnNode(g,r,alphaNode,p,pos,dl);
							}

						}
					} // betaturnloop
				}		
			}
		}


	 }
	  else      // straight road
      {

// ======== OLD ROAD CODE ==================================================
 		if(acp == 0)
		{
			sx1 = alphaNode.coord.x - (int)(width * 5);
			sy1 = alphaNode.coord.y - (int)(alphaMaxWidth * 5) - 1;
			sx2 = alphaNode.coord.x + (int)(width * 5);
			sy2 = alphaNode.coord.y - (int)(alphaMaxWidth * 5) - 1;
		}
		else if(acp == 1)
		{
			sx1 = alphaNode.coord.x + (int)(alphaMaxWidth * 5) + 1;
			sy1 = alphaNode.coord.y - (int)(width * 5);
			sx2 = alphaNode.coord.x + (int)(alphaMaxWidth * 5) + 1;
			sy2 = alphaNode.coord.y + (int)(width * 5);
		}
		else if(acp == 2)
		{
			sx1 = alphaNode.coord.x + (int)(width * 5);
			sy1 = alphaNode.coord.y + (int)(alphaMaxWidth * 5) + 1;
			sx2 = alphaNode.coord.x - (int)(width * 5);
			sy2 = alphaNode.coord.y + (int)(alphaMaxWidth * 5) + 1;
		}
		else
		{
			sx1 = alphaNode.coord.x - (int)(alphaMaxWidth * 5) - 1;
			sy1 = alphaNode.coord.y + (int)(width * 5);
			sx2 = alphaNode.coord.x - (int)(alphaMaxWidth * 5) - 1;
			sy2 = alphaNode.coord.y - (int)(width * 5);
		}

		if(bcp == 0)
		{
			ex2 = betaNode.coord.x - (int)(width * 5);
			ey2 = betaNode.coord.y - (int)(betaMaxWidth * 5) - 1;
			ex1 = betaNode.coord.x + (int)(width * 5);
			ey1 = betaNode.coord.y - (int)(betaMaxWidth * 5) - 1;
		}
		else if(bcp == 1)
		{
			ex2 = betaNode.coord.x + (int)(betaMaxWidth * 5) + 1;
			ey2 = betaNode.coord.y - (int)(width * 5);
			ex1 = betaNode.coord.x + (int)(betaMaxWidth * 5) + 1;
			ey1 = betaNode.coord.y + (int)(width * 5);
		}
		else if(bcp == 2)
		{
			ex2 = betaNode.coord.x + (int)(width * 5);
			ey2 = betaNode.coord.y + (int)(betaMaxWidth * 5) + 1;
			ex1 = betaNode.coord.x - (int)(width * 5);
			ey1 = betaNode.coord.y + (int)(betaMaxWidth * 5) + 1;
		}
		else
		{
			ex2 = betaNode.coord.x - (int)(betaMaxWidth * 5) - 1;
			ey2 = betaNode.coord.y + (int)(width * 5);
			ex1 = betaNode.coord.x - (int)(betaMaxWidth * 5) - 1;
			ey1 = betaNode.coord.y - (int)(width * 5);
		}


    // calculate angle between nodes
	double dlangle;
	if(alphaNode.coord.y == betaNode.coord.y)
	{
		if(alphaNode.coord.x > betaNode.coord.y)
			dlangle = Math.toRadians(-90.0);
		else
			dlangle = Math.toRadians(90.0);
	}
	else if(betaNode.coord.y < alphaNode.coord.y && betaNode.coord.x > alphaNode.coord.x)
		dlangle = Math.toRadians(180) + Math.atan((double)(betaNode.coord.x - alphaNode.coord.x) / (double)(betaNode.coord.y - alphaNode.coord.y));
	else
		dlangle = Math.atan((double)(betaNode.coord.x - alphaNode.coord.x) / (double)(betaNode.coord.y - alphaNode.coord.y));

	Drivelane[] dls = getAllLanes();
	for(int n = 1; n < width + 1; n++)
	{
		float part = (float) n/width;
		float part2 = (float) (n-1) / width;

		switch(dls[width - n].getType())
		{
			case 1: g.setColor(new Color(230,230,255)); break;	// Cars only
			case 2: g.setColor(new Color(255,230,230)); break;	// Busses only
			case 3: g.setColor(new Color(255,230,255)); break;	// Cars and Busses only (default, geloof ik)
			case 4: g.setColor(new Color(230,255,230)); break;	// Bicycles only
			case 5: g.setColor(new Color(230,255,255)); break;	// Cars and Bicycles only
			case 6: g.setColor(new Color(255,255,230)); break;	// Busses and Bicycles only
			default: g.setColor(Color.white); break;			// All Roadusers
		}

		int[] cx = new int[4];
		int[] cy = new int[4];

		cx[0] = (int)((x + sx1 * part + sx2 * (1-part)) * zf);
		cx[1] = (int)((x + sx1 * part2 + sx2 * (1-part2)) * zf);
		cx[2] = (int)((x + ex1 * part2 + ex2 * (1-part2)) * zf);
		cx[3] = (int)((x + ex1 * part + ex2 * (1-part)) * zf);

		cy[0] = (int)((y + sy1 * part + sy2 * (1-part)) * zf);
		cy[1] = (int)((y + sy1 * part2 + sy2 * (1-part2)) * zf);
		cy[2] = (int)((y + ey1 * part2 + ey2 * (1-part2)) * zf);
		cy[3] = (int)((y + ey1 * part + ey2 * (1-part)) * zf);

		g.fillPolygon(cx,cy,4);
  		dls[width-n].clearCurveBounds();
		dls[width-n].addCurveBounds(new Polygon(cx,cy,4));
	}

    g.setColor(Color.black);
	g.drawLine((int)((sx1 + x) * zf),(int)((sy1 + y) * zf),(int)((ex1 + x) * zf),(int)((ey1 + y) * zf));
	g.drawLine((int)((sx2 + x) * zf),(int)((sy2 + y) * zf),(int)((ex2 + x) * zf),(int)((ey2 + y) * zf));

	for(int n = 1; n < width; n++)
	{
		float part = (float) n/width;
		g.drawLine((int)((x + sx1 * part + sx2 * (1-part)) * zf),
						     (int)((y + sy1 * part + sy2 * (1-part)) * zf),
						     (int)((x + ex1 * part + ex2 * (1-part)) * zf),
						     (int)((y + ey1 * part + ey2 * (1-part)) * zf));
	}


		// teken Sign, direction en Roadusers voor alpha-inward Lanes
		Sign s = null;
		int num_alpha_lanes = alphaLanes.length;
		for(int n = 0; n < num_alpha_lanes; n++)
		{
			float part2 = (float) n / width;
			int dx = 6 * ((sx1 - alphaNode.coord.x) / Math.abs(sx1 - alphaNode.coord.x));
			int dy = 6 * ((sy1 - alphaNode.coord.y) / Math.abs(sy1 - alphaNode.coord.y));
			s = alphaLanes[n].getSign();
			if (s != null && s.getType() == Sign.TRAFFICLIGHT)
			{
				g.setColor(s.mayDrive() ? Color.green : Color.red);
        		g.fillOval((int)((x + sx2 * part2 + sx1 * (1-part2) - 3 - dx) * zf),
						(int)((y + sy2 * part2 + sy1 * (1-part2) - 3 - dy) * zf),
							6,6);
			}


			//calculate the points a and b
			int signx=0, signy=0;
			switch (acp){
				case 0: signx =  1;signy =  0;break;
				case 1: signx =  0;signy =  1;break;
				case 2: signx = -1;signy =  0;break;
				case 3: signx =  0;signy = -1;break;
			}
			float ax = (float) sx1 + (float) signx * (float) 10.0 * ((float) 0.5+n);
			float ay = (float) sy1 + (float) signy * (float) 10.0 * ((float) 0.5+n);

			switch (bcp){
				case 0: signx = -1;signy =  0;break;
				case 1: signx =  0;signy = -1;break;
				case 2: signx =  1;signy =  0;break;
				case 3: signx =  0;signy =  1;break;
			}

			float bx = (float) ex1 + (float) signx * (float) 10.0 * ((float) 0.5+n);
 			float by = (float) ey1 + (float) signy * (float) 10.0 * ((float) 0.5+n);


			g.setColor(Color.black);
			float p1 = (float)0.6 / (getLength()+1);
			float p2 = (float)1.5 / (getLength()+1);
			int ax1 = (int)((x+(ax==bx?(int)ax:Math.round(ax * (1-p1) + bx * p1))) * zf);
			int ay1 = (int)((y+(ay==by?(int)ay:Math.round(ay * (1-p1) + by * p1))) * zf);
			int ax2 = (int)((x+(ax==bx?(int)ax:Math.round(ax * (1-p2) + bx * p2))) * zf);
			int ay2 = (int)((y+(ay==by?(int)ay:Math.round(ay * (1-p2) + by * p2))) * zf);
			g.drawLine(ax1,ay1,ax2,ay2);

	    boolean[] tars = alphaLanes[n].getTargets();
			if(tars[0])	// left
			{
				g.drawLine(ax1,ay1,ax1+(int)(signx * 3 * zf),ay1+(int)(signy * 3 * zf));
			}
			if(tars[1]) // straight
			{
	  		float p3 = (float)0.2 / (getLength()+1);
	  		g.drawLine(ax1,ay1,(int)((x+(ax==bx?(int)ax:Math.round(ax * (1-p3) + bx * p3))) * zf),
	  	  		               (int)((y+(ay==by?(int)ay:Math.round(ay * (1-p3) + by * p3))) * zf));
			}
			if(tars[2]) // right
			{
				g.drawLine(ax1,ay1,ax1-(signx * 3),ay1-(signy * 3));
 			}


			LinkedList queue = alphaLanes[n].getQueue();	// get roadusers
			if (!queue.isEmpty())
			{
				int pos = -1;
				Roaduser ru = null;
				ListIterator li = queue.listIterator();
				Drivelane lane = alphaLanes[n];
				while(li.hasNext())
				{
					ru = (Roaduser) li.next();
					pos = ru.getPosition();
					if (pos >= length) {
								Point tp1 = new Point((int)(betaCoord.x - 5 * Math.cos(betaAngle)),
								                      (int)(betaCoord.y + 5 * Math.sin(betaAngle)));
								Point p = avg(getCurvePoint(tp1,betaAngle,-1,n,x,y),
								              getCurvePoint(tp1,betaAngle,-1,n+1,x,y));
                               paintRoaduserOnNode(g,ru,betaNode,p,pos-length,alphaLanes[n]);
 					}
					else {
						float part = (float) (pos+1) / (length+1);
						ru.paint(g, x+(ax==bx?(int)ax:Math.round(ax * (1-part) + bx * part)), y+(ay==by?(int)ay:Math.round(ay * (1-part) + by * part)), zf,Math.toRadians(180)+dlangle);
					}
		 		}
		 	}
    }
		// teken Sign en Roadusers voor beta-inward Lanes
		int num_beta_lanes = betaLanes.length;
    		for(int n = 0; n < num_beta_lanes; n++)
		{
			int dx = 6 * ((ex2 - betaNode.coord.x) / Math.abs(ex2 - betaNode.coord.x));
			int dy = 6 * ((ey2 - betaNode.coord.y) / Math.abs(ey2 - betaNode.coord.y));
			s = betaLanes[n].getSign();
			if(s != null && s.getType() == Sign.TRAFFICLIGHT)
			{

				float part = (float) (num_beta_lanes-n-1)/getWidth();
				g.setColor(s.mayDrive() ? Color.green : Color.red);
				g.fillOval((int)((x + ex1 * part + ex2 * (1-part) - 3 - dx) * zf),
							(int)((y + ey1 * part + ey2 * (1-part) - 3 - dy) * zf),
								6,6);
//							(int)(5 * zf),(int)(5 * zf));
     		}

			// calculate the points a and b
			int signx=0, signy=0;
			switch (acp){
				case 0: signx = -1;signy =  0;break;
				case 1: signx =  0;signy = -1;break;
				case 2: signx =  1;signy =  0;break;
				case 3: signx =  0;signy =  1;break;
			}
			float ax = (float) sx2 + (float) signx * (float) 10.0 * ((float) 0.5+(num_beta_lanes-n-1));
			float ay = (float) sy2 + (float) signy * (float) 10.0 * ((float) 0.5+(num_beta_lanes-n-1));

	switch (bcp){
		case 0: signx =  1;signy =  0;break;
		case 1: signx =  0;signy =  1;break;
		case 2: signx = -1;signy =  0;break;
		case 3: signx =  0;signy = -1;break;
	}
	float bx = (float) ex2 + (float) signx * (float) 10.0 * ((float) 0.5+(num_beta_lanes-n-1));
	float by = (float) ey2 + (float) signy * (float) 10.0 * ((float) 0.5+(num_beta_lanes-n-1));

	g.setColor(Color.black);
	float p1 = (float)0.6 / (getLength()+1);
	float p2 = (float)1.5 / (getLength()+1);
	int ax1 = (int)((x+(ax==bx?(int)ax:Math.round(ax * p1 + bx * (1-p1)))) * zf);
	int ay1 = (int)((y+(ay==by?(int)ay:Math.round(ay * p1 + by * (1-p1)))) * zf);
	int ax2 = (int)((x+(ax==bx?(int)ax:Math.round(ax * p2 + bx * (1-p2)))) * zf);
	int ay2 = (int)((y+(ay==by?(int)ay:Math.round(ay * p2 + by * (1-p2)))) * zf);
	g.drawLine(ax1,ay1,ax2,ay2);
	boolean[] tars = betaLanes[n].getTargets();
	if(tars[0])	// left
	{
		g.drawLine(ax1,ay1,ax1+(signx * 3),ay1+(signy * 3));
	}
	if(tars[1]) // straight
	{
 		float p3 = (float)0.2 / (getLength()+1);
  		g.drawLine(ax1,ay1,(int)((x+(ax==bx?(int)ax:Math.round(ax * p3 + bx * (1-p3)))) * zf),
  		             (int)((y+(ay==by?(int)ay:Math.round(ay * p3 + by * (1-p3)))) * zf));
	}
	if(tars[2]) // right
	{
		g.drawLine(ax1,ay1,ax1-(int)(signx * 3 * zf),ay1-(int)(signy * 3 * zf));
	}

			LinkedList queue = betaLanes[n].getQueue();	// get roadusers
			if (!queue.isEmpty())
			{
				Roaduser ru = null;
				int pos = -1;

				ListIterator li = queue.listIterator();
				while(li.hasNext())
				{
					ru = (Roaduser)li.next();
					pos = ru.getPosition();
					if (pos >= length) {
                              	Point tp1 = new Point((int)(alphaCoord.x - 5 * Math.cos(alphaAngle)),
					            					  (int)(alphaCoord.y + 5 * Math.sin(alphaAngle)));
								Point p = avg(getCurvePoint(tp1,alphaAngle,1,n,x,y),
								              getCurvePoint(tp1,alphaAngle,1,n+1,x,y));
                                paintRoaduserOnNode(g,ru,alphaNode,p,pos-length,betaLanes[n]);

					}
					else {
						float part = (float) (pos+1) / (getLength()+1);
						ru.paint(g, x+(ax==bx?(int)ax:Math.round(ax * part + bx * (1-part))), y+(ay==by?(int)ay:Math.round(ay * part + by * (1-part))), zf, dlangle);
					}
				}
			}
		}
        }
   }
}