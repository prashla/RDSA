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
import gld.GLDException;
import gld.xml.*;

import java.awt.Point;
import java.awt.Graphics;
import java.awt.Color;
import java.io.IOException;
import java.util.*;

/**
 * 
 * Basic junction. A Node either is a Junction or an EdgeNode.
 * 
 * @author Group Datastructures
 * @version 1.0
 */

public class Junction extends Node
{
    /** The type of this node */
    protected static final int type = Node.JUNCTION;
    /** The width/height of this node in drivelanes. Needed to draw the node */
    protected int width = 4;
    /** A vector containing all roads connected to this node */
    protected Road[] allRoads =
    { null, null, null, null };
    /** A vector containing all roads that have this node as their Alpha node */
    protected Road[] alphaRoads = {};
    /** A Vector containing all Signs on this node */
    protected Sign[] signs = {};
    /**
     * Contains all possible combinations of signs which may be turned green at
     * the same time
     */
    protected Sign[][] signconfigs =
    { {} };

    // Prashanth: for FixedExtTLC
    protected int configIndex;

    // Prashanth: actionIndices per junction for QTLCLinearFA
    protected int[][] actionIndices =
    { {} };
    
    // Prashanth: switchPeriod per junction for FixedPriority
    int switchPeriod = 12;
    int periodIndex;

    public int getSwitchPeriod()
    {
        return switchPeriod;
    }

    public void setSwitchPeriod(int k)
    {
        switchPeriod = k;
    }
    
    /**
     * Temporary data structure to tranfer info from the first stage loader to
     * the second stage loader
     */
    protected TwoStageLoaderData loadData = new TwoStageLoaderData();

    /** Creates an empty junction (for loading) */
    public Junction()
    {
    }

    /**
     * Creates a new standard Junction
     * 
     * @param _coord
     *            The coordinates of this node on the map in pixels.
     */
    public Junction(Point _coord)
    {
        super(_coord);
    }

    /* ============================================ */
    /* Basic GET and SET methods */
    /* ============================================ */

    /** Returns the type of this junction */
    public int getType()
    {
        return type;
    }

    /** Returns the name of this junction. */
    public String getName()
    {
        return "Junction " + nodeId;
    }

    /** Returns all roads connected to this node */
    public Road[] getAllRoads()
    {
        return allRoads;
    }

    /** Sets the array that contains all roads connected to this node */

    public void setAllRoads(Road[] r) throws InfraException
    {
        allRoads = r;
        updateLanes();
    }

    /** Returns the alpha roads connected to this node */
    public Road[] getAlphaRoads()
    {
        return alphaRoads;
    }

    /** Sets the alpha roads connected to this node */
    public void setAlphaRoads(Road[] r) throws InfraException
    {
        alphaRoads = r;
        updateLanes();
    }

    /** Returns all possible sign configurations */
    public Sign[][] getSignConfigs()
    {
        return signconfigs;
    }

    /** Returns configIndex: for FixedExtTLC */
    public int getConfigIndex()
    {
        return configIndex;
    }

    /** Returns all possible actionIndices */
    public int[][] getActionIndices()
    {
        return actionIndices;
    }

    /** Sets all possible sign configurations */
    public void setSignConfigs(Sign[][] s)
    {
        signconfigs = s;
    }

    /** Sets configIndex */
    public void setConfigIndex(int s)
    {
        configIndex = s;
    }

    /** Sets all possible sign configurations */
    public void setActionIndices(int[][] s)
    {
        actionIndices = s;
    }

    /** Returns the signs on this Node */
    public Sign[] getSigns()
    {
        return signs;
    }

    /** Sets the signs on this Node */
    public void setSigns(Sign[] s) throws InfraException
    {
        signs = s;
    }

    /** Returns the width of this Node in number of lanes */
    public int getWidth()
    {
        return width;
    }

    /** Sets the width of this Node in number of lanes */
    public void setWidth(int max)
    {
        width = max;
    }

    /* ============================================ */
    /* STATISTICS */
    /* ============================================ */

    protected int calcDelay(Roaduser ru, int stop, int distance)
    {
        // calculate the delay for the drivelane leading to this Junction
        int start = ru.getDrivelaneStartTime();
        int speed = ru.getSpeed();

        int min_steps = distance / speed;
        int num_steps = stop - start;

        int delay = num_steps - min_steps;

        ru.addDelay(delay);
        // System.out.println("Just crossed, and waited:"+delay+" start:"+start+" stop:"+stop+" dist:"+distance);
        ru.setDrivelaneStartTime(stop);

        // then return the delay
        return delay;
    }

    /* ============================================ */
    /* MODIFYING DATA */
    /* ============================================ */

    public void reset()
    {
        super.reset();
        for (int i = 0; i < alphaRoads.length; i++)
        {
            alphaRoads[i].reset();
        }
        for (int i = 0; i < signs.length; i++)
        {
            signs[i].reset();
        }
    }

    public void addRoad(Road r, int pos) throws InfraException
    {
        if (r == null)
            throw new InfraException("Parameter r is null");
        if (pos > 3 || pos < 0)
            throw new InfraException("Position out of range");
        if (allRoads[pos] != null)
            throw new InfraException("Road already connected to position "
                    + pos);
        allRoads[pos] = r;
        Node other = r.getOtherNode(this);
        if (other == null || !other.isAlphaRoad(r))
            alphaRoads = (Road[]) Arrayutils.addElement(alphaRoads, r);
        updateLanes();
        calculateWidth();
    }

    public void setAlphaRoad(int pos) throws InfraException
    {
        if (pos > 3 || pos < 0)
            throw new InfraException("Position out of range");
        if (allRoads[pos] == null)
            throw new InfraException("No road is conencted at position " + pos);
        alphaRoads = (Road[]) Arrayutils.addElementUnique(alphaRoads,
                allRoads[pos]);
        updateLanes();
    }

    public void remRoad(int pos) throws InfraException
    {
        if (pos > 3 || pos < 0)
            throw new InfraException("Position out of range");
        Road road = allRoads[pos];
        if (road == null)
            throw new InfraException("Road at position " + pos
                    + " does not exist");
        allRoads[pos] = null;
        alphaRoads = (Road[]) Arrayutils.remElement(alphaRoads, road);
        updateLanes();
        calculateWidth();
    }

    public void remRoad(Road r) throws InfraException
    {
        if (r == null)
            throw new InfraException("Parameter r is null");
        alphaRoads = (Road[]) Arrayutils.remElement(alphaRoads, r);
        for (int i = 0; i < 4; i++)
        {
            if (allRoads[i] == r)
            {
                allRoads[i] = null;
                updateLanes();
                calculateWidth();
                return;
            }
        }
        throw new InfraException("Road not found in this node");
    }

    public void remAllRoads() throws InfraException
    {
        for (int i = 0; i < allRoads.length; i++)
            allRoads[i] = null;
        alphaRoads = new Road[0];
        updateLanes();
        calculateWidth();
    }

    /** Adds a sign configuration */
    public void addSignconfig(Sign[] conf) throws InfraException
    {
        if (conf == null)
            throw new InfraException("Parameter conf is null");
        signconfigs = (Sign[][]) Arrayutils.addElement(signconfigs, conf);
    }

    /** Removes a sign configuration */
    public void remSignconfig(Sign[] conf) throws InfraException
    {
        if (conf == null)
            throw new InfraException("Parameter conf is null");
        int i = Arrayutils.findElement(signconfigs, conf);
        if (i == -1)
            throw new InfraException("Sign configuration is not in the list");
        signconfigs = (Sign[][]) Arrayutils.remElement(signconfigs, i);
    }

    /* ============================================ */
    /* SMALL GET */
    /* ============================================ */

    public boolean isAlphaRoad(Road r) throws InfraException
    {
        if (r == null)
            throw new InfraException("Parameter r == null");
        for (int i = 0; i < alphaRoads.length; i++)
        {
            if (alphaRoads[i] == r)
                return true;
        }
        return false;
    }

    public boolean isConnected(Road r) throws InfraException
    {
        return isConnectedAt(r) != -1;
    }

    public int isConnectedAt(Road r) throws InfraException
    {
        if (r == null)
            throw new InfraException("Parameter r == null");
        if (allRoads[0] == r)
            return 0;
        if (allRoads[1] == r)
            return 1;
        if (allRoads[2] == r)
            return 2;
        if (allRoads[3] == r)
            return 3;
        throw new InfraException("Road is not connected to this node");
    }

    public boolean isConnectionPosFree(int pos) throws InfraException
    {
        if (pos > 3 || pos < 0)
            throw new InfraException("Position out of range");
        return (allRoads[pos] == null);
    }

    public int getNumRoads()
    {
        int i = 0;
        if (allRoads[0] != null)
            i++;
        if (allRoads[1] != null)
            i++;
        if (allRoads[2] != null)
            i++;
        if (allRoads[3] != null)
            i++;
        return i;
    }

    public int getNumAlphaRoads()
    {
        return alphaRoads.length;
    }

    public int getNumInboundLanes() throws InfraException
    {
        int num = 0;
        for (int i = 0; i < allRoads.length; i++)
        {
            if (allRoads[i] != null)
                num += allRoads[i].getNumInboundLanes(this);
        }
        return num;
    }

    public int getNumOutboundLanes() throws InfraException
    {
        int num = 0;
        for (int i = 0; i < allRoads.length; i++)
        {
            if (allRoads[i] != null)
                num += allRoads[i].getNumOutboundLanes(this);
        }
        return num;
    }

    public int getNumAllLanes()
    {
        int num = 0;
        for (int i = 0; i < allRoads.length; i++)
        {
            if (allRoads[i] != null)
                num += allRoads[i].getNumAllLanes();
        }
        return num;
    }

    public int getNumSigns()
    {
        return signs.length;
    }

    public int getNumRealSigns()
    {
        int c = 0;
        for (int i = 0; i < signs.length; i++)
            if (signs[i].getType() != Sign.NO_SIGN)
                c++;
        return c;
    }

    public int getDesiredSignType() throws InfraException
    {
        return getNumRoads() > 2 ? Sign.TRAFFICLIGHT : Sign.NO_SIGN;
    }

    /* ============================================ */
    /* LARGE GET */
    /* ============================================ */

    public Drivelane[] getLanesLeadingTo(Drivelane lane, int ruType)
            throws InfraException
    {
        Road road = lane.getRoad();
        // Road[] which will contain the Roads of this Node in a sorted fashion:
        // [0] == the drivelanes on this Road will have to turn left to get to
        // 'road', ..
        Road[] srt_rarr = new Road[3];

        if (allRoads[0] == road)
        {
            srt_rarr[0] = allRoads[3]; // Must turn left
            srt_rarr[1] = allRoads[2]; // Must go straight on
            srt_rarr[2] = allRoads[1]; // Must turn right
        }
        else if (allRoads[1] == road)
        {
            srt_rarr[0] = allRoads[0];
            srt_rarr[1] = allRoads[3];
            srt_rarr[2] = allRoads[2];
        }
        else if (allRoads[2] == road)
        {
            srt_rarr[0] = allRoads[1];
            srt_rarr[1] = allRoads[0];
            srt_rarr[2] = allRoads[3];
        }
        else
        {
            srt_rarr[0] = allRoads[2];
            srt_rarr[1] = allRoads[1];
            srt_rarr[2] = allRoads[0];
        }

        Vector v = new Vector();
        Drivelane[] lanes;
        int num_lanes;
        int cnt_lanes = 0;
        boolean[] targets;

        for (int i = 0; i < 3; i++)
        {
            if (srt_rarr[i] != null)
            {
                lanes = srt_rarr[i].getInboundLanes(this);
                num_lanes = lanes.length;
                for (int j = 0; j < num_lanes; j++)
                {
                    Drivelane l = lanes[j];
                    targets = l.getTargets();
                    if (targets[i] == true && l.mayUse(ruType))
                    {
                        v.addElement(l);
                        cnt_lanes++;
                    }
                }
            }
        }

        return (Drivelane[]) v.toArray(new Drivelane[cnt_lanes]);
    }

    /* Needs Testing! */
    public Drivelane[] getLanesLeadingFrom(Drivelane lane, int ruType)
            throws InfraException
    {
        Road road = lane.getRoad();
        // Road[] which will contain the Roads of this Node in a sorted fashion:
        // [0] == the drivelanes on this Road will have to turn left to get to
        // 'road', ..
        Road[] srt_rarr = new Road[3];

        if (allRoads[0] == road)
        {
            srt_rarr[0] = allRoads[1]; // Must turn left
            srt_rarr[1] = allRoads[2]; // Must go straight on
            srt_rarr[2] = allRoads[3]; // Must turn right
        }
        else if (allRoads[1] == road)
        {
            srt_rarr[0] = allRoads[2];
            srt_rarr[1] = allRoads[3];
            srt_rarr[2] = allRoads[0];
        }
        else if (allRoads[2] == road)
        {
            srt_rarr[0] = allRoads[3];
            srt_rarr[1] = allRoads[0];
            srt_rarr[2] = allRoads[1];
        }
        else
        {
            srt_rarr[0] = allRoads[0];
            srt_rarr[1] = allRoads[1];
            srt_rarr[2] = allRoads[2];
        }

        // System.out.println("Junction getLanesLeadingFrom "+nodeId);
        Vector v = new Vector();
        Drivelane[] lanes;
        int num_lanes;
        int cnt_lanes = 0;
        boolean[] targets = lane.getTargets();

        for (int i = 0; i < 3; i++)
        {
            if (srt_rarr[i] != null && targets[i] == true)
            {
                // System.out.println("Road at target:"+i+" isnt null, getting Outboundlanes");
                lanes = srt_rarr[i].getOutboundLanes(this);
                num_lanes = lanes.length;
                // System.out.println("Num lanes :"+num_lanes);
                for (int j = 0; j < num_lanes; j++)
                {
                    Drivelane l = lanes[j];
                    // System.out.println("Lane"+j+" being checked now. Has type:"+l.getType());
                    if (l.mayUse(ruType))
                    {
                        v.addElement(l);
                        cnt_lanes++;
                    }
                }
            }
        }
        return (Drivelane[]) v.toArray(new Drivelane[cnt_lanes]);
    }

    public Drivelane[] getLanesLeadingFrom(Drivelane lane)
            throws InfraException
    {
        Road road = lane.getRoad();
        // Road[] which will contain the Roads of this Node in a sorted fashion:
        // [0] == the drivelanes on this Road will have to turn left to get to
        // 'road', ..
        Road[] srt_rarr = new Road[3];

        if (allRoads[0] == road)
        {
            srt_rarr[0] = allRoads[1]; // Must turn left
            srt_rarr[1] = allRoads[2]; // Must go straight on
            srt_rarr[2] = allRoads[3]; // Must turn right
        }
        else if (allRoads[1] == road)
        {
            srt_rarr[0] = allRoads[2];
            srt_rarr[1] = allRoads[3];
            srt_rarr[2] = allRoads[0];
        }
        else if (allRoads[2] == road)
        {
            srt_rarr[0] = allRoads[3];
            srt_rarr[1] = allRoads[0];
            srt_rarr[2] = allRoads[1];
        }
        else
        {
            srt_rarr[0] = allRoads[0];
            srt_rarr[1] = allRoads[1];
            srt_rarr[2] = allRoads[2];
        }

        // System.out.println("Junction getLanesLeadingFrom "+nodeId);
        Vector v = new Vector();
        Drivelane[] lanes;
        int num_lanes;
        int cnt_lanes = 0;
        boolean[] targets = lane.getTargets();

        for (int i = 0; i < 3; i++)
        {
            if (srt_rarr[i] != null && targets[i] == true)
            {
                // System.out.println("Road at target:"+i+" isnt null, getting Outboundlanes");
                lanes = srt_rarr[i].getOutboundLanes(this);
                num_lanes = lanes.length;
                // System.out.println("Num lanes :"+num_lanes);
                for (int j = 0; j < num_lanes; j++)
                {
                    v.addElement(lanes[j]);
                    cnt_lanes++;
                }
            }
        }
        return (Drivelane[]) v.toArray(new Drivelane[cnt_lanes]);
    }

    /** Returns an array of all outbound lanes on this junction */
    public Drivelane[] getOutboundLanes() throws InfraException
    {
        int pointer = 0;
        // System.out.println("NewNumOutboundLanes: "+getNumOutboundLanes());
        Drivelane[] lanes = new Drivelane[getNumOutboundLanes()];
        Drivelane[] temp;
        for (int i = 0; i < allRoads.length; i++)
        {
            if (allRoads[i] != null)
            {
                temp = allRoads[i].getOutboundLanes(this);
                System.arraycopy(temp, 0, lanes, pointer, temp.length);
                pointer += temp.length;
            }
        }
        return lanes;
    }

    /** Returns an array of all inbound lanes on this junction */
    public Drivelane[] getInboundLanes() throws InfraException
    {
        // System.out.println("Junction.getInboundLanes()");
        int pointer = 0;
        Drivelane[] lanes = new Drivelane[getNumInboundLanes()];
        Drivelane[] temp;
        for (int i = 0; i < allRoads.length; i++)
        {
            if (allRoads[i] != null)
            {
                temp = allRoads[i].getInboundLanes(this);
                System.arraycopy(temp, 0, lanes, pointer, temp.length);
                pointer += temp.length;
            }
        }
        return lanes;
    }

    /* clockwise order guaranteed */
    public Drivelane[] getAllLanes() throws InfraException
    {
        int pointer = 0;
        Drivelane[] lanes = new Drivelane[getNumAllLanes()];
        Drivelane[] temp;
        Road road;
        for (int i = 0; i < allRoads.length; i++)
        {
            road = allRoads[i];
            if (road != null)
            {
                temp = road.getInboundLanes(this);
                System.arraycopy(temp, 0, lanes, pointer, temp.length);
                pointer += temp.length;
                temp = road.getOutboundLanes(this);
                System.arraycopy(temp, 0, lanes, pointer, temp.length);
                pointer += temp.length;
            }
        }
        return lanes;
    }

    /**
     * Returns an array of all lanes connected to this node, in clock-wise
     * order, starting at the given lane
     */
    public Drivelane[] getAllLanesCW(Drivelane lane) throws InfraException
    {
        Drivelane[] lanes = getAllLanes(); // in clockwise order starting at
        // road 0, lane 0

        // find the starting-lane
        int i = Arrayutils.findElement(lanes, lane);
        if (i == -1)
            throw new InfraException("Lane is not on this node");

        // shift all the lanes i places and remove the i-th element
        Drivelane[] result = new Drivelane[lanes.length - 1];
        System.arraycopy(lanes, i + 1, result, 0, lanes.length - i - 1);
        System.arraycopy(lanes, 0, result, lanes.length - i - 1, i);

        return result;
    }

    /*
     * Returns whether or not the tails of the lanes not being accesible by the
     * given array of Signs are free
     */
    public boolean areOtherTailsFree(Sign[] mayUse)
    {
        boolean[] pos = new boolean[4];
        Road[] roads = getAllRoads();
        Drivelane lane;
        int thisRoad;
        boolean[] targets;

        for (int i = 0; i < 4; i++)
            pos[i] = true;

        int num_mayuse = mayUse.length;

        for (int i = 0; i < num_mayuse; i++)
        {
            lane = mayUse[i].getLane();
            try
            {
                thisRoad = isConnectedAt(lane.getRoad());
            } catch (InfraException e)
            {
                thisRoad = 0;
                System.out
                        .println("Something went wrong in areOtherTailsFree()");
            }
            targets = lane.getTargets();

            if (targets[0])
            {
                pos[(thisRoad + 1) % 4] = false;
            }
            else if (targets[1])
            {
                pos[(thisRoad + 2) % 4] = false;
            }
            else if (targets[2])
            {
                pos[(thisRoad + 3) % 4] = false;
            }
        }

        int num_check = 0;
        for (int i = 0; i < 4; i++)
        {
            if (pos[i] && roads[i] != null)
            {
                Drivelane[] check = new Drivelane[0];
                try
                {
                    check = roads[i].getOutboundLanes(this);
                } catch (Exception e)
                {
                    check = new Drivelane[0];
                    System.out
                            .println("Something went wrong in areOtherTailsFree() 2");
                }

                num_check = check.length;
                for (int j = 0; j < num_check; j++)
                {
                    if (!check[j].isTailFree())
                    {
                        return false;
                    }
                }
            }
        }
        return true;

    }

    public void paint(Graphics g) throws GLDException
    {
        paint(g, 0, 0, 1.0f, 0.0);
    }

    public void paint(Graphics g, int x, int y, float zf) throws GLDException
    {
        paint(g, x, y, zf, 0.0);
    }

    public void paint(Graphics g, int x, int y, float zf, double bogus)
            throws GLDException
    {
        // TODO: * tekenen status stoplichten
        int width = getWidth();
        g.setColor(Color.black);
        g.drawRect((int) ((coord.x + x - 5 * width) * zf),
                (int) ((coord.y + y - 5 * width) * zf),
                (int) (10 * width * zf), (int) (10 * width * zf));
        if (nodeId != -1)
            g.drawString("" + nodeId,
                    (int) ((coord.x + x - 5 * width) * zf) - 10,
                    (int) ((coord.y + y - 5 * width) * zf) - 3);
    }

    /** (Re)Calculates the width of this junction */
    public void calculateWidth()
    {
        Road road;
        width = 4;
        for (int i = 0; i < 4; i++)
        {
            road = allRoads[i];
            if (road != null && road.getWidth() > width)
                width = road.getWidth();
        }
    }

    /* ============================================ */
    /* LOAD and SAVE */
    /* ============================================ */

    public void load(XMLElement myElement, XMLLoader loader)
            throws XMLTreeException, IOException, XMLInvalidInputException
    {
        super.load(myElement, loader);
        width = myElement.getAttribute("width").getIntValue();
        alphaRoads = (Road[]) XMLArray.loadArray(this, loader);
        loadData.roads = (int[]) XMLArray.loadArray(this, loader);
        loadData.signconfigs = (int[][]) XMLArray.loadArray(this, loader);
        loadData.signs = (int[]) XMLArray.loadArray(this, loader);
    }

    public XMLElement saveSelf() throws XMLCannotSaveException
    {
        XMLElement result = super.saveSelf();
        result.setName("node-junction");
        result.addAttribute(new XMLAttribute("width", width));
        return result;
    }

    public void saveChilds(XMLSaver saver) throws XMLTreeException,
            IOException, XMLCannotSaveException
    {
        super.saveChilds(saver);
        XMLArray.saveArray(alphaRoads, this, saver, "alpha-roads");
        XMLArray.saveArray(getRoadIdArray(), this, saver, "roads");
        XMLArray.saveArray(getSignConfigIdArray(), this, saver, "sign-configs");
        XMLArray.saveArray(getSignIdArray(), this, saver, "signs");
    }

    protected int[] getSignIdArray()
    {
        int[] result = new int[signs.length];
        for (int t = 0; t < signs.length; t++)
        {
            if (signs[t] == null)
                result[t] = -1;
            else
                result[t] = signs[t].getId();
        }
        return result;
    }

    protected int[] getRoadIdArray()
    {
        int[] result = new int[allRoads.length];
        for (int t = 0; t < allRoads.length; t++)
        {
            if (allRoads[t] == null)
                result[t] = -1;
            else
                result[t] = allRoads[t].getId();
        }
        return result;
    }

    protected int[][] getSignConfigIdArray()
    {
        int[][] result;
        if (signconfigs.length == 0)
            result = new int[0][0];
        else
            result = new int[signconfigs.length][signconfigs[0].length];
        for (int t = 0; t < signconfigs.length; t++)
        {
            result[t] = new int[signconfigs[t].length];
            for (int u = 0; u < signconfigs[t].length; u++)
            {
                if (signconfigs[t][u] == null)
                    result[t][u] = -1;
                else
                    result[t][u] = signconfigs[t][u].getId();
            }
        }
        return result;
    }

    public String getXMLName()
    {
        return parentName + ".node-junction";
    }

    private class TwoStageLoaderData
    {
        int[] roads; // Storage for road-id's
        int[] signs; // For sign ids
        int[][] signconfigs; // For the sign configs
    }

    public void loadSecondStage(Dictionary dictionaries)
            throws XMLInvalidInputException, XMLTreeException
    {
        super.loadSecondStage(dictionaries);
        // Load roads
        Dictionary roadDictionary = (Dictionary) (dictionaries.get("road"));
        allRoads = new Road[loadData.roads.length];
        for (int t = 0; t < loadData.roads.length; t++)
        {
            allRoads[t] = (Road) (roadDictionary.get(new Integer(
                    loadData.roads[t])));
            if (allRoads[t] == null && loadData.roads[t] != -1)
                System.out.println("Warning : " + getName()
                        + " could not find road " + loadData.roads[t]);
        }
        // Load normal signs
        Dictionary laneDictionary = (Dictionary) (dictionaries.get("lane"));
        signs = new Sign[loadData.signs.length];
        for (int t = 0; t < loadData.signs.length; t++)
            signs[t] = getSign(laneDictionary, loadData.signs[t]);
        // Load Signconfigurations
        signconfigs = new Sign[loadData.signconfigs.length][2];
        for (int t = 0; t < signconfigs.length; t++)
        {
            signconfigs[t] = new Sign[loadData.signconfigs[t].length];
            for (int u = 0; u < signconfigs[t].length; u++)
            {
                signconfigs[t][u] = getSign(laneDictionary,
                        loadData.signconfigs[t][u]);
            }
        }
        // Tell *all* roads to load themselves
        // It's possible that this Node has a BetaLane that has not been
        // SecondStageLoaded
        // And so we cant do an UpdateLanes() as that one needs secondStageData
        // to proceed.
        // Hence, we need to 2ndStage all Roads.
        Enumeration e = new ArrayEnumeration(allRoads);
        Road tmpRoad;
        while (e.hasMoreElements())
        {
            tmpRoad = (Road) e.nextElement();
            if (tmpRoad != null)
                tmpRoad.loadSecondStage(dictionaries);
        }
        try
        { // System.out.println("Trying to updateLanes()");
            updateLanes();
        } catch (InfraException x)
        {
            throw new XMLInvalidInputException(
                    "Cannot initialize lanes of node " + nodeId);
        }
    }

    protected Sign getSign(Dictionary laneDictionary, int id)
    {
        Drivelane tmp = (Drivelane) (laneDictionary.get(new Integer(id)));
        if (tmp == null)
        {
            return null;
        }
        else
        {
            return tmp.getSign();
        }
    }
}