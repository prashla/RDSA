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

package gld.sim;

import gld.Model;
import gld.GLDException;

import gld.algo.dp.*;
import gld.algo.tlc.*;
import gld.infra.*;
import gld.sim.stats.TrackerFactory;
import gld.utils.Arrayutils;
import gld.utils.NumberDispenser;
import gld.xml.*;

import java.awt.Point;
import java.awt.Color;
import java.io.IOException;
import java.util.*;

/**
 * 
 * The heart of the simulation.
 * 
 * @author Group Model
 * @version 1.0
 */

public class SimModel extends Model implements XMLSerializable
{
    /** The pseudo-random-number-generator we need in this simulation */
    protected Random generator;
    /** The second thread that runs the actual simulation */
    protected SimModelThread thread;
    /** The SimController */
    protected SimController controller;
    /** The current cycle we're in */
    protected int curCycle;
    /** The Driving Policy in this Simulation */
    protected static DrivingPolicy dp;
    /** The TrafficLightControlling Algorithm */
    protected TLController tlc;
    /** The Thing that makes all Trafficlights shudder */
    protected SignController sgnctrl;
    /** Name of the simulation */
    protected String simName = "untitled";
    /** A boolean to keep track if this sim has already run (ivm initialization) */
    protected boolean hasRun = false;
    /** Indicates if roadusers cross nodes or jump over them. */
    public static boolean CrossNodes = true;
    /** Indicates whether we are running a series of simulations */
    protected boolean runSeries = false;
    protected boolean locked = false;
    /** The number of steps each of the simulations in a series should make */
    protected static int numSeriesSteps = 50000;
    protected static int LOCK_THRESHOLD = 10000;
    protected static int numSeries = 10;
    protected int curSeries = 0;

    // prashanth:
    protected int simulationDuration = 150;
    protected static int UNIFORM = 0;
    protected static int POISSON = 1;
    int inputDistribution = UNIFORM;

    /**
     * Creates second thread
     */
    public SimModel()
    {
        thread = new SimModelThread();
        thread.start();
        curCycle = 0;
        generator = new Random(1015);
        sgnctrl = new SignController(tlc, infra);
    }

    public void setSimController(SimController sc)
    {
        controller = sc;
    }

    public void setInfrastructure(Infrastructure i)
    {
        pause();

        super.setInfrastructure(i);
        if (tlc != null)
            tlc.setInfrastructure(i);
        if (sgnctrl != null)
            sgnctrl.setInfrastructure(i);
    }

    /** Returns the SIMULATION DURATION */
    public int getSimulationDuration()
    {
        return simulationDuration;
    }
    
    /** Returns the SIMULATION DURATION */
    public void setSimulationDuration(int simDur)
    {
        simulationDuration = simDur;
    }


    /** Returns the current cycle */
    public int getCurCycle()
    {
        return curCycle;
    }

    /** Sets the current cycle */
    public void setCurCycle(int c)
    {
        curCycle = c;
        tlc.setCurCycle(c);
    }

    /** Returns the current Driving Policy */
    public static DrivingPolicy getDrivingPolicy()
    {
        return dp;
    }

    /** Sets the current DrivTLController */
    public void setDrivingPolicy(DrivingPolicy _dp)
    {
        dp = _dp;
    }

    /** Returns the current TLController */
    public TLController getTLController()
    {
        return tlc;
    }

    /** Sets the current TLController */
    public void setTLController(TLController _tlc)
    {
        tlc = _tlc;
        sgnctrl.setTLC(tlc);
    }

    /** Returns the random number generator */
    public Random getRandom()
    {
        return generator;
    }

    /** Sets the random number generator */
    public void setRandom(Random r)
    {
        generator = r;
    }

    /** Returns the name of the simulation */
    public String getSimName()
    {
        return simName;
    }

    /** Sets the name of the simulation */
    public void setSimName(String s)
    {
        simName = s;
    }

    /** Returns the pseudo-random-number generator of this Model */
    public Random getRNGen()
    {
        return generator;
    }

    /** Sets spawn frequency for given node and ru type. */
    public void setSpawnFrequency(EdgeNode en, int rutype, float newspawn)
    {
        en.setSpawnFrequency(rutype, newspawn);
        setChanged();
        notifyObservers();
    }

    /**
     * Stops the simulation. This should only be called when the program exits.
     * To start a new simulation, the simulation should be paused with a call to
     * pause(), then followed by a call to reset(), and finally resumed with
     * unpause().
     */
    public void stop()
    {
        sgnctrl.callFinish();
        thread.die();
    }

    /**
     * Pauses the simulation
     */
    public void pause()
    {
        thread.pause();
    }

    /**
     * Unpauses the simulation
     */
    public void unpause()
    {
        thread.unpause();
    }

    public boolean isRunning()
    {
        return thread.isRunning();
    }

    public void runSeries()
    {
        curSeries = 0;
        runSeries = true;
        nextSeries();
    }

    public void nextSeries()
    {
        controller.nextSeries();
    }

    public void lockedSeries()
    {
        pause();
        for (; curCycle < numSeriesSteps; curCycle++)
        {
            setChanged();
            notifyObservers();
        }
        locked = false;
        nextSeries();
    }

    public void stopSeries()
    {
        curSeries = 0;
        runSeries = false;
    }

    public void nextCurSeries()
    {
        curSeries++;
    }

    public int getCurSeries()
    {
        return curSeries;
    }

    public boolean isRunSeries()
    {
        return runSeries;
    }

    public int getNumSeries()
    {
        return numSeries;
    }

    /**
     * Resets data
     */
    public void reset() throws SimulationRunningException
    {
        if (thread.isRunning())
            throw new SimulationRunningException(
                    "Cannot reset data while simulation is running.");
        infra.reset();
        tlc.reset();
        curCycle = 0;
        generator = new Random(1015);
        TrackerFactory.resetTrackers();

        setChanged();
        notifyObservers();
    }

    /**
     * Does 1 step in the simulation. All cars move, pedestrians get squashed
     * etc...
     */
    public void doStep()
    {
        curCycle++;
        tlc.setCurCycle(curCycle);

        if (!hasRun)
        {
            initialize();
            hasRun = true;
        }

        cityDoStep();

        setChanged();
        notifyObservers();
        // ///////////////////////////////////////////////////////////////////////
        if (curCycle == simulationDuration)
            pause();
        // ///////////////////////////////////////////////////////////////////////
        if (runSeries && curCycle >= numSeriesSteps)
        {
            nextSeries();
        }
        if (locked && runSeries)
        {
            lockedSeries();
        }
    }

    public void initialize()
    {
        SAVE_STATS = true;
        Enumeration e = Arrayutils.getEnumeration(infra.getSpecialNodes());
        while (e.hasMoreElements())
            ((SpecialNode) (e.nextElement())).start();
    }

    /** Gets the speed of the simulation */
    public int getSpeed()
    {
        return thread.getSleepTime();
    }

    /** Sets the speed of the simulation */
    public void setSpeed(int s)
    {
        thread.setSleepTime(s);
    }

    protected void cityDoStep()
    {
        try
        {
            specialNodesDoStep();
            moveAllRoadusers();
            spawnNewRoadusers();
            if (curCycle < simulationDuration)
                sgnctrl.callUpdateQ();
            else
                System.out.println("Reached cycle limit " + curCycle);
            sgnctrl.switchSigns();
        } catch (Exception e)
        {
            System.out.println("The simulator made a booboo:");
            System.out.println(e.getMessage());
            e.printStackTrace();
        }
    }

    protected void specialNodesDoStep() throws Exception
    {
        Enumeration specialNodes = Arrayutils.getEnumeration(infra
                .getSpecialNodes());
        while (specialNodes.hasMoreElements())
            ((SpecialNode) (specialNodes.nextElement())).doStep(this);
    }

    /**
     * 
     * moving all Roadusers to their new places
     * 
     * @author Blah... why put an author tag at every 10-line piece of code? ->
     *         Just because we can! MUHAHAahahahahah!
     */
    public void moveAllRoadusers() throws InfraException
    { // Line below is faster than the obvious alternative
        Enumeration lanes = Arrayutils.getEnumeration(infra
                .getAllInboundLanes().toArray());
        Drivelane lane;
        while (lanes.hasMoreElements())
        {
            lane = (Drivelane) (lanes.nextElement());
            // First you should check wheter they are already moved ......
            if (lane.getCycleMoved() != curCycle)
                moveLane(lane);
        }
    }

    /**
     * moving all roadusers from one lane to their new places
     * 
     * @author Jilles V, Arne K, Chaim Z and Siets el S
     * @param lane
     *            The lane whose roadusers should be moved
     * @param callingLanes
     *            Vector of drivelanes, for recursive internal use only, this
     *            parameter should have the value null, when called from the
     *            outside
     * @version 1.0
     */
    protected void moveLane(Drivelane lane) throws InfraException
    {
        LinkedList queue;
        ListIterator li;
        Drivelane sourceLane, destLane;
        Node node;
        Sign sign;
        Roaduser ru;
        int ru_pos, ru_des, ru_speed, ru_type, ru_len;

        sign = lane.getSign();
        queue = lane.getQueue();
        li = queue.listIterator();

        while (li.hasNext())
        {
            try
            {
                ru = (Roaduser) li.next();
            }
            // When this exception is thrown you removed the first element of
            // the queue, therefore re-create the iterator.
            catch (Exception e)
            {
                li = queue.listIterator();
                continue;
            }

            if (!ru.didMove(curCycle))
            { // Only attempt to move this RU when it hasnt already
                ru.setCycleAsked(curCycle);
                node = sign.getNode();
                ru_pos = ru.getPosition();
                ru_speed = ru.getSpeed();
                ru_len = ru.getLength();
                ru_type = ru.getType();
                ru_des = ru.getDestNode().getId();

                PosMov[] posMovs = calcPosMovs(node, sign, lane, ru, li);

                ru.setInQueueForSign(false);
                if (lane.getFreeUnitsInFront(ru) < ru_speed)
                    ru.setInQueueForSign(true);

                if (ru_pos - ru_speed < 0)
                { // Handle Roadusers that possibly can cross a Node
                    if (node instanceof SpecialNode)
                    { // Handle Roadusers that get to Special Nodes
                        if (ru_pos == 0
                                || 0 == lane.getPosFree(li, 0, ru_len,
                                        ru_speed, ru))
                        {
                            node.processStats(ru, curCycle, sign);
                            ru.setPosition(-1);
                            ru.setPrevSign(-1);
                            li.remove();
                            ru
                                    .setWaitPos(sign.getId(), sign.mayDrive(),
                                            ru_pos);
                            ru.setInQueueForSign(false);
                            tlc.updateRoaduserMove(ru, lane, sign, ru_pos,
                                    null, null, 0, posMovs, null);
                            ((SpecialNode) (node)).enter(ru);
                            ru = null;
                        }
                    }
                    else if (lane.getSign().mayDrive())
                    { // Handle Roadusers that are (or nearly) at a Sign
                        if (ru_pos == 0
                                || 0 == lane.getPosFree(li, 0, ru_len,
                                        ru_speed, ru))
                        { // Can cross-check
                            destLane = dp.getDirection(ru, lane, node);
                            if (destLane.isLastPosFree(ru_len))
                            { // Check if there is room on the node
                                try
                                {
                                    node.processStats(ru, curCycle, sign); // Let
                                    // the
                                    // RU
                                    // touch
                                    // the
                                    // Sign
                                    // to
                                    // refresh/unload
                                    // some
                                    // statistical
                                    // data
                                    destLane.addRoaduserAtEnd(ru);
                                    ru.setPrevSign(lane.getSign().getId());
                                    li.remove(); // Remove the RU from the
                                    // present lane, and place it
                                    // on the destination lane
                                    ru.setWaitPos(sign.getId(),
                                            sign.mayDrive(), ru_pos);
                                    ru.setInQueueForSign(false);
                                    tlc.updateRoaduserMove(ru, lane, sign,
                                            ru_pos, destLane, destLane
                                                    .getSign(), ru
                                                    .getPosition(), posMovs,
                                            destLane);
                                } catch (Exception e)
                                {
                                    System.out
                                            .println("Something screwd up in SimModel.moveLane where a Roaduser is about to cross:"
                                                    + e);
                                }
                            }
                            else
                            { // Otherwise, check if the next lane should move,
                                // and then do just that
                                if (curCycle != destLane.getCycleAsked()
                                        && curCycle != destLane.getCycleMoved())
                                { // If the position is not free, then check if
                                    // it already moved this turn, if not:
                                    // moveLane(destLane); //
                                    // System.out.println("Waiting for another lane to move..");
                                }
                                if (destLane.isLastPosFree(ru_len))
                                { // Ok now the lane that should have moved,
                                    // moved so try again .........
                                    try
                                    {
                                        node.processStats(ru, curCycle, sign);
                                        destLane.addRoaduserAtEnd(ru);
                                        ru.setPrevSign(lane.getSign().getId());
                                        li.remove();
                                        ru.setWaitPos(sign.getId(), sign
                                                .mayDrive(), ru_pos);
                                        ru.setInQueueForSign(false);
                                        tlc.updateRoaduserMove(ru, lane, sign,
                                                ru_pos, destLane, destLane
                                                        .getSign(), ru
                                                        .getPosition(),
                                                posMovs, destLane);
                                    } catch (Exception e)
                                    {
                                    }
                                }
                                else
                                { // Apparently no space was created, so we're
                                    // still here.
                                    if (moveRoaduserOnLane(li, ru, ru_speed,
                                            lane) > 0)
                                    {
                                        ru.setWaitPos(sign.getId(), sign
                                                .mayDrive(), ru_pos);
                                    }
                                    tlc.updateRoaduserMove(ru, lane, sign,
                                            ru_pos, lane, sign, ru
                                                    .getPosition(), posMovs,
                                            null);
                                }
                            }
                        }
                    }
                    else
                    { /*
                       * Light==red, Try to move user as far as it can go on
                       * this lane. Update it's move.
                       */
                        if (moveRoaduserOnLane(li, ru, ru_speed, lane) > 0)
                        {
                            ru
                                    .setWaitPos(sign.getId(), sign.mayDrive(),
                                            ru_pos);
                        }
                        tlc.updateRoaduserMove(ru, lane, sign, ru_pos, lane,
                                sign, ru.getPosition(), posMovs, null);
                    }
                }
                else
                { /*
                   * Roaduser impossibly can cross a sign. The maximum amount of
                   * space per speed is travelled
                   */
                    if (moveRoaduserOnLane(li, ru, ru_speed, lane) > 0)
                    {
                        ru.setWaitPos(sign.getId(), sign.mayDrive(), ru_pos);
                    }
                    tlc.updateRoaduserMove(ru, lane, sign, ru_pos, lane, sign,
                            ru.getPosition(), posMovs, null);
                }

                if (ru != null)
                {
                    ru.setCycleMoved(curCycle);
                }
            }
        }
        lane.setCycleAsked(curCycle);
        lane.setCycleMoved(curCycle);
    }

    /* Moves a roaduser on it's present lane as far as it can go. */
    protected int moveRoaduserOnLane(ListIterator li, Roaduser ru,
            int speed_left, Drivelane lane)
    {
        int ru_pos = ru.getPosition();
        int ru_len = ru.getLength();
        int best_pos = ru_pos;
        int max_pos = ru_pos;
        int target_pos = (ru_pos - speed_left > 0) ? ru_pos - speed_left : 0;
        int waitsteps;
        // System.out.println("Targetpos:"+target_pos+" and hasPrev:"+li.hasPrevious());

        // Previous should be 'ru'
        Roaduser prv = (Roaduser) li.previous();

        if (prv == ru && li.hasPrevious())
        { /* has car in front */
            prv = (Roaduser) li.previous(); /* named prv */
            int prv_pos = prv.getPosition();
            max_pos = prv_pos + prv.getLength();
            if (max_pos < target_pos)
                best_pos = target_pos;
            else
                best_pos = max_pos;
            li.next();
            // System.out.println("RU had previous, now bestpos ="+best_pos);
        }
        else
            best_pos = target_pos;

        li.next();
        if (best_pos != ru_pos)
        { /* has no car in front, advance to your best pos */
            // The Roaduser can advance some positions
            ru.setPosition(best_pos);
            return (speed_left - (ru_pos - best_pos));
        }
        else
        { /* best_pos == ru_pos, or, you cant move. */
            return speed_left;
        }
    }

    protected PosMov[] calcPosMovs(Node node, Sign sign, Drivelane lane,
            Roaduser ru, ListIterator li)
    {
        // =======================================
        // Calculating the ranges per drivelane to where roaduser could get to
        // =======================================
        int ru_pos = ru.getPosition();
        int ru_speed = ru.getSpeed();
        int ru_len = ru.getLength();
        int ru_type = ru.getType();
        int ru_des = ru.getDestNode().getId();

        Vector vPosMovs = new Vector();
        int tlId = sign.getId();

        // Get the position closest to the Sign the RU can reach
        int bestPos = lane.getPosFree(li, ru_pos, ru_len, ru_speed, ru);
        for (int z = ru_pos; z >= bestPos; z--)
            vPosMovs.addElement(new PosMov(tlId, z));

        int speedLeft = ru_speed - ru_pos; // ru_pos as that is the number of
        // units to be moven to the Sign

        // Now figure out the other possible lanes
        if (bestPos == 0 && speedLeft > 0)
        {
            Drivelane[] possiblelanes = node.getShortestPaths(ru_des, ru_type);
            int lanes = possiblelanes.length;
            for (int j = 0; j < lanes; j++)
            {
                // For each possible lane
                Drivelane testLane = possiblelanes[j];

                if (testLane.isLastPosFree(ru_len))
                {
                    bestPos = -1;
                    speedLeft = speedLeft > ru_len ? speedLeft - ru_len : 0;
                    int worstPos = testLane.getCompleteLength() - ru_len;
                    int tltlId = testLane.getId();

                    // We kunnen ervanuitgaan dat we nooit 'echt' op de
                    // drivelane springen
                    // We kunnen wel 'naar' deze drivelane springen
                    // dwz, de posities op de node tot
                    // max(lane.length-ru_len,lane.clength-speedleft) zijn vrij
                    // om te gaan
                    bestPos = Math.max(testLane.getCompleteLength() - ru_len
                            - speedLeft, testLane.getLength() - ru_len);
                    for (int k = worstPos; k >= bestPos; k--)
                        vPosMovs.addElement(new PosMov(tltlId, k));
                }
            }
        }
        // Fuck it, we aint got the power to cross, so don't even bother
        // calculating further..
        PosMov[] posMovs = new PosMov[vPosMovs.size()];
        vPosMovs.copyInto(posMovs);
        vPosMovs = null;

        return posMovs;
    }

    /**
     * New road users are placed on the roads when necessary. When roads are
     * full, new road users are queued.
     */
    public void spawnNewRoadusers() throws InfraException,
            ClassNotFoundException
    {
        SpecialNode[] specialNodes = infra.getSpecialNodes();
        LinkedList wqueue;
        ListIterator list;
        EdgeNode edge;
        Roaduser r;
        int num_edges = specialNodes.length;
        int total_queue = 0;

        for (int i = 0; i < num_edges; i++)
        {
            if (!(specialNodes[i] instanceof EdgeNode))
                break;
            else
                edge = (EdgeNode) (specialNodes[i]);
            boolean placed = false;
            wqueue = edge.getWaitingQueue();
            int wqsize = wqueue.size();
            list = wqueue.listIterator();
            while (list.hasNext())
            {
                total_queue++;
                r = (Roaduser) list.next();
                if (placeRoaduser(r, edge))
                    list.remove();
            }

            SpawnFrequency[] freqs = edge.getSpawnFrequencies();
            DestFrequency[][] destfreqs = edge.getDestFrequencies();
            int num_freqs = freqs.length;
            int cur_index;
            int[] freqIndexes = new int[num_freqs];

            for (int nrs = 0; nrs < num_freqs; nrs++)
                freqIndexes[nrs] = nrs; // Shuffle the indexes

            Arrayutils.randomizeIntArray(freqIndexes, generator);

            for (int j = 0; j < num_freqs; j++)
            {
                // First try to place new road users on the road.
                cur_index = freqIndexes[j];
                if (inputDistribution == POISSON)
                {
                    double L = Math.exp(-freqs[cur_index].freq);
                    double p = 1.0;
                    int k = 0;

                    do
                    {
                        k++;
                        p *= Math.random();
                    } while (p > L);
                    for (int i1=1;i1<=k-1;i1++)
                    {
                        int ruType = freqs[cur_index].ruType;
                        /*
                         * Spawn road user of type freqs[i].ruType to a random
                         * destination. When all drivelanes are full the road
                         * users are queued.
                         */
                        SpecialNode dest = getRandomDestination(specialNodes,
                                edge, ruType, destfreqs);
                        r = RoaduserFactory.genRoaduser(ruType, edge, dest, 0);
                        r.setDrivelaneStartTime(curCycle);

                        // Add R in queue if there is no place
                        if (!placeRoaduser(r, edge))
                        {
                            if (wqsize < 374)
                            {
                                list.add(r);
                                wqsize++;
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                }
                else // Uniform spawn rates
                {
                    if (freqs[cur_index].freq >= generator.nextFloat())
                    {
                        int ruType = freqs[cur_index].ruType;
                        /*
                         * Spawn road user of type freqs[i].ruType to a random
                         * destination. When all drivelanes are full the road
                         * users are queued.
                         */
                        SpecialNode dest = getRandomDestination(specialNodes,
                                edge, ruType, destfreqs);
                        r = RoaduserFactory.genRoaduser(ruType, edge, dest, 0);
                        r.setDrivelaneStartTime(curCycle);

                        // Add R in queue if there is no place
                        if (!placeRoaduser(r, edge))
                        {
                            if (wqsize < 374)
                            {
                                list.add(r);
                                wqsize++;
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                }
            }
        }

        if (total_queue >= LOCK_THRESHOLD)
            locked = true;
    }

    /**
     * A road user is placed on the given edge node. When road is full the ru is
     * queued
     */
    private boolean placeRoaduser(Roaduser r, SpecialNode edge)
    {
        Drivelane found = findDrivelaneForRU(r, edge);
        if (found == null)
            return false;
        else
        {
            // There is room for me!
            try
            {
                // System.out.println("Adding RU with type:"+r.getType()+" to lane:"+found.getSign().getId()+" going to Node:"+found.getNodeLeadsTo().getId()+" at pos:"+found.getNodeLeadsTo().isConnectedAt(found.getRoad())+" with type:"+found.getType());
                found.addRoaduserAtEnd(r, found.getLength() - r.getLength());
                r.addDelay(curCycle - r.getDrivelaneStartTime());
                r.setDrivelaneStartTime(curCycle);
                return true;
            } catch (Exception e)
            {
                return false;
            }
        }
    }

    private Drivelane findDrivelaneForRU(Roaduser r, SpecialNode e)
    {
        SpecialNode dest = (SpecialNode) r.getDestNode();
        Drivelane[] lanes = (Drivelane[]) e.getShortestPaths(dest.getId(),
                r.getType()).clone();
        Arrayutils.randomizeArray(lanes);
        int num_lanes = lanes.length;
        for (int i = 0; i < num_lanes; i++)
        {
            if (lanes[i].isLastPosFree(r.getLength()))
                return lanes[i];
        }
        // System.out.println("Couldnt place RU");
        return null;
    }

    /** Get a completely random destination, don't choose moi */
    public SpecialNode getRandomDestination(SpecialNode moi)
            throws InfraException
    {
        SpecialNode[] dests = infra.getSpecialNodes();
        if (dests.length < 2)
            throw new InfraException(
                    "Cannot choose random destination. Not enough special nodes.");
        SpecialNode result;
        while (moi == (result = dests[(int) (generator.nextFloat() * dests.length)]))
            ;
        return result;
    }

    /* Choose a destination */
    private SpecialNode getRandomDestination(SpecialNode[] dests,
            SpecialNode here, int ruType, DestFrequency[][] destfreqs)
    {
        int[] destIds = here.getShortestPathDestinations(ruType);
        float choice = generator.nextFloat();
        float total = 0f;

        /*
         * All frequencies are between 0 and 1, but their total can be greater
         * than 1
         */
        for (int i = 0; i < destIds.length; i++)
            for (int j = 0; j < destfreqs[i].length; j++)
                if (destfreqs[destIds[i]][j].ruType == ruType)
                    total += destfreqs[destIds[i]][j].freq;

        float sumSoFar = 0f;
        int j = 0;
        int index = 0;
        boolean foundIndex = false;
        while (j < destIds.length && !foundIndex)
        {
            for (int i = 0; i < destfreqs[j].length; i++)
                if (destfreqs[destIds[j]][i].ruType == ruType)
                {
                    float now = (destfreqs[destIds[j]][i].freq) / total;
                    if (now + sumSoFar >= choice)
                    {
                        foundIndex = true;
                        index = j;
                    }
                    else
                        sumSoFar += now;
                }
            j++;
        }
        return dests[destIds[index]];
    }

    /* Get a random index out of the lanes */
    private int getRandomLaneNr(Drivelane[] lanes)
    {
        int ind = (int) Math.floor(generator.nextFloat() * (lanes.length));
        while (ind != lanes.length)
            ind = (int) Math.floor(generator.nextFloat() * (lanes.length));
        return ind;
    }

    /**
     * 
     * The second thread that runs the simulation.
     * 
     * @author Joep Moritz
     * @version 1.0
     */
    public class SimModelThread extends Thread
    {
        /** Is the thread suspended? */
        private volatile boolean suspended;
        /**
         * Is the thread alive? If this is set to false, the thread will die
         * gracefully
         */
        private volatile boolean alive;
        /** The time in milliseconds this thread sleeps after a call to doStep() */
        private int sleepTime = 100;

        /** Returns the current sleep time */
        public int getSleepTime()
        {
            return sleepTime;
        }

        /** Sets the sleep time */
        public void setSleepTime(int s)
        {
            sleepTime = s;
        }

        /**
         * Starts the thread.
         */
        public SimModelThread()
        {
            alive = true;
            suspended = true;
        }

        /**
         * Suspends the thread.
         */
        public synchronized void pause()
        {
            suspended = true;
        }

        /**
         * Resumes the thread.
         */
        public synchronized void unpause()
        {
            suspended = false;
            notify();
        }

        /**
         * Stops the thread. Invoked when the program exitst. This method cannot
         * be named stop().
         */
        public synchronized void die()
        {
            alive = false;
            interrupt();
        }

        /**
         * Returns true if the thread is not suspended and not dead
         */
        public boolean isRunning()
        {
            return !suspended && alive;
        }

        /**
         * Invokes Model.doStep() and sleeps for sleepTime milliseconds
         */

        public void run()
        {
            while (alive)
            {
                try
                {
                    sleep(sleepTime);
                    synchronized (this)
                    {
                        while (suspended && alive)
                            wait();
                    }
                    doStep();
                } catch (InterruptedException e)
                {
                }
            }
        }
    }

    // Some XMLSerializable stuff

    public void load(XMLElement myElement, XMLLoader loader)
            throws XMLTreeException, IOException, XMLInvalidInputException
    {
        super.load(myElement, loader);
        setInfrastructure(infra);
        Dictionary loadDictionary;
        try
        {
            loadDictionary = infra.getMainDictionary();
        } catch (InfraException ohNo)
        {
            throw new XMLInvalidInputException(
                    "This is weird. The infra can't make a dictionary for the second "
                            + "stage loader of the algorithms. Like : " + ohNo);
        }
        Dictionary infraDictionary = new Hashtable();
        infraDictionary.put("infra", infra);
        loadDictionary.put("infra", infraDictionary);
        boolean savedBySim = ("simulator").equals(myElement.getAttribute(
                "saved-by").getValue());
        if (savedBySim)
        {
            thread.setSleepTime(myElement.getAttribute("speed").getIntValue());
            simName = myElement.getAttribute("sim-name").getValue();
            curCycle = myElement.getAttribute("current-cycle").getIntValue();
            TLCFactory factory = new TLCFactory(infra);
            tlc = null;

            try
            {
                tlc = factory.getInstanceForLoad(factory
                        .getNumberByXMLTagName(loader.getNextElementName()));
                loader.load(this, tlc);
                System.out.println("Loaded TLC " + tlc.getXMLName());
            } catch (InfraException e2)
            {
                throw new XMLInvalidInputException(
                        "Problem while TLC algorithm was processing infrastructure :"
                                + e2);
            }
            tlc.loadSecondStage(loadDictionary);
            DPFactory dpFactory = new DPFactory(this, tlc);
            try
            {
                dp = dpFactory.getInstance(dpFactory
                        .getNumberByXMLTagName(loader.getNextElementName()));
                loader.load(this, dp);
                System.out.println("Loaded DP " + dp.getXMLName());
            } catch (ClassNotFoundException e)
            {
                throw new XMLInvalidInputException(
                        "Problem with creating DP in SimModel."
                                + "Could not generate instance of DP type :"
                                + e);
            }
            dp.loadSecondStage(loadDictionary);
            loader.load(this, sgnctrl);
            sgnctrl.setTLC(tlc);
        }
        else
        {
            curCycle = 0;
        }
        while (loader.getNextElementName().equals("dispenser"))
            loader.load(this, new NumberDispenser());
    }

    public XMLElement saveSelf()
    {
        XMLElement result = super.saveSelf();
        result.addAttribute(new XMLAttribute("sim-name", simName));
        result.addAttribute(new XMLAttribute("saved-by", "simulator"));
        result.addAttribute(new XMLAttribute("speed", thread.getSleepTime()));
        result.addAttribute(new XMLAttribute("current-cycle", curCycle));
        return result;
    }

    public void saveChilds(XMLSaver saver) throws IOException,
            XMLTreeException, XMLCannotSaveException
    {
        super.saveChilds(saver);
        System.out.println("Saving TLC " + tlc.getXMLName());
        saver.saveObject(tlc);
        System.out.println("Saving DP " + dp.getXMLName());
        saver.saveObject(dp);
        saver.saveObject(sgnctrl);
    }
}
