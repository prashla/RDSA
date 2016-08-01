/*-----------------------------------------------------------------------
 *
 * Qlearning with linear FA
 *
 *------------------------------------------------------------------------*/

package gld.algo.tlc;

import gld.*;
import gld.sim.*;
import gld.algo.tlc.*;
import gld.infra.*;
import gld.utils.*;
import gld.xml.*;
import java.io.IOException;
import java.util.*;
import java.awt.Point;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * At any instant, PTLC assigns a priority value for each lane $k$ in a graded
 * fashion such that lanes with high congestion and longer elapsed times get the
 * highest priority, those with low congestion/elapsed time get lowest priority
 * and the configurations in-between are assigned intermediate priorities. The
 * sign configuration chosen is the one that maximizes the sum of priorities
 * over all lanes.
 * 
 * @version 1.0 31 Jul 2016
 * @author Prashanth L.A.
 * @see gld.tt.2RDSA_OuterLoop
 * @see gld.tt.2SPSA_OuterLoop
 * @see gld.tt.RDSARunner
 */

public class PTLCL12345T123 extends TLCBase
{
    public final static String shortXMLName = "PTLCL12345T123";

    // TLC vars
    protected int num_nodes;
    public static boolean first = true;

    // 1/n type stepsize for smoothing rewards
    double stepsize;

    /**
     * The thresholds: L1..L5 are for queue lengths and T1..T3 are for elapsed
     * times
     */
    double L1, L2, L3, L4, L5, T1, T2, T3;

    /**
     * Set the thresholds
     */
    public void setTheta(double theta[])
    {
        L1 = theta[0];
        L2 = theta[1];
        L3 = theta[2];
        L4 = theta[3];
        L5 = theta[4];

        T1 = theta[5];
        T2 = theta[6];
        T3 = theta[7];
    }

    /**
     * Get the thresholds
     */
    public double[] getTheta()
    {
        double thresholds[] = { L1, L2, L3, L4, L5, T1, T2, T3 };
        return thresholds;
    }

    /**
     * The constructor for TL controllers
     * 
     * @param The
     *            model being used.
     */
    public PTLCL12345T123(Infrastructure i)
    {
        super(i);
        System.out.println("PTLC initialized");
    }

    /**
     * Initialize the TLC algorithm
     */
    public void setInfrastructure(Infrastructure infra)
    {
        super.setInfrastructure(infra);
        try
        {
            num_nodes = infra.getAllNodes().length;
            averageCost = 0;
            stepsize = 1.0;
        }
        catch (Exception e)
        {
            System.out.println("PTLC failed to initialise.");
        }
    }

    /**
     * Calculates how every traffic light should be switched Per node, per sign
     * the waiting roadusers are passed and per each roaduser the gain is
     * calculated.
     * 
     * @param The
     *            TLDecision is a tuple consisting of a traffic light and a
     *            reward (Q) value, for it to be green
     * @see gld.algo.tlc.TLDecision
     */
    public TLDecision[][] decideTLs()
    {
        int num_lanes, numCycles, numUsers, gain = 0, lGain = 0, tGain = 0;

        // System.out.println("--------------------------Decide TLs:"
        // + getCurCycle()
        // + "----------------------------------------------------");

        updateNumElapsedCycles();
        clearNumGreenCycles();
        clearNumRedCycles();

        for (int i = 0; i < num_nodes; i++)
        {
            num_lanes = tld[i].length;
            for (int j = 0; j < num_lanes; j++)
            {
                numUsers = tld[i][j].getTL().getLane().getNumRoadusersWaiting();
                numCycles = tld[i][j].getNumRedCycles();
                if (numUsers < L1)
                {
                    lGain = 1;
                }
                else if (numUsers < L2)
                {
                    lGain = 2;
                }

                else if (numUsers < L3)
                {
                    lGain = 3;
                }
                else if (numUsers < L4)
                {
                    lGain = 4;
                }

                else if (numUsers < L5)
                {
                    lGain = 5;
                }
                else
                {
                    lGain = 6;
                }
                // Elapsed time - calculate gain
                if (numUsers < T1)
                {
                    tGain = 1;
                }
                else if (numUsers < T2)
                {
                    tGain = 2;
                }

                else if (numUsers < T3)
                {
                    tGain = 3;
                }
                else
                {
                    tGain = 4;
                }
                gain = lGain * tGain;
                tld[i][j].setGain(gain);
            }
        }

        return tld;
    }

    public void computeCost()
    {
        int num_lanes;
        cost = 0;

        double elapsedCost = 0.0;
        double queueCost = 0.0;
        int lpos = 0;
        for (int i = 0; i < num_nodes; i++)
        {
            num_lanes = tld[i].length;
            for (int j = 0; j < num_lanes; j++)
            {
                if (tld[i][j].getTL().getType() == Sign.TRAFFICLIGHT)
                {
                    elapsedCost += laneWeights[lpos]
                            * getCoarseElapsedTime(tld[i][j].getNumRedCycles());
                    queueCost += laneWeights[lpos]
                            * getCoarseQueueLength(tld[i][j].getTL().getLane()
                                    .getNumRoadusersWaiting());
                }
                lpos++;
                // reward -= tld[i][j].getNumGreenCycles();
            }
        }

//        System.out.println("QueueCost: " + queueCost + " ElapsedCost: " + elapsedCost);
        cost = queuesWeightInReward * queueCost
                + (1 - queuesWeightInReward) * elapsedCost;
        cost*=10;
    }

    double getCoarseQueueLength(int queueLength)
    {
        double coarseQueueLength;

        if (queueLength < L1)
        {
            coarseQueueLength = 0.0;
        }
        else if (queueLength < L2)
        {
            coarseQueueLength = 0.2;
        }

        else if (queueLength < L3)
        {
            coarseQueueLength = 0.4;
        }
        else if (queueLength < L4)
        {
            coarseQueueLength = 0.6;
        }

        else if (queueLength < L5)
        {
            coarseQueueLength = 0.8;
        }
        else
        {
            coarseQueueLength = 1.0;
        }
        return coarseQueueLength;
    }

    double getCoarseElapsedTime(int elapsedTime)
    {
        double coarseElapsedTime;

        if (elapsedTime < T1)
        {
            coarseElapsedTime = 0.0;
        }
        else if (elapsedTime < T2)
        {
            coarseElapsedTime = 0.33;
        }

        else if (elapsedTime < T3)
        {
            coarseElapsedTime = 0.66;
        }
        else
        {
            coarseElapsedTime = 1.0;
        }

        return coarseElapsedTime;
    }

    public void updateQ()
    {
        if (first)
        {
            first = false;
            return;
        }
        computeCost();
        averageCost = (1 - stepsize) * averageCost + stepsize * cost;

        // update step size
        stepsize = Math.pow(getCurCycle(), -1);
    }

    public void updateRoaduserMove(Roaduser _ru, Drivelane _prevlane,
            Sign _prevsign, int _prevpos, Drivelane _dlanenow, Sign _signnow,
            int _posnow, PosMov[] posMovs, Drivelane desired)
    { // No needed
    }

    // XMLSerializable implementation

    public XMLElement saveSelf() throws XMLCannotSaveException
    {
        XMLElement result = super.saveSelf();
        result.setName(shortXMLName);
        return result;
    }

    public String getXMLName()
    {
        return "model." + shortXMLName;
    }

}
