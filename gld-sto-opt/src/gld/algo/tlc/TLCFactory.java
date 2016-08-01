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

package gld.algo.tlc;

/** This class can be used to create instances of Traffic Light Controllers
 * for a specific infrastructure.
 */

import gld.infra.InfraException;
import gld.infra.Infrastructure;
import gld.utils.StringUtils;
import java.util.Dictionary;
import java.util.Hashtable;
import java.util.Random;

public class TLCFactory
{
    protected Infrastructure infra;
    protected Random random;

    protected static final int  
            PTLC = 0;

    protected static final String[] tlcDescs = { "Priority based traffic light controller" };

    protected static final String[] xmlNames = {
            PTLCL12345T123.shortXMLName
            };

    protected static final String[] categoryDescs = {
            "Graded TLC" };
    protected static final int[][] categoryTLCs = { { PTLC }
 };

    /**
     * Makes a new TLCFactory for a specific infrastructure with a new random
     * number generator.
     * 
     * @param infra
     *            The infrastructure
     */
    public TLCFactory(Infrastructure infra)
    {
        this.infra = infra;
        random = new Random(1015);

    }

    /**
     * Makes a new TLCFactory for a specific infrastructure
     * 
     * @param random
     *            The random number generator which some algorithms use
     * @param infra
     *            The infrastructure
     */
    public TLCFactory(Infrastructure infra, Random random)
    {
        this.infra = infra;
        this.random = random;
    }

    /**
     * Looks up the id of a TLC algorithm by its description
     * 
     * @param algoDesc
     *            The description of the algorithm
     * @returns The id of the algorithm
     * @throws NoSuchElementException
     *             If there is no algorithm with that description.
     */
    public static int getId(String algoDesc)
    {
        return StringUtils.getIndexObject(tlcDescs, algoDesc);
    }

    /** Returns an array of TLC descriptions */
    public static String[] getTLCDescriptions()
    {
        return tlcDescs;
    }

    /**
     * Look up the description of a TLC algorithm by its id
     * 
     * @param algoId
     *            The id of the algorithm
     * @returns The description
     * @throws NoSuchElementException
     *             If there is no algorithm with the specified id.
     */
    public static String getDescription(int algoId)
    {
        return (String) (StringUtils.lookUpNumber(tlcDescs, algoId));
    }

    /** Returns an array containing the TLC category descriptions. */
    public static String[] getCategoryDescs()
    {
        return categoryDescs;
    }

    /** Returns an array of TLC numbers for each TLC category. */
    public static int[][] getCategoryTLCs()
    {
        return categoryTLCs;
    }

    /** Returns the total number of TLCs currently available. */
    public static int getNumberOfTLCs()
    {
        return tlcDescs.length;
    }

    /** Gets the number of an algorithm from its XML tag name */
    public static int getNumberByXMLTagName(String tagName)
    {
        return StringUtils.getIndexObject(xmlNames, tagName);
    }

    /** Returns an instance of a TLC by its description. */
    public TLController genTLC(String tlcDesc) throws InfraException
    {
        return getInstanceForLoad(getId(tlcDesc));
    }

    public TLController genTLC(int cat, int tlc) throws InfraException
    {
        return getInstanceForLoad(categoryTLCs[cat][tlc]);
    }

    /**
     * Gets a new instance of an algorithm by its number. This method is meant
     * to be used for loading.
     */
    public TLController getInstanceForLoad(int algoId) throws InfraException
    {
        switch (algoId)
        {
        case PTLC:
            return new PTLCL12345T123(infra);
        }
        throw new InfraException(
                "The TLCFactory can't make TLC's of type " + algoId);
    }
}