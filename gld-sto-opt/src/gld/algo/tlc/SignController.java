
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

import gld.*;
import gld.algo.tlc.*;
import gld.infra.*;
import gld.xml.*;
import gld.utils.*;

import java.io.IOException;
import java.util.Random;

/**
 *
 * This is the class for the sign controller. Here is decided how each sign should work.
 * Reward values are gathered for each sign, when set on 'Green' asking the applied TLController.
 * Non-TrafficLight signs can also be implemented when a Q-value algorithm is implemented for the appropiate
 * type of sign. After those values are gathered, this class will select the most rewarding 
 * trafficlightconfiguration for each Node and will set the Signs accordingly.
 * @see gld.infra.Sign
 *
 * @author Group Algorithms
 * @version 1.0
 */
public class SignController implements XMLSerializable
{
	/** The TLcontroller is used to gather Q-values and use them to set all the TLs */
	protected TLController tlcontroller;
	/**
	 * Indicates if trafficlights must be switched safely or not.
	 * Safely means no 2 trafficlights may be turned to green at the same
	 * time, if roadusers moving from the drivelanes behind those signs
	 * would collide on the node.
	 */
	public static boolean CrossNodesSafely = true;
	protected static final String shortXMLName = "signcontroller";
	protected Infrastructure infra;
	protected int num_nodes;
	//protected int[] chosenConfigs;
	protected Sign[][] currentSC;
	Random generator;
		
	
	public SignController(TLController t, Infrastructure i) {
		tlcontroller = t;
		infra = i;
		num_nodes = i.getNumNodes();
		generator = new Random(1015);	    

		/*chosenConfigs = new int[i.getAllNodes().length];
		for (int j=0; j<chosenConfigs.length; j++) chosenConfigs[j]=-1; // This should be done, otherwise all chosenConfigs are 0. This leads to problems, when config 0 is the first time the best config.*/
		currentSC = new Sign[num_nodes][0];
	}
	
	public Infrastructure getInfrastructure() { return infra; }

	public void setInfrastructure(Infrastructure i) 
	{ 
		infra = i; 
		num_nodes = i.getNumNodes();
		/*chosenConfigs = new int[i.getAllNodes().length];
		for (int j=0; j<chosenConfigs.length; j++) chosenConfigs[j]=-1;*/
		currentSC = new Sign[num_nodes][0];
	}
	public TLController getTLC() { return tlcontroller; }
	public void setTLC(TLController t) { tlcontroller = t; }
	
	/**
	 * Switch all the signs to their appropiate value.
	 */
	public void switchSigns()
	{
		TLDecision[][] decisions = tlcontroller.decideTLs();
		Node node = null;
		Node[] nodes = infra.getAllNodes();
		if (num_nodes>decisions.length)
			System.out.println("SignController switchSigns WARNING : "+"Less decisions than nodes !!!");
		for (int i=0; i < num_nodes; i++) {
			node = nodes[i];
			if (node.getType() == Node.JUNCTION) {
		 		if (decisions[i].length > 0) {
                    ///////////////////////////////////////////////////////////////////
		 			if (tlcontroller.getKeepSwitchControl()) {
		 				if (!(node.getKeepTLDFlag()))
		 					switchTrafficLights((Junction)node, decisions[i]);
		 			}
		 			else { 
		 				switchTrafficLights((Junction)node, decisions[i]);
		 				//System.out.println("SignController: switchcontrol off");
		 			}
		 			///////////////////////////////////////////////////////////////////
		 		}
				else
				if (node.getType() == Node.NON_TL )
				{
					switchNonTrafficLights((NonTLJunction)node, decisions[i]);
				}
			}
		}
	}
	
	/** 
	 * Switch the the non-TLsigns to their appropiate value according to normal traffic rules
	 * On a normal junction traffic from the right gets priority
	 * @param node The node involved.
	 */
	private void switchNonTrafficLights(NonTLJunction node, TLDecision[] dec)
	{
		//System.out.println("NonTLSwitch");
		//Sign[][] signConfigs = node.getSignConfigs() ;
		Sign[] signs = node.getSigns() ;
		boolean[] mayDrive = { false, false, false, false } ;
		//Check for traffic from the right
		Road[] roads = node.getAllRoads(); //In clockwise order I hope/think ;)
		Road right = roads[3] ;//Start with last road as being right
		Drivelane[] lanes ;
		boolean alpha = false ;
		
		for ( int i=0 ; i<roads.length ; i++ )//length is always 4 
		{
			alpha = right.getAlphaNode() == node ;//Check which are the incoming lanes
			lanes = alpha ? right.getAlphaLanes() : right.getBetaLanes() ;
			for ( int  j=0; j<lanes.length ; j++ ) {
				if ( lanes[j].getNumRoadusersWaiting() > 0 ) {
					mayDrive[i] = false ;//Traffic is coming from the right
				}
				else {
					mayDrive[i] = true ;
				}
			}	
		}
		
		Random random = new Random(1015) ;		
		boolean deadlock = true ;//All roads are waiting for eachother
		int choosenRoad = -1 ;
		//length is always 4
		for (int i=0; i<roads.length; i++)	{
			if(deadlock && mayDrive[i])	{
				deadlock = false ;
				choosenRoad = i ;
			}
		}
		if (deadlock)	{
			choosenRoad = (int)Math.floor(random.nextFloat()*4) ;
		}
		//Now switch the lights
		lanes = null ;
		try {
			lanes = node.getInboundLanes() ;
		}
		catch (InfraException e) {
			System.out.println("An error occured while setting non-Traffic lights") ;
		}
		
		for ( int i=0 ; i<lanes.length ; i++ ) {
			if ( lanes[i].getRoad() == roads[i] ) {
				signs[i].setState( true ) ;
			}
			else {
				signs[i].setState( false ) ;
			}
		}
	}
	
	/** 
	 * Switch the the TLsigns to their appropiate value
	 * The config with highest total gain is chosen
	 * @param dec The decision array consists of the generated values.
	 */
	public void switchTrafficLights(Junction node, TLDecision[] dec)
	{	//System.out.println("Switching TL's on junction "+node.getId());
		Sign[][] signConfs = node.getSignConfigs();
		Sign[] signs = node.getSigns();
		
		int num_sc = signConfs.length;
		Sign[] thisConfig;
		/*boolean flag = true;
		if(flag)
		{
			for (int i=0; i < num_sc; i++) {
				thisConfig = signConfs[i];
				System.out.print("[" + i + "] ->  ");
				for(int j=0; j < thisConfig.length; j++) {				
					System.out.print(thisConfig[j].getId() + ",");
				}			
				System.out.print("\n");
			}
			flag=false;
		}*/
		
		Sign[][] possibleSC = new Sign[num_sc][];
		int p_index = 0;
		
		float maxGain = Float.MIN_VALUE;
		float gain;
		
		Sign[] thisSC;
		int num_thissc;
		
		int num_dec = dec.length;		
		
		for (int i=0; i < num_sc; i++) {
			gain = 0;
			thisSC = signConfs[i];
			num_thissc = thisSC.length;
			
			// Summation of all gains in this SignConfig
			for (int j=0; j < num_thissc; j++) {
				for (int k=0; k < num_dec; k++) {
					if (dec[k].getTL() == thisSC[j]) {
						gain += dec[k].getGain();
					}
				}
			}
			
			// If the gain of this SignConfig is better than the Max till now...
			if (gain == maxGain) {
				possibleSC[p_index] = thisSC;
				p_index++;
			}
			else if (gain > maxGain) {
				possibleSC[0] = thisSC;
				p_index = 1;
				maxGain = gain;
			}
		}
		if(p_index == 0) {
			possibleSC[0] = signConfs[(int) Math.round(generator.nextFloat()*(num_sc-1))];
			p_index = 1;
		}
		
		int desSCId = (int) Math.round(generator.nextFloat()*(p_index-1));
		if(desSCId<0) {
			desSCId = 0;
			System.out.println("Dear sir.\n It seems your absolutely fabulous TrafficLightController algorithm caused quite a mess in this program. Could you be so kind to clean it up?\nThank you very much, sincerely\nThe Green Light District Team");
		}
		Sign[] desiredSC = possibleSC[desSCId];
		
//		System.out.print("Current SC is:");
//		for(int i=0; i<currentSC[node.getId()].length; i++)
//		{
//			System.out.print(currentSC[node.getId()][i].getLane().getId() + "-" + currentSC[node.getId()][i].getState() + " ");			
//		}
//		System.out.print("\n");
//		System.out.print("Desired SC is:");
//		for(int i=0; i<desiredSC.length; i++)
//		{			
//			System.out.print(desiredSC[i].getLane().getId() + "-" + desiredSC[i].getState() + " ");
//		}
//		System.out.print("\n");
		
		Sign[] nonOverlapSC = calcNonOverlapInFirst(currentSC[node.getId()], desiredSC);
		for(int i=0;i<nonOverlapSC.length;i++) {			
			nonOverlapSC[i].setCycleSwitchedToRed(tlcontroller.getCurCycle());
		}
		/*
		System.out.print("Nonoverlap SC in CurrentSC is:");
		for(int i=0; i<nonOverlapSC.length; i++)
		{			
			System.out.print(nonOverlapSC[i].getLane().getId() + " ");
		}
		System.out.print("\n");*/
		
		nonOverlapSC = calcNonOverlapInFirst(desiredSC,currentSC[node.getId()]);
		for(int i=0;i<nonOverlapSC.length;i++) {			
			nonOverlapSC[i].setCycleSwitchedToGreen(tlcontroller.getCurCycle());
		}
		/*
		System.out.print("Nonoverlap SC in DesiredSC is:");
		for(int i=0; i<nonOverlapSC.length; i++)
		{			
			System.out.print(nonOverlapSC[i].getLane().getId() + " ");
		}
		System.out.print("\n");*/
		
		if(CrossNodesSafely) {
			Sign[] overlapSC = calcOverlap(desiredSC, currentSC[node.getId()]);		
			
			int num_signs = signs.length;
			for (int i=0; i < num_signs; i++)
				signs[i].setState(false);
			if(true || node.areOtherTailsFree(overlapSC)) {
				currentSC[node.getId()] = desiredSC;
			}
			else {
				currentSC[node.getId()] = overlapSC;
			}
			
			int num_overlap = overlapSC.length;
			for(int i=0;i<num_overlap;i++) {
				overlapSC[i].setState(true);
				//overlapSC[i].setCycleSwitchedToRed(tlcontroller.getCurCycle());
			}
		}
	
		for(int i=0;i<num_dec;i++) {
			dec[i].getTL().setState(false);
		}
		
		int num_destl = possibleSC[desSCId].length;
		for(int i=0;i<num_destl;i++) {
			desiredSC[i].setState(true);
		}
	}
	
	protected Sign[] calcOverlap(Sign[] ar1, Sign[] ar2) {
		Sign[] outp = null;
		if(ar1 != null && ar2 != null) {
			int num_ar1 = ar1.length;
			int num_ar2 = ar2.length;
			outp = new Sign[(num_ar1<num_ar2)?num_ar1:num_ar2];
			int outp_index = 0;
			for(int i=0;i<num_ar1;i++) {
				for(int j=0;j<num_ar2;j++) {
					if(ar1[i].getId() == ar2[j].getId()) {
						outp[outp_index] = ar1[i];
						outp_index++;
					}
				}
			}
			outp = (Sign[]) Arrayutils.cropArray(outp,outp_index);
		}
		else {
			outp = new Sign[0];
		}		
		return outp;
	}
	
	// Used for setting cycleSwitched: Transition from green to red
	protected Sign[] calcNonOverlapInFirst(Sign[] ar1, Sign[] ar2) {
		Sign[] outp = null;
		if(ar1 != null && ar2 != null) {
			int num_ar1 = ar1.length;
			int num_ar2 = ar2.length;
			outp = new Sign[num_ar1];
			int outp_index = 0;
			for(int i=0;i<num_ar1;i++) {
				boolean flag=false;
				for(int j=0;j<num_ar2;j++) {
					if(ar1[i].getId() == ar2[j].getId()) { flag=true; break; }
				}
				if(!flag)
				{
					outp[outp_index] = ar1[i];
					outp_index++;
				}
			}
			outp = (Sign[]) Arrayutils.cropArray(outp,outp_index);
		}
		else {
			outp = new Sign[0];
		}		
		return outp;
	}
	
	public void load (XMLElement myElement,XMLLoader loader) throws XMLTreeException,IOException,XMLInvalidInputException
	{	currentSC=(Sign[][])XMLArray.loadArray(this,loader);
	}
	
	public XMLElement saveSelf () throws XMLCannotSaveException
	{ 	XMLElement result=new XMLElement(shortXMLName);
		result.setName(shortXMLName);	
	  	return result;
	}
	
	public void saveChilds (XMLSaver saver) throws XMLTreeException,IOException,XMLCannotSaveException
	{ 	for ( int x=0; x< currentSC.length ; x++)
		{	for ( int y=0;y < currentSC[x].length ; y++ )
			{	currentSC[x][y].setParentName(getXMLName());
			}
		}
		XMLArray.saveArray(currentSC,this,saver,"current-sc");
	}
	
	public String getXMLName ()
	{ 	return "model."+shortXMLName;
	}
	
	public void setParentName (String parentName) throws XMLTreeException
	{	throw new XMLTreeException
		("Attempt to change fixed parentName of a SignController class.");
	}
	
	public void callUpdateQ()
	{
		tlcontroller.updateQ();
	}
	
	public void callFinish()
	{
		tlcontroller.finish();
	}
	
}