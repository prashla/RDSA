
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

import java.awt.Point;
import java.awt.Dimension;
import gld.infra.*;

/**
 *
 * Simple infrastructure used for testing
 *
 * @author Group Datastructures
 * @version 1.0
 */

public class LessSimpleInfra extends Infrastructure
{
	public LessSimpleInfra()
	{
		super(new Dimension(1000, 800));
		
		try {
			EdgeNode edge0 = new EdgeNode(new Point(-200,-200));
			EdgeNode edge1 = new EdgeNode(new Point(100,-200));
			EdgeNode edge2 = new EdgeNode(new Point(-200,0));
			EdgeNode edge3 = new EdgeNode(new Point(100,0));
			Junction node4 = new Junction(new Point(-100,-200));
			Junction node5 = new Junction(new Point(0,-200));
			Junction node6 = new Junction(new Point(-100,0));
			Junction node7 = new Junction(new Point(0,0));
			Junction node8 = new Junction(new Point(-100,-100));
			Junction node9 = new Junction(new Point(0,-100));
		
			edge0.setId(0);
			edge1.setId(1);
			edge2.setId(2);
			edge3.setId(3);
			node4.setId(4);
			node5.setId(5);
			node6.setId(6);
			node7.setId(7);
			node8.setId(8);
			node9.setId(9);

			Road road1 = new Road(node4, edge0, 10);
			Road road2 = new Road(node5, node4, 10);
			Road road3 = new Road(edge1, node5, 10);
			Road road4 = new Road(node8, node4, 10);
			Road road5 = new Road(node9, node5, 10);
			Road road6 = new Road(node6, node8, 10);
			Road road7 = new Road(node7, node9, 10);
			Road road8 = new Road(node6, edge2, 10);
			Road road9 = new Road(node7, node6, 10);
			Road road10 = new Road(edge3, node7, 10);
			Road road11 = new Road(node9, node8, 10);
		
			Drivelane d11 = new Drivelane(road1);
			d11.setType(1);
			boolean[] d11targets = {false, true, true};
			d11.setTargets(d11targets);
			TrafficLight s11 = new TrafficLight(node4,d11);
			d11.setSign(s11);
		
			Drivelane d12 = new Drivelane(road1);
			d12.setType(1);
			boolean[] d12targets = {false, true, false};
			d12.setTargets(d12targets);
			TrafficLight s12 = new TrafficLight(node4,d12);
			d12.setSign(s12);

			Drivelane d13 = new Drivelane(road1);
			d13.setType(1);
			boolean[] d13targets = {false, true, false};
			d13.setTargets(d13targets);
			TrafficLight s13 = new TrafficLight(edge0,d13);
			d13.setSign(s13);

			Drivelane d14 = new Drivelane(road1);
			d14.setType(1);
			boolean[] d14targets = {false, true, false};
			d14.setTargets(d14targets);
			TrafficLight s14 = new TrafficLight(edge0,d14);
			d14.setSign(s14);
		
		
			Drivelane d21 = new Drivelane(road2);
			d21.setType(1);
			boolean[] d21targets = {false, true, true};
			d21.setTargets(d21targets);
			TrafficLight s21 = new TrafficLight(node5,d21);
			d21.setSign(s21);
		
			Drivelane d22 = new Drivelane(road2);
			d22.setType(1);
			boolean[] d22targets = {false, true, false};
			d22.setTargets(d22targets);
			TrafficLight s22 = new TrafficLight(node5,d22);
			d22.setSign(s22);		
		
			Drivelane d23 = new Drivelane(road2);
			d23.setType(1);
			boolean[] d23targets = {true, false, false};
			d23.setTargets(d23targets);
			TrafficLight s23 = new TrafficLight(node4,d23);
			d23.setSign(s23);
		
			Drivelane d24 = new Drivelane(road2);
			d24.setType(1);
			boolean[] d24targets = {false, true, false};
			d24.setTargets(d24targets);
			TrafficLight s24 = new TrafficLight(node4,d24);
			d24.setSign(s24);		

		
			Drivelane d31 = new Drivelane(road3);
			d31.setType(1);
			boolean[] d31targets = {false, true, false};
			d31.setTargets(d31targets);
			TrafficLight s31 = new TrafficLight(edge1,d31);
			d31.setSign(s31);
		
			Drivelane d32 = new Drivelane(road3);
			d32.setType(1);
			boolean[] d32targets = {false, true, false};
			d32.setTargets(d32targets);
			TrafficLight s32 = new TrafficLight(edge1,d32);
			d32.setSign(s32);		

			Drivelane d33 = new Drivelane(road3);
			d33.setType(1);
			boolean[] d33targets = {true, false, false};
			d33.setTargets(d33targets);
			TrafficLight s33 = new TrafficLight(node5,d33);
			d33.setSign(s33);
		
			Drivelane d34 = new Drivelane(road3);
			d34.setType(1);
			boolean[] d34targets = {false, true, false};
			d34.setTargets(d34targets);
			TrafficLight s34 = new TrafficLight(node5,d34);
			d34.setSign(s34);		


			Drivelane d41 = new Drivelane(road4);
			d41.setType(1);
			boolean[] d41targets = {false, true, false};
			d41.setTargets(d41targets);
			TrafficLight s41 = new TrafficLight(node8, d41);
			d41.setSign(s41);
		
			Drivelane d42 = new Drivelane(road4);
			d42.setType(1);
			boolean[] d42targets = {true, false, false};
			d42.setTargets(d42targets);
			TrafficLight s42 = new TrafficLight(node8, d42);
			d42.setSign(s42);

			Drivelane d43 = new Drivelane(road4);
			d43.setType(1);
			boolean[] d43targets = {true, false, false};
			d43.setTargets(d43targets);
			TrafficLight s43 = new TrafficLight(node4, d43);
			d43.setSign(s43);
		
			Drivelane d44 = new Drivelane(road4);
			d44.setType(1);
			boolean[] d44targets = {false, false, true};
			d44.setTargets(d44targets);
			TrafficLight s44 = new TrafficLight(node4, d44);
			d44.setSign(s44);


			Drivelane d51 = new Drivelane(road5);
			d51.setType(1);
			boolean[] d51targets = {false, true, true};
			d51.setTargets(d51targets);
			TrafficLight s51 = new TrafficLight(node9, d51);
			d51.setSign(s51);
		
			Drivelane d52 = new Drivelane(road5);
			d52.setType(1);
			boolean[] d52targets = {false, true, false};
			d52.setTargets(d52targets);
			TrafficLight s52 = new TrafficLight(node9, d52);
			d52.setSign(s52);		

			Drivelane d53 = new Drivelane(road5);
			d53.setType(1);
			boolean[] d53targets = {true, false, false};
			d53.setTargets(d53targets);
			TrafficLight s53 = new TrafficLight(node5, d53);
			d53.setSign(s53);
		
			Drivelane d54 = new Drivelane(road5);
			d54.setType(1);
			boolean[] d54targets = {false, false, true};
			d54.setTargets(d54targets);
			TrafficLight s54 = new TrafficLight(node5, d54);
			d54.setSign(s54);		
	
		
			Drivelane d61 = new Drivelane(road6);
			d61.setType(1);
			boolean[] d61targets = {false, false, true};
			d61.setTargets(d61targets);
			TrafficLight s61 = new TrafficLight(node6, d61);
			d61.setSign(s61);
		
			Drivelane d62 = new Drivelane(road6);
			d62.setType(1);
			boolean[] d62targets = {true, false, false};
			d62.setTargets(d62targets);
			TrafficLight s62 = new TrafficLight(node6, d62);
			d62.setSign(s62);		

			Drivelane d63 = new Drivelane(road6);
			d63.setType(1);
			boolean[] d63targets = {false, true, false};
			d63.setTargets(d63targets);
			TrafficLight s63 = new TrafficLight(node8, d63);
			d63.setSign(s63);
		
			Drivelane d64 = new Drivelane(road6);
			d64.setType(1);
			boolean[] d64targets = {false, true, true};
			d64.setTargets(d64targets);
			TrafficLight s64 = new TrafficLight(node8, d64);
			d64.setSign(s64);
		
		
			Drivelane d71 = new Drivelane(road7);
			d71.setType(1);
			boolean[] d71targets = {false, false, true};
			d71.setTargets(d71targets);
			TrafficLight s71 = new TrafficLight(node7,d71);
			d71.setSign(s71);
		
			Drivelane d72 = new Drivelane(road7);
			d72.setType(1);
			boolean[] d72targets = {true, false, false};
			d72.setTargets(d72targets);
			TrafficLight s72 = new TrafficLight(node7,d72);
			d72.setSign(s72);		

			Drivelane d73 = new Drivelane(road7);
			d73.setType(1);
			boolean[] d73targets = {true, false, false};
			d73.setTargets(d73targets);
			TrafficLight s73 = new TrafficLight(node9,d73);
			d73.setSign(s73);
		
			Drivelane d74 = new Drivelane(road7);
			d74.setType(1);
			boolean[] d74targets = {false, true, false};
			d74.setTargets(d74targets);
			TrafficLight s74 = new TrafficLight(node9,d74);
			d74.setSign(s74);
		
		
			Drivelane d81 = new Drivelane(road8);
			d81.setType(1);
			boolean[] d81targets = {false, true, false};
			d81.setTargets(d81targets);
			TrafficLight s81 = new TrafficLight(node6,d81);
			d81.setSign(s81);
		
			Drivelane d82 = new Drivelane(road8);
			d82.setType(1);
			boolean[] d82targets = {true, false, false};
			d82.setTargets(d82targets);
			TrafficLight s82 = new TrafficLight(node6,d82);
			d82.setSign(s82);		

			Drivelane d83 = new Drivelane(road8);
			d83.setType(1);
			boolean[] d83targets = {false, true, false};
			d83.setTargets(d83targets);
			TrafficLight s83 = new TrafficLight(edge2,d83);
			d83.setSign(s83);
		
			Drivelane d84 = new Drivelane(road8);
			d84.setType(1);
			boolean[] d84targets = {false, true, false};
			d84.setTargets(d84targets);
			TrafficLight s84 = new TrafficLight(edge2,d84);
			d84.setSign(s84);	
		
		
			Drivelane d91 = new Drivelane(road9);
			d91.setType(1);
			boolean[] d91targets = {false, true, false};
			d91.setTargets(d91targets);
			TrafficLight s91 = new TrafficLight(node7,d91);
			d91.setSign(s91);
		
			Drivelane d92 = new Drivelane(road9);
			d92.setType(1);
			boolean[] d92targets = {true, false, false};
			d92.setTargets(d92targets);
			TrafficLight s92 = new TrafficLight(node7,d92);
			d92.setSign(s92);		

			Drivelane d93 = new Drivelane(road9);
			d93.setType(1);
			boolean[] d93targets = {false, true, false};
			d93.setTargets(d93targets);
			TrafficLight s93 = new TrafficLight(node6,d93);
			d93.setSign(s93);
		
			Drivelane d94 = new Drivelane(road9);
			d94.setType(1);
			boolean[] d94targets = {false, true, true};
			d94.setTargets(d94targets);
			TrafficLight s94 = new TrafficLight(node6,d94);
			d94.setSign(s94);
		
		
			Drivelane d101 = new Drivelane(road10);
			d101.setType(1);
			boolean[] d101targets = {false, true, false};
			d101.setTargets(d101targets);
			TrafficLight s101 = new TrafficLight(edge3,d101);
			d101.setSign(s101);
		
			Drivelane d102 = new Drivelane(road10);
			d102.setType(1);
			boolean[] d102targets = {false, true, false};
			d102.setTargets(d102targets);
			TrafficLight s102 = new TrafficLight(edge3,d102);
			d102.setSign(s102);		

			Drivelane d103 = new Drivelane(road10);
			d103.setType(1);
			boolean[] d103targets = {false, true, false};
			d103.setTargets(d103targets);
			TrafficLight s103 = new TrafficLight(node7,d103);
			d103.setSign(s103);
		
			Drivelane d104 = new Drivelane(road10);
			d104.setType(1);
			boolean[] d104targets = {false, true, true};
			d104.setTargets(d104targets);
			TrafficLight s104 = new TrafficLight(node7,d104);
			d104.setSign(s104);	
		
		
			Drivelane d111 = new Drivelane(road11);
			d111.setType(1);
			boolean[] d111targets = {false, false, true};
			d111.setTargets(d111targets);
			TrafficLight s111 = new TrafficLight(node9,d111);
			d111.setSign(s111);
		
			Drivelane d112 = new Drivelane(road11);
			d112.setType(1);
			boolean[] d112targets = {true, false, false};
			d112.setTargets(d112targets);
			TrafficLight s112 = new TrafficLight(node9,d112);
			d112.setSign(s112);		

			Drivelane d113 = new Drivelane(road11);
			d113.setType(1);
			boolean[] d113targets = {true, false, false};
			d113.setTargets(d113targets);
			TrafficLight s113 = new TrafficLight(node8,d113);
			d113.setSign(s113);
		
			Drivelane d114 = new Drivelane(road11);
			d114.setType(1);
			boolean[] d114targets = {false, false, true};
			d114.setTargets(d114targets);
			TrafficLight s114 = new TrafficLight(node8,d114);
			d114.setSign(s114);



			road1.addAlphaLane(d11);
			road1.addAlphaLane(d12);
			road1.addBetaLane(d13);
			road1.addBetaLane(d14);
			road2.addAlphaLane(d21);
			road2.addAlphaLane(d22);
			road2.addBetaLane(d23);
			road2.addBetaLane(d24);
			road3.addAlphaLane(d31);
			road3.addAlphaLane(d32);
			road3.addBetaLane(d33);
			road3.addBetaLane(d34);
			road4.addAlphaLane(d41);
			road4.addAlphaLane(d42);
			road4.addBetaLane(d43);
			road4.addBetaLane(d44);
			road5.addAlphaLane(d51);
			road5.addAlphaLane(d52);
			road5.addBetaLane(d53);
			road5.addBetaLane(d54);
			road6.addAlphaLane(d61);
			road6.addAlphaLane(d62);
			road6.addBetaLane(d63);
			road6.addBetaLane(d64);
			road7.addAlphaLane(d71);
			road7.addAlphaLane(d72);
			road7.addBetaLane(d73);
			road7.addBetaLane(d74);
			road8.addAlphaLane(d81);
			road8.addAlphaLane(d82);
			road8.addBetaLane(d83);
			road8.addBetaLane(d84);
			road9.addAlphaLane(d91);
			road9.addAlphaLane(d92);
			road9.addBetaLane(d93);
			road9.addBetaLane(d94);	
			road10.addAlphaLane(d101);
			road10.addAlphaLane(d102);
			road10.addBetaLane(d103);
			road10.addBetaLane(d104);
			road11.addAlphaLane(d111);
			road11.addAlphaLane(d112);
			road11.addBetaLane(d113);
			road11.addBetaLane(d114);
		
			edge0.addRoad(road1, 1);
			edge1.addRoad(road3, 3);
			edge2.addRoad(road8, 1);
			edge3.addRoad(road10, 3);
			node4.addRoad(road2,1);
			node4.addRoad(road4,2);
			node4.addRoad(road1,3);
			node5.addRoad(road3,1);
			node5.addRoad(road5,2);
			node5.addRoad(road2,3);
			node6.addRoad(road6,0);
			node6.addRoad(road9,1);
			node6.addRoad(road8,3);
			node7.addRoad(road7,0);
			node7.addRoad(road10,1);
			node7.addRoad(road9,3);
			node8.addRoad(road4,0);
			node8.addRoad(road11,1);
			node8.addRoad(road6,2);
			node9.addRoad(road5,0);
			node9.addRoad(road7,2);
			node9.addRoad(road11,3);
		
			Sign[] signs4 = { s11, s12, s23, s24, s43, s44 };
			node4.setSigns(signs4);
			
			Sign[] signs5 = { s21, s22, s33, s34, s53, s54 };
			node5.setSigns(signs5);
			
			Sign[] signs6 = { s61, s62, s81, s82, s93, s94 };
			node6.setSigns(signs6);
			
			Sign[] signs7 = { s71, s72, s91, s92, s103, s104 };
			node7.setSigns(signs7);
			
			Sign[] signs8 = { s41, s42, s63, s64, s113, s114 };
			node8.setSigns(signs8);
			
			Sign[] signs9 = { s51, s52, s73, s74, s111, s112 };
			node9.setSigns(signs9);


			Node[] nodes = {edge0,edge1,edge2,edge3,node4,node5,node6,node7,node8,node9};	
			EdgeNode[] edges = {edge0,edge1,edge2,edge3};
			Dimension d = new Dimension(1000,800);
			
			allNodes = nodes;
			specialNodes = edges;

	
			{SpawnFrequency sf = new SpawnFrequency(1,(float) 0.2);
			SpawnFrequency[] sfa = {sf};
			edge1.setSpawnFrequencies(sfa);}
	
			{SpawnFrequency sf = new SpawnFrequency(1,(float) 0.1);
			SpawnFrequency[] sfa = {sf};
			edge2.setSpawnFrequencies(sfa);}
	
			{SpawnFrequency sf = new SpawnFrequency(1,(float) 0.3);
			SpawnFrequency[] sfa = {sf};
			edge3.setSpawnFrequencies(sfa);}
			
		}
		catch (InfraException e) {
			e.printStackTrace();
		}
	}
}
			