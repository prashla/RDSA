package gld.tt;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import gld.*;
import gld.algo.tlc.Key;
import gld.infra.EdgeNode;
import gld.infra.Infrastructure;
import gld.infra.Junction;
import gld.infra.Node;
import gld.infra.SpecialNode;
import gld.sim.*;
import gld.sim.stats.*;
import gld.utils.ListCombinationGenerator;

public class GLDRunner
{
    String tlController, mapdir, map;
    int tlcID, tlcCat;
    String[] params;

    GLDStarter gldStarter;
    SimController simController;
    TotalWaitTrackingView viewATWT;
    TotalRoadusersTrackingView viewTAR;

    GLDRunner(String[] params, String mapdir, String map, String tlController,
            int tlcCat, int tlcID)
    {
        this.params = params;
        this.mapdir = mapdir;
        this.map = map;
        this.tlController = tlController;
        this.tlcCat = tlcCat;
        this.tlcID = tlcID;
    }


    public void init()
    {
        gldStarter = new GLDStarter(params, GLDStarter.SIMULATOR);
        simController = (SimController) gldStarter.getController();
        simController.setTLC(tlcCat, tlcID); // RS_SPSA

        try
        {
            simController.doLoad(mapdir + map + ".sim");
        }
        catch (Exception e)
        {
        }

        Infrastructure infra = simController.getSimModel().getInfrastructure();
        infra.setTitle(mapdir + map + ".sim");

        viewATWT = new TotalWaitTrackingView(
                simController.getSimModel().getCurCycle(),
                simController.getSimModel());
        TrackerFactory.genExtTracker(simController.getSimModel(), simController,
                viewATWT);
        viewTAR = new TotalRoadusersTrackingView(
                simController.getSimModel().getCurCycle(),
                simController.getSimModel());
        TrackerFactory.genExtTracker(simController.getSimModel(), simController,
                viewTAR);
    }

    public void resetSimulator()
    {
        try
        {
            simController.getSimModel().reset();
        }
        catch (SimulationRunningException e)
        {
            e.printStackTrace();
        }
    }
    
    public double run(double[] theta, int simDuration, boolean getTrafficStats, String filePrefix)
    {
        simController.setSpeed(3); // speed to maximum
        simController.setTLC(tlcCat, tlcID); // RS_SPSA

        simController.getSimModel().setSimName("CPT SPSA");
        simController.getSimModel().getTLController().setThresholds(6, 14, 130);

        simController.getSimModel().getTLController()
                .setLanesFile(mapdir + map + "_pri.txt");
        simController.getSimModel().getTLController().setTheta(theta);

        simController.unpause(); 
        simController.getSimModel().setSimulationDuration(simDuration);

        while (simController.getSimModel().getCurCycle() < simController
                .getSimModel().getSimulationDuration())
        {
            try
            {
                Thread.sleep(10000);
                System.out.println(
                        "CYCLE:" + simController.getSimModel().getCurCycle());
            }
            catch (InterruptedException e)
            {

            }
        }
        
        if(getTrafficStats == true)
        {
            try
            {
                // simController.getSimModel().getTLController().finish();
                viewATWT.saveData(
                        filePrefix + "_ATWT_" + map  + ".txt",
                        simController.getSimModel());
                viewTAR.saveData(
                        filePrefix + "_TAR_" + map + ".txt",
                        simController.getSimModel());
            }
            catch (Exception e)
            {
            }
        }

        return simController.getSimModel().getTLController().getAverageCost();
    }
}
