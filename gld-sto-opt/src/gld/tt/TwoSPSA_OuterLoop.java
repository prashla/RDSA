package gld.tt;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;

import gld.infra.EdgeNode;
import gld.infra.Infrastructure;
import gld.infra.Junction;
import gld.infra.Node;
import gld.infra.SpecialNode;
import gld.sim.*;
import gld.sim.stats.*;
import gld.utils.ListCombinationGenerator;
import Jama.Matrix;

/**
 * Thresholds parameter optimizing outer loop that is based on 2RDSA
 * 
 * @see RDSARunner, TwoRDSA_OuterLoop
 * @version 1.0 31 Jul 2016
 * @author Prashanth L.A.
 */
public class TwoSPSA_OuterLoop
{
    /**
     * Training phase executes PTLC and this class is a wrapper that takes in a
     * threshold parameter, runs the GLD simulator with PTLC for traffic lights
     * and returns the average cost
     */
    GLDRunner gldRunnerObject;

    int dimTheta;
    PrintWriter thetaFile, traceFile, resultsFile;

    // See the constructor for a description of these quantities
    int trainingBudget, numReplicationsForTesting,
            trajectoryLengthTrainingPhase, trajectoryLengthTestingPhase;

    // Projection limits for the thresholds while they are tuned in training
    // phase
    double THETA_MIN[] = { 0, 8, 16, 24, 32, 40, 60, 120 };
    double THETA_MAX[] = { 8, 16, 24, 32, 40, 60, 120, 180 };

    Random randGen;

    // For projecting Hessian onto the set of positive definite matrices, each
    // diagonal entry is projected onto the positive side using the variable
    // below
    double eta = 0.15;

    // Random perturbations
    int Delta[], DeltaCap[];

    // deterministic perturbation constants
    double delta;

    // Step sizes
    double cn = 1.0, an = 1.0;

    // Cost function
    double Zplus, Zminus, Zplusplus, Zminusplus;

    // Threshold parameter
    double theta[];
    // Perturbed threshold parameters
    double thetaPlusDelta[], thetaMinusDelta[], thetaPlusDeltaDeltacap[],
    thetaMinusDeltaPlusDeltacap[];
    
    // Hessian estimate
    double HessianEstimate[][];

    /**
     * Initialize the PTLC+2RDSA algorithm
     * 
     * @param mapdir
     *            +map Specify the path to the map file containing the road
     *            network to be tested
     * @param trainingBudget
     *            the number of function evaluations in the training phase
     * @param numReplicationsForTesting
     *            After the training phase, the thresholds are fixed and then a
     *            number of independent simulations (that this variable is set
     *            to) are run and the empirical average cost from each
     *            simulation is recorded
     * @param trajectoryLengthTrainingPhase
     *            length of each simulated trajectory during training phase
     * @param trajectoryLengthTestingPhase
     *            length of each simulated trajectory during testing phase
     */
    TwoSPSA_OuterLoop(String[] params, String mapdir, String map,
            int trainingBudget, int numReplicationsForTesting,
            int trajectoryLengthTrainingPhase,
            int trajectoryLengthTestingPhase)
    {
        System.out
                .println("Initialize TwoSPSA");

        this.trainingBudget = trainingBudget;
        this.numReplicationsForTesting = numReplicationsForTesting;
        this.trajectoryLengthTrainingPhase = trajectoryLengthTrainingPhase;
        this.trajectoryLengthTestingPhase = trajectoryLengthTestingPhase;

        String tlController = "PTLC_2SPSA";
        int tlControllerId = 0;
        int tlcCat = 0;

        gldRunnerObject = new GLDRunner(params, mapdir, map, tlController,
                tlcCat, tlControllerId);
        gldRunnerObject.init();

        randGen = new Random(1444);

        // traceFile.println("numRealSigns: " + numRealSigns);
        dimTheta = 8;

        // Initialize phi vars
        theta = new double[dimTheta];
        theta[0] = 8.0;
        theta[1] = 16.0;
        theta[2] = 24;
        theta[3] = 32.0;
        theta[4] = 40.0;
        theta[5] = 60;
        theta[6] = 120;
        theta[7] = 180;

        // Hessian estimate matrix
        HessianEstimate = new double[dimTheta][dimTheta];
        for (int i = 0; i < HessianEstimate.length; i++)
        {
            HessianEstimate[i][i] = 1;
        }

        // Initialize phi vars
        thetaPlusDelta = new double[dimTheta];
        thetaMinusDelta = new double[dimTheta];
        thetaPlusDeltaDeltacap = new double[dimTheta];
        thetaMinusDeltaPlusDeltacap = new double[dimTheta];
        Delta = new int[dimTheta];
        DeltaCap = new int[dimTheta];


        delta = 3.8;
        Zplus = Zminus = Zplusplus = Zminusplus = 0;

    }

    public void setPerturbedThetas()
    {
        setDeltas();

        // Perturb theta now
        for (int i = 0; i < theta.length; i++)
        {
            thetaPlusDelta[i] = theta[i] + delta * Delta[i];
            thetaMinusDelta[i] = theta[i] - delta * Delta[i];
            thetaPlusDeltaDeltacap[i] = theta[i] + delta * Delta[i]
                    + delta * DeltaCap[i];
            thetaMinusDeltaPlusDeltacap[i] = theta[i] - delta * Delta[i]
                    + delta * DeltaCap[i];
        }
    }

    void setDeltas()
    {
        for (int i = 0; i < Delta.length; i++)
        {
            Delta[i] = genBernoulli();
            DeltaCap[i] = genBernoulli();
        }
    }

    int genBernoulli()
    {
        int rv;
        double rndValue = randGen.nextDouble();
        if (rndValue < 0.5)
            rv = 1;
        else
            rv = -1;
        return rv;
    }

    double[] getTheta(int thetaTypeforSimulation)
    {
        if (thetaTypeforSimulation == 0)
        {
            System.out.println("Type: " + thetaTypeforSimulation + " theta: "
                    + Arrays.toString(thetaPlusDelta));
            return thetaPlusDelta;
        }
        else if (thetaTypeforSimulation == 1)
        {
            traceFile.println("Type: " + thetaTypeforSimulation
                    + " thetaPlusDelta: " + Arrays.toString(thetaMinusDelta));
            return thetaMinusDelta;
        }
        else if (thetaTypeforSimulation == 2)
        {
            traceFile.println("Type: " + thetaTypeforSimulation
                    + " thetaMinusDeltaPlusDeltacap: "
                    + Arrays.toString(thetaMinusDeltaPlusDeltacap));
            return thetaMinusDeltaPlusDeltacap;
        }
        else
        {
            traceFile.println(
                    "Type: " + thetaTypeforSimulation + " thetaMinusDelta: "
                            + Arrays.toString(thetaPlusDeltaDeltacap));
            return thetaPlusDeltaDeltacap;
        }
    }

    public void initFiles(String netw)
    {
        try
        {
            thetaFile = new PrintWriter(new FileWriter(new File(
                    "2SPSA_theta_" + netw + ".txt")));
            thetaFile.println("TLC Algo: PTLC");
            thetaFile.println("# theta[1] .... theta[dimFeatures]");
            thetaFile.println("#");
            // trace log
            traceFile = new PrintWriter(new FileWriter(new File(
                    "2SPSA_trace_" + netw + ".log")));
            resultsFile = new PrintWriter(new FileWriter(new File(
                    "2SPSA_results_" + netw + ".log")));
        }
        catch (Exception e)
        {
            System.out.println("PTLC_2SPSA failed to initialise.");
        }
    }

    public void finish()
    {
        thetaFile.close();
        traceFile.close();
        resultsFile.close();
    }

    void updateParams(int simNr)
    {
        for (int i = 0; i < HessianEstimate.length; i++)
        {
            double zplusfactor = (Zplusplus - Zplus) / (delta * DeltaCap[i]);
            double zminusfactor = (Zminusplus - Zminus) / (delta * DeltaCap[i]);
            HessianEstimate[i][i] = HessianEstimate[i][i] + cn * ((zplusfactor - zminusfactor)
                    / (2 * delta * delta * Delta[i]) - HessianEstimate[i][i]);

            if (HessianEstimate[i][i] <= eta)
            {
                HessianEstimate[i][i] = eta;
            }
        }

        Matrix Pij = new Matrix(HessianEstimate);
        Pij = Pij.inverse();
        // Upper bound Pkl
        for (int i = 0; i < HessianEstimate.length; i++)
        {
            if (Pij.get(i, i) < eta)
            {
                Pij.set(i, i, eta);
            }
            if (Pij.get(i, i) > 1 / eta)
            {
                Pij.set(i, i, 1 / eta);
            }
        }
        traceFile.println("Pij(n): " + Arrays.deepToString(Pij.getArray()));

        for (int i = 0; i < theta.length; i++)
        {
            double prod = 0;
            for (int j = 0; j < HessianEstimate.length; j++)
            {
                prod += Pij.get(i, j) * (Zminus - Zplus)
                        / (2 * delta * Delta[j]);
            }
            theta[i] -= an * prod;

            if (theta[i] <= THETA_MIN[i])
                theta[i] = THETA_MIN[i];
            if (theta[i] >= THETA_MAX[i])
                theta[i] = THETA_MAX[i];

        }

        cn = Math.pow(simNr, -1); // Hessian medium ts
        an = Math.pow(simNr, -0.6); // theta slowest ts
        delta = Math.pow(simNr, -0.101);

        System.out.println("theta: " + Arrays.toString(theta));
        traceFile.println("theta: " + Arrays.toString(theta));
        String ln = simNr + "\t";
        thetaFile.print(ln);
        for (int i = 0; i < theta.length; i++)
        {
            // Record theta
            thetaFile.printf("%.2f\t", theta[i]);
        }
        thetaFile.println();
    }

    /**
     * This constitutes the training phase where the thresholds are tuned
     */
    void optimizeThresholds()
    {
        traceFile.println("############# Training 2SPSA ################");
        System.out.println("############# Training 2SPSA ################");

        for (int simNr = 1; simNr < trainingBudget + 1; simNr++)
        {
            traceFile.println("############################ " + simNr
                    + " ################################");
            System.out.println("############################ " + simNr
                    + " ################################");

            int thetaTypeforSimulation = simNr % 4;
            if (thetaTypeforSimulation == 1)
            {
                setPerturbedThetas();
            }
            double thetaSim[];
            thetaSim = getTheta(thetaTypeforSimulation);
            double avgCostSample = gldRunnerObject.run(thetaSim,
                    trajectoryLengthTrainingPhase, false, "");

            gldRunnerObject.resetSimulator();

            if (simNr % 4 == 0)
            {
                Zplusplus = avgCostSample;
                traceFile.println("Zplusplus:" + Zplusplus);
                    updateParams(simNr);
            }
            else if (simNr % 4 == 3)
            {
                Zminusplus = avgCostSample;
                traceFile.println("Zminusplus:" + Zminusplus);
            }
            else if (simNr % 4 == 2)
            {
                Zminus = avgCostSample;
                traceFile.println("Zminus:" + Zminus);
            }
            else if (simNr % 4 == 1)
            {
                Zplus = avgCostSample;
                traceFile.println("Zplus:" + Zplus);
            }
        }
    }

    /**
     * This constitutes the test phase where the thresholds from the end of
     * training phase are taken, multiple simulations are performed and the
     * empirical average cost is recorded
     */
    void testThresholds()
    {
        traceFile.println("############# Testing 2SPSA ################");
        System.out.println("############# Testing 2SPSA ################");

        for (int simNr = 1; simNr < numReplicationsForTesting + 1; simNr++)
        {
            traceFile.println("############################ " + simNr
                    + " ################################");
            System.out.println("############################ " + simNr
                    + " ################################");

            String filePrefix = "2SPSA_" + Integer.toString(simNr);
            double avgCostSample = gldRunnerObject.run(theta,
                    trajectoryLengthTestingPhase, true, filePrefix);

            resultsFile.println(avgCostSample);
            System.out.println("Avg Value: " + avgCostSample);

            gldRunnerObject.resetSimulator();
        }
        System.out.println("******************************************************************************");
    }
}
