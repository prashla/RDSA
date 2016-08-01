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
 * @see RDSARunner, TwoSPSA_OuterLoop
 * @version 1.0 31 Jul 2016
 * @author Prashanth L.A.
 */
public class TwoRDSA_OuterLoop
{
    /**
     * The different performance criterion: (i) AVG -> no gain/loss separation
     * via utilities, no probability distortion, (ii) EUT -> has gain/loss
     * separation via utilities, no probability distortion, (iii) CPT -> has
     * both gain/loss separation and probability distortion
     */
    protected static final int UNIF = 0, ASYMBER = 1;
    protected static final String[] typeDescs = { "2RDSA_Unif", "2RDSA_AsymBer" };
    int perturbationType;

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

    // Asym Bernoulli distribution parameter
    double epsilon = 1;

    // Random perturbations
    double Delta[];

    // deterministic perturbation constants
    double delta;

    // Step sizes
    double cn = 1.0, an = 1.0;

    // Cost function samples corresponding to three simulations in each
    // iteration during training phase
    double Zplus, Zminus, Z;

    // Threshold parameter
    double theta[];
    // Perturbed threshold parameters
    double thetaPlusDelta[], thetaMinusDelta[];

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
     * @param perturbationType
     *            UNIF for uniform perturbations, ASYMBER for asymmetric Bernoulli  
     */
    TwoRDSA_OuterLoop(String[] params, String mapdir, String map,
            int trainingBudget, int numReplicationsForTesting,
            int trajectoryLengthTrainingPhase,
            int trajectoryLengthTestingPhase, int perturbationType)
    {

        this.trainingBudget = trainingBudget;
        this.numReplicationsForTesting = numReplicationsForTesting;
        this.trajectoryLengthTrainingPhase = trajectoryLengthTrainingPhase;
        this.trajectoryLengthTestingPhase = trajectoryLengthTestingPhase;
        this.perturbationType = perturbationType;

        String tlController = "PTLC_" + typeDescs[perturbationType];
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
        Delta = new double[dimTheta];

        delta = 3.8;
        Zplus = Zminus = Z = 0;

        System.out.println("Initialized " + typeDescs[perturbationType]);

    }

    public void setPerturbedThetas()
    {
        setDeltas();

        // Perturb theta now
        for (int i = 0; i < theta.length; i++)
        {
            // Perturb thetas - will be used in the system
            thetaPlusDelta[i] = theta[i] + delta * Delta[i];
            thetaMinusDelta[i] = theta[i] - delta * Delta[i];
        }
    }

    void setDeltas()
    {
        for (int i = 0; i < Delta.length; i++)
        {
            Delta[i] = genRandomPerturbations();
        }
    }

    double genRandomPerturbations()
    {
        double rv = 0.0;
        double rndValue = randGen.nextDouble();

        switch (perturbationType)
        {
        case UNIF:
            rv = -1.0 + rndValue * 2.0;
            break;

        case ASYMBER:
            if (rndValue < (1 + epsilon) / (2 + epsilon))
                rv = -1;
            else
                rv = 1 + epsilon;
            break;
        }

        return rv;
    }

    /**
     * Return a perturbed thresholds parameter
     * 
     * @param thetaTypeforSimulation
     *            0, 1, 2 correspond to theta, (theta+deltaDelta) and (theta-deltaDelta)
     */
    double[] getTheta(int thetaTypeforSimulation)
    {
        if (thetaTypeforSimulation == 0)
        {
            System.out.println("Type: " + thetaTypeforSimulation + " theta: "
                    + Arrays.toString(theta));
            return theta;
        }
        else if (thetaTypeforSimulation == 1)
        {
            traceFile.println("Type: " + thetaTypeforSimulation
                    + " thetaPlusDelta: " + Arrays.toString(thetaPlusDelta));
            return thetaPlusDelta;
        }
        else
        {
            traceFile.println("Type: " + thetaTypeforSimulation
                    + " thetaMinusDelta: " + Arrays.toString(thetaMinusDelta));
            return thetaMinusDelta;
        }
    }

    /**
     * Open files for trace prints, theta in each iteration and results
     * (average cost-value) from testing phase
     * 
     * @param netw
     *            Specify the path to the map file containing the road network
     *            to be tested
     */
    public void initFiles(String netw)
    {
        try
        {
            thetaFile = new PrintWriter(new FileWriter(new File(
                    typeDescs[perturbationType] + "_theta_" + netw + ".txt")));
            thetaFile.println("TLC Algo: PTLC");
            thetaFile.println("# theta[1] .... theta[dimFeatures]");
            thetaFile.println("#");
            // trace log
            traceFile = new PrintWriter(new FileWriter(new File(
                    typeDescs[perturbationType] + "_trace_" + netw + ".log")));
            resultsFile = new PrintWriter(new FileWriter(new File(
                    typeDescs[perturbationType] + "_results_" + netw + ".log")));
            resultsFile.println("# Epsilon: " + epsilon);
        }
        catch (Exception e)
        {
            System.out.println(typeDescs[perturbationType]
                    + " failed to initialise.");
        }
    }

    /**
     * Close file handles
     */
    public void finish()
    {
        thetaFile.close();
        traceFile.close();
        resultsFile.close();
    }

    /**
     * Update thresholds parameter theta in descent direction using 2RDSA based
     * update iteration that uses 3 function evaluations
     * This function is for the case of Uniform perturbations based 2RDSA
     * 
     * @param simNr
     *            2RDSA based outer loop iteration index
     */    
    void updateParamsUnif(int simNr)
    {
         // Generate the Hessian estimate 
        for (int i = 0; i < HessianEstimate.length; i++)
        {
            double perturbationFactor = 2.5 * (Delta[i] * Delta[i] - Math.pow(
                    3, -1));
            HessianEstimate[i][i] = HessianEstimate[i][i]
                    + cn
                    * (4.5 * perturbationFactor * (Zplus + Zminus - 2 * Z)
                            / (delta * delta) - HessianEstimate[i][i]);

            if (HessianEstimate[i][i] <= eta)
            {
                HessianEstimate[i][i] = eta;
            }
        }

        // Invert the Hessian matrix
        Matrix Pij = new Matrix(HessianEstimate);
        Pij = Pij.inverse();
        // Project the Hessian inverse onto the set of positive defnite matrices
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

        // Update theta using a Newton step
        for (int i = 0; i < theta.length; i++)
        {
            double prod = 0;
            for (int j = 0; j < HessianEstimate.length; j++)
            {
                prod += Pij.get(i, j) * 3.0 * (Zminus - Zplus) / (2 * delta)
                        * Delta[j];
            }
            theta[i] -= an * prod;

            // Projection to ensure stability of theta recursion
            if (theta[i] <= THETA_MIN[i])
                theta[i] = THETA_MIN[i];
            if (theta[i] >= THETA_MAX[i])
                theta[i] = THETA_MAX[i];
        }

        // Update stepsizes and perturbation constants
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
     * This function is for the case of asymmetric Bernoulli perturbations based 2RDSA
     * The overall structure is very similar to updateParamsUnif
     * 
     * @param simNr
     *            2RDSA based outer loop iteration index
     */    
    void updateParamsAsymBer(int simNr)
    {
        double tau = (1 + epsilon) * (1 + Math.pow(1 + epsilon, 3))
                / (2 + epsilon);
        double kappa = tau * (1 - Math.pow(1 + epsilon, 2) / tau);
        for (int i = 0; i < HessianEstimate.length; i++)
        {
            double perturbationFactor = Math.pow(kappa, -1)
                    * (Delta[i] * Delta[i] - (1 + epsilon));
            HessianEstimate[i][i] = HessianEstimate[i][i]
                    + cn
                    * ((perturbationFactor * (Zplus + Zminus - 2 * Z) / (delta * delta)) - HessianEstimate[i][i]);

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
                prod += Pij.get(i, j) * (1 / (1 + epsilon)) * (Zminus - Zplus)
                        / (2 * delta) * Delta[j];
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
        traceFile.println("############# Training "
                + typeDescs[perturbationType] + " ################");
        System.out.println("############# Training "
                + typeDescs[perturbationType] + " ################");

        for (int simNr = 1; simNr < trainingBudget + 1; simNr++)
        {
            traceFile.println("############################ " + simNr
                    + " ################################");
            System.out.println("############################ " + simNr
                    + " ################################");

            int thetaTypeforSimulation = simNr % 3;
            if (thetaTypeforSimulation == 0)
            {
                setPerturbedThetas();
            }
            double thetaSim[];
            thetaSim = getTheta(thetaTypeforSimulation);
            double avgCostSample = gldRunnerObject.run(thetaSim,
                    trajectoryLengthTrainingPhase, false, "");

            gldRunnerObject.resetSimulator();

            if (simNr % 3 == 0)
            {
                Zminus = avgCostSample;
                traceFile.println("Zminus:" + Zminus);
                if (perturbationType == UNIF)
                {
                    updateParamsUnif(simNr);
                }
                else
                {
                    updateParamsAsymBer(simNr);
                }
            }
            else if (simNr % 3 == 2)
            {
                Zplus = avgCostSample;
                traceFile.println("Zplus:" + Zplus);
            }
            else if (simNr % 3 == 1)
            {
                Z = avgCostSample;
                traceFile.println("Z:" + Z);
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
        traceFile.println("############# Testing "
                + typeDescs[perturbationType] + " ################");
        System.out.println("############# Testing "
                + typeDescs[perturbationType] + " ################");

        for (int simNr = 1; simNr < numReplicationsForTesting + 1; simNr++)
        {
            traceFile.println("############################ " + simNr
                    + " ################################");
            System.out.println("############################ " + simNr
                    + " ################################");

            String filePrefix = typeDescs[perturbationType] + "_"
                    + Integer.toString(simNr);
            double avgCostSample = gldRunnerObject.run(theta,
                    trajectoryLengthTestingPhase, true, filePrefix);

            resultsFile.println(avgCostSample);
            System.out.println("Avg Value: " + avgCostSample);

            gldRunnerObject.resetSimulator();
        }
    }
}
