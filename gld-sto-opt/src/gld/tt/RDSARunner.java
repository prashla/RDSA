package gld.tt;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

import gld.*;

/**
 * Wrapper for running 2SPSA as well as 2RDSA (both AsymmetricBernoulli as well
 * as Uniform perturbation variants)
 * 
 * @see TwoRDSA_OuterLoop, TwoSPSA_OuterLoop
 * @version 1.0 31 Jul 2016
 * @author Prashanth L.A.
 */
public class RDSARunner
{
    protected static final int UNIF = 0, ASYMBER = 1;

    public static void main(String[] params)
    {
        /**
         * the number of function evaluations in the training phase
         */
        int trainingBudget = 2000;
        /**
         * After the training phase, the thresholds are fixed and then a number
         * of independent simulations (that this variable is set to) are run and
         * the empirical average cost from each simulation is recorded
         */
        int numReplicationsForTesting = 50;
        /**
         * length of each simulated trajectory during training phase
         */
        int trajectoryLengthTrainingPhase = 250;
        /**
         * length of each simulated trajectory during testing phase
         */
        int trajectoryLengthTestingPhase = 10000;

        String mapdir = "/home/prashanth/Documents/github/rdsa-office/gld-sto-opt/maps/";
        String map = "corridor10";

        /*
         * Run 2SPSA with PTLC as the underlying TLC algorithm
         */
        // Create an object that run the outer loop of 2SPSA+PTLC
        TwoSPSA_OuterLoop spsaObject = new TwoSPSA_OuterLoop(params, mapdir,
                map, trainingBudget, numReplicationsForTesting,
                trajectoryLengthTrainingPhase, trajectoryLengthTestingPhase);

        // open file handles for trace, theta and results log - these are stored
        // in the gld folder
        spsaObject.initFiles(map);

        // This constitutes the training phase where 2SPSA is run to optimize
        // the thresholds
        spsaObject.optimizeThresholds();

        // Test the parameter obtained at the end of training phase and log the
        // avg Cost samples for each replication
        spsaObject.testThresholds();

        // close all file handles
        spsaObject.finish();

        /*
         * Wrapper for running 2RDSA with Asymmetric Bernoulli perturbations For
         * running the Uniform perturbations variant, change the last parameter
         * in the constructor below to UNIF
         */
        TwoRDSA_OuterLoop rdsaObject = new TwoRDSA_OuterLoop(params, mapdir,
                map, trainingBudget, numReplicationsForTesting,
                trajectoryLengthTrainingPhase, trajectoryLengthTestingPhase,
                ASYMBER);

        // open file handles for trace, theta and results log
        rdsaObject.initFiles(map);

        // This constitutes the training phase where 2RDSA is run to optimize
        // the thresholds
        rdsaObject.optimizeThresholds();

        // Test the parameter obtained at the end of training phase and log the
        // avg Cost samples for each replication
        rdsaObject.testThresholds();

        // close all file handles
        rdsaObject.finish();
    }
}
