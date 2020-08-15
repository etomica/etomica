/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2SoftSphere;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.util.Arrays;

/**
 * LJ simulation using Mayer sampling to integrals that show up in series
 * expansion.
 */
public class VirialLJSeries {


    public static void main(String[] args) {
        VirialLJParam params = new VirialLJParam();
        if (args.length > 0) {
            params.writeRefPref = true;
            ReadParameters readParameters = new ReadParameters(args[0], params);
            readParameters.readParameters();
        }
        
        runVirial(params);
    }
    
    public static void runVirial(VirialLJParam params) {
        final int nPoints = params.nPoints;
        double temperature = 1;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        int[] bondList = params.bondList;

        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        System.out.println("Lennard Jones series overlap sampling B"+nPoints+" at T="+temperature);
        System.out.println("bond list: "+Arrays.toString(bondList));
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        Potential2Spherical pTarget = new P2SoftSphere(space,1.0,4.0,12);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        boolean[] used = new boolean[10];
        int[] bondMap = new int[10];
        for (int i=1; i<10; i++) {
            bondMap[i] = -1;
        }
        int iBond = 0;
        for (int i=0; i<bondList.length; i++) {
            if (!used[bondList[i]]) {
                used[bondList[i]] = true;
                iBond++;
            }
        }
        MayerFunction[] allF = new MayerFunction[iBond];
        iBond = 0;
        if (used[0]) {
            allF[0] = fTarget;
            bondMap[0] = 0;
            iBond++;
        }
        for (int i=1; i<10; i++) {
            if (used[i]) {
                bondMap[i] = iBond;
                allF[iBond] = new MayerSSSeries(6*i);
                iBond++;
            }
        }
        int[] newBondList = new int[bondList.length];
        for (int i=0; i<bondList.length; i++) {
            newBondList[i] = bondMap[bondList[i]];
        }
        ClusterAbstract targetCluster = Standard.virialSeriesCluster(nPoints, allF, newBondList);
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("A")), temperature,refCluster,targetCluster);
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20);
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            return;
        }

        
        sim.integratorOS.setNumSubSteps(1000);
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
//        sim.setRefPref(1.0082398078547523);
        sim.initRefPref(refFileName, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/40);
        
        System.out.println("equilibration finished");

//        IAction progressReport = new IAction() {
//            public void actionPerformed() {
//                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
//                double ratio = sim.dsvo.getDataAsScalar();
//                double error = sim.dsvo.getError();
//                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
//            }
//        };
//        sim.integratorOS.addIntervalAction(progressReport);
//        sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));
        
//        VirialHistogram virialHistogram = new VirialHistogram(new P2HardSphere(space, sigmaHSRef, false), pTarget, sim.box[1]);
//        virialHistogram.setBinFac(100);
//        sim.integrators[1].addIntervalAction(virialHistogram);

        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorOS), steps);
        
//        long[][] histogram = virialHistogram.getHistogram();
//        int numNegBins = virialHistogram.getNumNegBins();
//        double zeroOffset = virialHistogram.getZeroOffset();
//        double binFac = virialHistogram.getBinFac();
//        if (false) {
//        try {
//            FileWriter writer = new FileWriter("nover.dat");
//            for (int i=0; i<numNegBins; i++) {
//                if (histogram[0][i] == 0) continue;
//                writer.write(-(Math.exp(i/binFac+Math.log(zeroOffset))-zeroOffset)+" "+histogram[0][i]+"\n");
//            }
//            for (int i=numNegBins; i<histogram[0].length; i++) {
//                if (histogram[0][i] == 0) continue;
//                writer.write((Math.exp((i-numNegBins)/binFac+Math.log(zeroOffset))-zeroOffset)+" "+histogram[0][i]+"\n");
//            }
//            writer.close();
//            writer = new FileWriter("over.dat");
//            for (int i=0; i<numNegBins; i++) {
//                if (histogram[1][i] == 0) continue;
//                writer.write(-(Math.exp(i/binFac+Math.log(zeroOffset))-zeroOffset)+" "+histogram[1][i]+"\n");
//            }
//            for (int i=numNegBins; i<histogram[1].length; i++) {
//                if (histogram[1][i] == 0) continue;
//                writer.write((Math.exp((i-numNegBins)/binFac+Math.log(zeroOffset))-zeroOffset)+" "+histogram[1][i]+"\n");
//            }
//            writer.close();
//        }
//        catch (IOException e) {
//            throw new RuntimeException(e);
//        }
//        }

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());

        sim.printResults(HSB[nPoints]);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
        public int nPoints = 3;
        public int[] bondList = new int[]{0,9,9};
        public long numSteps = 10000;
        public double sigmaHSRef = 1.5;
        public boolean writeRefPref = false;
    }
}
