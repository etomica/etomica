/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import java.util.List;
import java.util.Set;

import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;
import etomica.chem.elements.ElementSimple;
import etomica.graph.model.Graph;
import etomica.graph.property.FFTDecomposition;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.ConfigurationClusterMove;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

/**
 * MSMC computations (via overlap sampling) of the bridge diagram(s) within the
 * f-bond-only formulation.  Diagrams are filtered from the full set of
 * biconnected diagrams.  It is also possible to select diagrams with
 * un-FFT-able biconnected components of a certain size.
 * 
 * Kate Shaul
 * Andrew Schultz
 */
public class VirialLJBridge {


    public static void main(String[] args) {

        VirialLJParam params = new VirialLJParam();
        ParseArgs.doParseArgs(params, args);
        boolean isCommandLine = args.length > 0;
        double temperature = params.temperature;
        final int nPoints = params.numMolecules;
        final double sigmaHSRef = params.sigmaHSRef;
        long steps = params.numSteps;
        final int iGraph = params.iGraph;

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        final double refIntegral = Standard.BHS(nPoints, sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+refIntegral);
        System.out.println("Lennard Jones overlap sampling f-bond decomposition of B"+nPoints+" at T= "+temperature);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        Potential2Spherical pTarget = new P2LennardJones(space,1.0,1.0);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        
        if (iGraph == -1) {
            System.out.println("all bridge diagrams");
        }
        else if (iGraph < 0) {
            System.out.println("all bridge diagrams with segment size >= "+(-iGraph));
        }
        else {
            System.out.println("bridge diagram "+iGraph);
        }
        VirialDiagrams diagrams = new VirialDiagrams(nPoints, false, false);
        diagrams.setDoReeHoover(false);
        diagrams.setDoShortcut(true);
        Set<Graph> graphs = diagrams.getMSMCGraphs(true, false);
        FFTDecomposition isFFT = new FFTDecomposition();
        List<Integer> segments = isFFT.getSegments();
        Set<Graph> bridgeGraphs = VirialDiagrams.makeGraphList();
        // 4th order: 31
        // 5th order: around 500
        // 6th order: around 16,000
        // 7th order: around 1,000,000
        int[] segmentDividers = new int[]{0,0,0,30,32,1000,2000000};
        int gCount = 0;
        for (Graph g : graphs) {
            isFFT.check(g);
            if (segments.size() == 0) continue;
            long maxSegment = segments.get(0);
            segments.clear();
            int segmentSize = -1;
            for (int i=3; i<nPoints+1; i++) {
                if (segmentDividers[i-1] < maxSegment && segmentDividers[i] > maxSegment) {
                    segmentSize = i;
                }
            }
//            System.out.println(maxSegment+" "+segmentSize);
            if (iGraph < 0 && segmentSize == -iGraph) {
                bridgeGraphs.add(g);
                System.out.println(g);
            }
            else if (gCount == iGraph) {
                bridgeGraphs.add(g);
                System.out.println(g);
            }
            gCount++;
        }

        ClusterSum targetCluster = diagrams.makeVirialCluster(bridgeGraphs, fTarget, null);
        System.out.println();

        int[][] refBondList = new int[nPoints-1][2];
        for (int i=0; i<nPoints-1; i++) {
            int first = i;
            int second = i+1;
            if (first == 1) first++;
            if (second == 2) first++;
            if (second == nPoints) second = 1;
            refBondList[i][0] = first;
            refBondList[i][1] = second;
        }
        
        
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println(steps+" steps (1000 blocks of "+(steps/1000)+")");
		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("LJ")), temperature,refCluster,targetCluster);
        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        
        ConfigurationClusterMove config = new ConfigurationClusterMove(space, sim.getRandom());
        config.initializeCoordinates(sim.box[1]);
        
        steps /= 1000;
        long t1 = System.currentTimeMillis();
        // if running interactively, don't use the file
        String refFileName = isCommandLine ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
//        sim.setRefPref(1.0082398078547523);
        sim.initRefPref(refFileName, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/40);
        
        System.out.println("equilibration finished");


        
        if (!isCommandLine) {
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.getController().getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    System.out.println("abs average: "+ratio*refIntegral+" error: "+error*refIntegral);
                    if (ratio == 0 || Double.isNaN(ratio)) {
                        throw new RuntimeException("oops");
                    }
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.integratorOS.setNumSubSteps((int)steps);
        sim.setAccumulatorBlockSize(steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorOS), 1000);
        long t2 = System.currentTimeMillis();

        System.out.println("ideal reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        sim.printResults(refIntegral);
        System.out.println("time: "+(t2-t1)/1000.0);
	}

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
    	
    	// number of molecules in simulation (e.g., 2 for B2 calculation)

        public int numMolecules = 6;
        public double temperature = 1.0;
        public long numSteps = 1000000000;
        public double sigmaHSRef = 1.5;
        public int iGraph = -6;
    }
    
}
