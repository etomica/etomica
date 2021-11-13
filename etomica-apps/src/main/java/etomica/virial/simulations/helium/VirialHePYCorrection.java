/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.helium;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P2HeSimplified;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.*;
import etomica.virial.integralequation.PYGenerator;
import etomica.virial.simulations.SimulationVirialOverlap2;

import java.util.Set;

/**
 * Computes corrections to Percus-Yevick approximations of B4 and B5 for a helium pair potential.
 * 
 * Select the compressibility or virial route via boolean variable compressibility.
 * 
 * Use semiclassical boolean to change pair potential to the quadratic Feynman-Hibbs potential.
 * 
 * Use calcApprox boolean to calculate quantities with the approximate potential.
 * Use subtractApprox boolean to compute the correction to that approximate potential.
 * 
 * If determining which option is most efficient via short calculations to estimate standard error, 
 * maintain a 50-50 split of steps between reference and target during data collection.
 * 
 * @author kate
 * @author Andrew Schultz
 */
public class VirialHePYCorrection {


    public static void main(String[] args) {

        VirialHePYCParam params = new VirialHePYCParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize parameters here
            params.nPoints = 6;
            params.semiClassical = true;
            params.calcApprox = false;
            params.subtractApprox = true;
            params.temperature = 1000;
            params.numSteps = 1000000;
        }
        
    	final int nPoints = params.nPoints;
        final double temperatureK = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        if (sigmaHSRef < 0) {
            sigmaHSRef = 2.4 + 120/(100+temperatureK);
        }
        final boolean semiClassical = params.semiClassical;
        final double refFrac = params.refFrac;
        final boolean subtractApprox = params.subtractApprox;
        final boolean calcApprox = !subtractApprox && params.calcApprox;

        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);

        System.out.println("Overlap sampling for He pair potential of Przybytek et al. (2010) at " + temperatureK + " K");
        if (semiClassical) {
        	System.out.println("Quadratic Feymann-Hibbs version of potential employed.");
        }
        
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("  B"+nPoints+"HS: "+HSB[nPoints]);
        if (calcApprox) System.out.println("Calculating coefficients for approximate potential");
        if (subtractApprox) {
            System.out.println("computing difference from approximate He");
        }
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        MayerGeneralSpherical fTarget;
        MayerGeneralSpherical fTargetApprox;
        if (semiClassical) {
            P2HeSimplified p2cApprox = new P2HeSimplified(space);
            Potential2Spherical p2Approx = p2cApprox.makeQFH(temperature);
            
            P2HePCKLJS p2c = new P2HePCKLJS(space);
            Potential2Spherical p2 = p2c.makeQFH(temperature);

            fTarget = new MayerGeneralSpherical(calcApprox ? p2Approx : p2);
            fTargetApprox = new MayerGeneralSpherical(p2Approx);
        } else {
            P2HeSimplified p2Approx = new P2HeSimplified(space);
            
            P2HePCKLJS p2 = new P2HePCKLJS(space);

            fTarget = new MayerGeneralSpherical(calcApprox ? p2Approx : p2);
            fTargetApprox = new MayerGeneralSpherical(p2Approx);
        }

        IsomorphismFilter.DEBUG_MODE = false;
        Set<Graph> correction = VirialDiagrams.makeGraphList();
        System.out.println("starting");
        correction.addAll(PYGenerator.getPYCorrection((byte)nPoints));
        System.out.println("here");
        MetadataImpl.rootPointsSpecial = false;
        IsoFree isoFree = new IsoFree();
        MulScalar mulScalar = new MulScalar();
        MulScalarParameters msp = new MulScalarParameters(-1, nPoints);
        Set<Graph> correctionB = VirialDiagrams.makeGraphList();
        correctionB.addAll(isoFree.apply(mulScalar.apply(correction, msp), null));

        VirialDiagrams diagrams = new VirialDiagrams(nPoints, false, false);
        diagrams.setDoReeHoover(true);
        diagrams.setAllPermutations(false);
        diagrams.setDoShortcut(true);
        ClusterSum fullTargetCluster = diagrams.makeVirialCluster(correctionB, fTarget, null);

        ClusterAbstract targetCluster = null;
        if (subtractApprox) {
            final ClusterSum[] targetSubtract = new ClusterSum[1];
            ClusterBonds[] minusBonds = fullTargetCluster.getClusters();
            double[] wMinus = fullTargetCluster.getWeights();
            targetSubtract[0] = new ClusterSum(minusBonds, wMinus, new MayerFunction[]{fTargetApprox});
            targetCluster = new ClusterDifference(fullTargetCluster, targetSubtract);
        }
        else {
            targetCluster = fullTargetCluster;
        }

        
        targetCluster.setTemperature(temperature);

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoReeHoover(true);
        rigidDiagrams.setDoShortcut(true);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);

        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");

        ISpecies species = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("He")));
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints}, temperature, refCluster, targetCluster);
        sim.init();
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        
        long t1 = System.currentTimeMillis();
        // The diagram which constitutes B4-B4PY is zero for an overlapped configuration.  Without ConfigurationClusterMove, the initial configuration will be overlapped, and gamma/pi will be zero.
        // Such configurations will not be visited later, precisely because pi is zero.
        double r = 4;
        for (int i=1; i<nPoints; i++) {
            Vector v = sim.box[1].getLeafList().get(i).getPosition();
            v.setX(0, r*Math.cos(2*(i-1)*Math.PI/(nPoints-1)));
            v.setX(1, r*Math.sin(2*(i-1)*Math.PI/(nPoints-1)));
        }
        sim.box[1].trialNotify();
        sim.box[1].acceptNotify();
        if (sim.box[1].getSampleCluster().value(sim.box[1]) == 0) {
            throw new RuntimeException("oops");
        }
        
        sim.integratorOS.setNumSubSteps(1000);
        
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        steps /= 1000;
        sim.setAccumulatorBlockSize(steps);
        
        System.out.println();
        String refFileName = null;
        if (isCommandline) {
            // if running interactively, don't use the file
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints+"_2b_"+tempString;
            refFileName += semiClassical ? "_sc" : "_c";
            if (calcApprox) {
                refFileName += "a";
            }
            else if (subtractApprox) {
                refFileName += "sa";
            }
            refFileName += "PY";
            refFileName += "C";
        }

        sim.initRefPref(refFileName, steps/40);
        sim.equilibrate(refFileName, steps/20);
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
System.out.println("equilibration finished");

        sim.integratorOS.setNumSubSteps((int)steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
sim.getController().runActivityBlocking(ai);

        System.out.println();
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        System.out.println();
        
        sim.printResults(HSB[nPoints]);

        DataGroup allYourBase = (DataGroup)sim.accumulators[1].getData();
        IData averageData = allYourBase.getData(sim.accumulators[1].AVERAGE.index);
        IData errorData = allYourBase.getData(sim.accumulators[1].ERROR.index);
        IData covarianceData = allYourBase.getData(sim.accumulators[1].BLOCK_COVARIANCE.index);
        int n = 0;
        double correlationCoef = covarianceData.getValue(n+1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue((n+2)*(n+2)-1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
        System.out.print(String.format("diagram "+nPoints+"bc average: %20.15e error: %9.4e ocor: %6.4f\n",
                averageData.getValue(0), errorData.getValue(0), correlationCoef));

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHePYCParam extends ParameterBase {
        // don't change these
        public int nPoints = 4;
        public double temperature = 100;
        public long numSteps = 10000000;
        public double refFrac = -1;
        public double sigmaHSRef = -1;
        public boolean semiClassical = false;
        public boolean calcApprox = false;
        public boolean subtractApprox = false;
    }
}
