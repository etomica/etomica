/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.potential.P2EffectiveFeynmanHibbs;
import etomica.potential.P2HydrogenPatkowskiIso;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

import java.util.Set;

/**
 * Computes corrections to Percus-Yevick approximations of B4 and B5 for a hydrogen pair potential.
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
public class VirialH2PYCorrection {


    public static void main(String[] args) {

        VirialH2PYCParam params = new VirialH2PYCParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize parameters here
            params.nPoints = 2;            
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
        final double refFrac = params.refFrac;
        
        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);

//        System.out.println("Overlap sampling for He pair potential of Przybytek et al. (2010) at " + temperatureK + " K");
        
        System.out.println("Quadratic Feymann-Hibbs version of potential employed.");        
        
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("  B"+nPoints+"HS: "+HSB[nPoints]);        
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        MayerGeneralSpherical fTarget;
        
        P2EffectiveFeynmanHibbs p2SemiClassical = null;
        final P2HydrogenPatkowskiIso p2patIso = new P2HydrogenPatkowskiIso(space);        

        p2SemiClassical = new P2EffectiveFeynmanHibbs(space, p2patIso);
        p2SemiClassical.setMass(Hydrogen.INSTANCE.getMass()*2);
        p2SemiClassical.setTemperature(temperature);
        Potential2Spherical p2 = p2SemiClassical;
        fTarget = new MayerGeneralSpherical(p2);

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
        targetCluster = fullTargetCluster;        
        targetCluster.setTemperature(temperature);

        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        rigidDiagrams.setDoReeHoover(true);
        rigidDiagrams.setDoShortcut(true);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);

        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("He")), temperature, refCluster, targetCluster);
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
            refFileName += "_scPYC";
        }

        sim.initRefPref(refFileName, steps/40);
        sim.equilibrate(refFileName, steps/20);
        
        System.out.println("equilibration finished");
        
        sim.integratorOS.setNumSubSteps((int)steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
        }
sim.getController2().runActivityBlocking(new etomica.action.activity.ActivityIntegrate2(sim.integratorOS), 1000);

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
    public static class VirialH2PYCParam extends ParameterBase {
        // don't change these
        public int nPoints = 4;
        public double temperature = 100;
        public long numSteps = 10000000;
        public double refFrac = -1;
        public double sigmaHSRef = -1;        
    }
}
