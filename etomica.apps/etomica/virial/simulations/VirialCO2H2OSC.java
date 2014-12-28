/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.api.IAtomList;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.models.co2.P2CO2H2OWheatley;
import etomica.models.co2.P2CO2Hellmann;
import etomica.models.water.P2WaterSzalewicz;
import etomica.models.water.P2WaterSzalewicz.Component;
import etomica.potential.P2SemiclassicalAtomic;
import etomica.potential.P2SemiclassicalAtomic.AtomInfo;
import etomica.potential.PotentialMolecularMonatomic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Kelvin;
import etomica.util.DoubleRange;
import etomica.util.HistogramNotSoSimple;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterCoupledAtomFlipped;
import etomica.virial.ClusterWheatleyHS;
import etomica.virial.ClusterWheatleyMultibodyMix;
import etomica.virial.ClusterWheatleySoftMix;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;
import etomica.virial.MayerFunctionMolecularThreeBody;
import etomica.virial.MayerFunctionNonAdditive;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.PotentialNonAdditive;
import etomica.virial.cluster.Standard;

/**
 * Computes CO2-H2O mixture virial coefficients using ab-initio potentials
 * for both components.  Classical and semiclassical coefficients can be
 * computed.
 * 
 * For water, only pairwise-additive contributions are computed.
 * 
 * @author Andrew Schultz
 */
public class VirialCO2H2OSC {


    public static void main(String[] args) {

        VirialCO2SCParam params = new VirialCO2SCParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize parameters here
            params.nTypes = new int[]{1,1};
            params.level = Level.SEMICLASSICAL_TI;
            params.temperature = 300;
            params.numSteps = 1000000;
            params.sigmaHSRef = 5;
//            params.nonAdditive = true;
        }

    	final int[] nTypes = params.nTypes;
        final int nPoints = params.nTypes[0] + params.nTypes[1];
        final double temperatureK = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        if (sigmaHSRef < 0 && false) {
            sigmaHSRef = 2.4 + 120/(100+temperatureK);
        }
        final Level level = params.level;
        final double refFrac = params.refFrac;
        final boolean nonAdditive = params.nonAdditive;

        final double HSB = Standard.BHS(nPoints, sigmaHSRef);

        System.out.println("Overlap sampling for CO2 pair potential of Hellmann (2014) at " + temperatureK + " K");
        if (level == Level.SEMICLASSICAL_FH) {
        	System.out.println("Quadratic Feymann-Hibbs effective potential employed.");
        }
        if (level == Level.SEMICLASSICAL_TI) {
            System.out.println("Takahashi-Imada effecitve potential employed.");
        }

        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("  B"+nPoints+"HS: "+HSB);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        SpeciesSpheresRotating speciesCO2 = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+Oxygen.INSTANCE.getMass()*2));
        SpeciesSpheresRotating speciesH2O = new SpeciesSpheresRotating(space, new ElementSimple("H2O", Oxygen.INSTANCE.getMass()+Hydrogen.INSTANCE.getMass()*2));
        speciesH2O.setAxisSymmetric(false);
        P2CO2Hellmann p2cCO2 = new P2CO2Hellmann(space, P2CO2Hellmann.Parameters.B);
        IPotentialAtomic p2aCO2 = level == Level.CLASSICAL ? null : (level == Level.SEMICLASSICAL_FH ? p2cCO2.makeSemiclassical(temperature) : new P2SemiclassicalAtomic(space, p2cCO2, temperature));

        P2CO2H2OWheatley p2cCO2H2O = new P2CO2H2OWheatley(space);
        IPotentialAtomic p2aCO2H2O = level == Level.CLASSICAL ? null : (level == Level.SEMICLASSICAL_FH ? p2cCO2H2O.makeSemiclassical(temperature) : p2cCO2H2O.makeSemiclassicalTI(temperature));
        
        P2WaterSzalewicz p2cH2O = new P2WaterSzalewicz(space, 2);
        IPotentialAtomic p2aH2O = level == Level.CLASSICAL ? null : (level == Level.SEMICLASSICAL_FH ? p2cH2O.makeSemiclassical(temperature) : new P2SemiclassicalAtomic(space, p2cH2O, temperature));

        PotentialMolecularMonatomic p2CO2 = new PotentialMolecularMonatomic(space, level==Level.CLASSICAL ? p2cCO2 : p2aCO2);
        PotentialMolecularMonatomic p2H2O = new PotentialMolecularMonatomic(space, level==Level.CLASSICAL ? p2cH2O : p2aH2O);
        PotentialMolecularMonatomic p2CO2H2O = new PotentialMolecularMonatomic(space, level==Level.CLASSICAL ? p2cCO2H2O : p2aCO2H2O);
        MayerGeneral fCO2 = new MayerGeneral(p2CO2);
        MayerGeneral fCO2H2O = new MayerGeneral(p2CO2H2O);
        MayerGeneral fH2O = new MayerGeneral(p2H2O);
        MayerFunction[][] allF = new MayerFunction[][]{{fCO2,fCO2H2O},{fCO2H2O,fH2O}};

        ClusterAbstract targetCluster = new ClusterWheatleySoftMix(nPoints, nTypes, allF, 1e-12);
        if (nTypes[1] == 2 && nPoints == 2) {
            // pure B2 for water.  we need flipping.
            // additive B3 for water should be fine and biconnectivity will help with mixture coefficients.
            ((ClusterWheatleySoftMix)targetCluster).setDoCaching(false);
            targetCluster = new ClusterCoupledAtomFlipped(targetCluster, space, 20);
        }
        if (nonAdditive) {

            P2WaterSzalewicz p23cH2O = new P2WaterSzalewicz(space, 2);
            p23cH2O.setComponent(Component.NON_PAIR);
            IPotentialAtomic p23aH2O = level == Level.CLASSICAL ? null : (level == Level.SEMICLASSICAL_FH ? p23cH2O.makeSemiclassical(temperature) : new P2SemiclassicalAtomic(space, p23cH2O, temperature));
            PotentialMolecularMonatomic p23H2O = new PotentialMolecularMonatomic(space, level==Level.CLASSICAL ? p23cH2O : p23aH2O);

            P2WaterSzalewicz p3cH2O = new P2WaterSzalewicz(space, 3);
            p3cH2O.setComponent(Component.NON_PAIR);
            IPotentialAtomic p3aH2O = level == Level.CLASSICAL ? null : (level == Level.SEMICLASSICAL_FH ? p3cH2O.makeSemiclassical(temperature) : null);
            PotentialMolecularMonatomic p3H2O = new PotentialMolecularMonatomic(space, level==Level.CLASSICAL ? p3cH2O : p3aH2O);
            MayerFunctionMolecularThreeBody f3H2O = new MayerFunctionMolecularThreeBody(new PotentialNonAdditive(new IPotentialMolecular[]{p23H2O,p3H2O}));
            MayerFunctionNonAdditive[][][] allFNA = new MayerFunctionNonAdditive[2][2][2];
            allFNA[1][1][1] = f3H2O;
            targetCluster = new ClusterWheatleyMultibodyMix(nPoints, nTypes, allF, allFNA, 1e-12, true);
            ((ClusterWheatleyMultibodyMix)targetCluster).setDoCaching(false);
            targetCluster = new ClusterCoupledAtomFlipped(targetCluster, space, 20);
        }

        
        targetCluster.setTemperature(temperature);

        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);

        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
 		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{speciesCO2,speciesH2O}, nTypes, temperature, refCluster, targetCluster);
        sim.init();
        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        if (level == Level.SEMICLASSICAL_TI && false) {
            if (true) throw new RuntimeException("implement me for anything other than CO2-CO2");
            final IVectorMutable[] rv = new IVectorMutable[4];
            for (int i=0; i<4; i++) {
                rv[i] = space.makeVector();
            }
            double om = Oxygen.INSTANCE.getMass();
            double bondLength = p2cCO2.posB[5] - p2cCO2.posB[1];
            rv[0].setX(0, om*bondLength*bondLength*0.25);
            rv[0].setX(1, om*bondLength*bondLength*0.25);
            ((P2SemiclassicalAtomic)p2aCO2).setAtomInfo(speciesCO2.getLeafType(), new AtomInfo() {
                public IVector[] getMomentAndAxes(IAtomOriented molecule) {
                    // rv[0,2] = 0
                    // rv[3] is the orientation
                    rv[3].E(molecule.getOrientation().getDirection());
                    // rv[1] is an axis perpendicular to rv[3]
                    rv[1].E(0);
                    if (Math.abs(rv[3].getX(0)) < 0.5) {
                        rv[1].setX(0, 1);
                    }
                    else if (Math.abs(rv[3].getX(1)) < 0.5) {
                        rv[1].setX(1, 1);
                    }
                    else {
                        rv[1].setX(2, 1);
                    }
                    rv[2].Ea1Tv1(rv[1].dot(rv[3]), rv[3]);
                    rv[1].ME(rv[2]);
                    // rv[2] is an axis perpendicular to rv[3] and rv[1]
                    rv[2].E(rv[1]);
                    rv[2].XE(rv[3]);
                    return rv;
                }
            });
        }
        
        if (nonAdditive) {
            double r = 8;
            IAtomList atoms = sim.box[1].getLeafList();
            for (int i=1; i<nPoints; i++) {
                IVectorMutable pos = atoms.getAtom(i).getPosition();
                double theta = 2*(i-1)*Math.PI/(nPoints-1);
                pos.setX(0, r*Math.cos(theta));
                pos.setX(1, r*Math.sin(theta));
            }
        }

        
        long t1 = System.currentTimeMillis();
        
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
            refFileName += level == Level.CLASSICAL ? "_c" : (level == Level.SEMICLASSICAL_FH ? "_fh" : "_ti");
            refFileName += "C";
        }

        if (nonAdditive) {
            sim.integrators[1].getEventManager().addListener(new IIntegratorListener() {
                
                public void integratorStepStarted(IIntegratorEvent e) {}
                
                public void integratorStepFinished(IIntegratorEvent e) {
                    if (Math.random() > 0.0001) return;
                    CoordinatePairSet cpairs = sim.box[1].getCPairSet();
                    System.out.println(String.format("r2: %8.3f  %8.3f  %8.3f", Math.sqrt(cpairs.getr2(0,1)), Math.sqrt(cpairs.getr2(0,2)), Math.sqrt(cpairs.getr2(1,2))));
                    DataGroup allYourBase = (DataGroup)sim.accumulators[1].getData();
                    IData averageData = allYourBase.getData(sim.accumulators[1].AVERAGE.index);
                    IData errorData = allYourBase.getData(sim.accumulators[1].ERROR.index);
                    IData correlationData = allYourBase.getData(sim.accumulators[1].BLOCK_CORRELATION.index);
                    long nBlocks = sim.accumulators[1].getBlockCount();
                    System.out.print(String.format("target average(%d): %20.15e error: %9.4e cor: %6.4f\n", nBlocks,
                            averageData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
                }
                
                public void integratorInitialized(IIntegratorEvent e) {}
            });
        }

        sim.initRefPref(refFileName, steps/40);
        sim.equilibrate(refFileName, steps/20);
        
        System.out.println("equilibration finished");
        
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }

        final HistogramNotSoSimple targHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple targPiHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        int nBins = 100;
        double dx = sigmaHSRef/nBins;
        final HistogramNotSoSimple hist = new HistogramNotSoSimple(nBins, new DoubleRange(dx*0.5, sigmaHSRef+dx*0.5));
        final HistogramNotSoSimple piHist = new HistogramNotSoSimple(nBins, new DoubleRange(dx*0.5, sigmaHSRef+dx*0.5));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();
        IIntegratorListener histListenerRef = new IIntegratorListener() {
            public void integratorStepStarted(IIntegratorEvent e) {}
            
            public void integratorStepFinished(IIntegratorEvent e) {
                double r2Max = 0;
                CoordinatePairSet cPairs = sim.box[0].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }
                double v = finalTargetCluster.value(sim.box[0]);
                hist.addValue(Math.sqrt(r2Max), v);
                piHist.addValue(Math.sqrt(r2Max), Math.abs(v));
            }
            
            public void integratorInitialized(IIntegratorEvent e) {
            }
        };
        IIntegratorListener histListenerTarget = new IIntegratorListener() {
            public void integratorStepStarted(IIntegratorEvent e) {}
            
            public void integratorStepFinished(IIntegratorEvent e) {
                double r2Max = 0;
                double r2Min = Double.POSITIVE_INFINITY;
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij < r2Min) r2Min = r2ij;
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }

                double v = finalTargetCluster.value(sim.box[1]);
                double r = Math.sqrt(r2Max);
                if (r > 1) {
                    r = Math.log(r);
                }
                else {
                    r -= 1;
                }
                targHist.addValue(r, v);
                targPiHist.addValue(r, Math.abs(v));
            }

            public void integratorInitialized(IIntegratorEvent e) {}
        };

        if (params.doHist) {
            IIntegratorListener histReport = new IIntegratorListener() {
                public void integratorInitialized(IIntegratorEvent e) {}
                public void integratorStepStarted(IIntegratorEvent e) {}
                public void integratorStepFinished(IIntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.println("**** reference ****");
                    double[] xValues = hist.xValues();
                    double[] h = hist.getHistogram();
                    double[] piH = piHist.getHistogram();
                    for (int i=0; i<xValues.length; i++) {
                        if (!Double.isNaN(h[i])) {
                            System.out.println(xValues[i]+" "+h[i]+" "+piH[i]);
                        }
                    }
                    System.out.println("**** target ****");
                    xValues = targHist.xValues();
                    h = targHist.getHistogram();
                    piH = targPiHist.getHistogram();
                    for (int i=0; i<xValues.length; i++) {
                        if (!Double.isNaN(h[i])) {
                            double r = xValues[i];
                            if (r < 0) r += 1;
                            else r = Math.exp(r);
                            System.out.println(r+" "+h[i]+" "+piH[i]);
                        }
                    }
                }
            };
            sim.integratorOS.getEventManager().addListener(histReport);

            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[0].getEventManager().addListener(histListenerRef);
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }

        sim.getController().actionPerformed();
        
        if (params.doHist) {
            double[] xValues = hist.xValues();
            double[] h = hist.getHistogram();
            
            System.out.println("final ref histogram");
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
//                    System.out.println(xValues[i]+" "+(-2*h[i]+1)+" "+Math.exp(-u/temperature));
                    System.out.println(xValues[i]+" "+(-2*h[i]+1));
                }
            }
        }


        System.out.println();
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        System.out.println();
        
        sim.printResults(HSB);

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
    
    enum Level {
        CLASSICAL, SEMICLASSICAL_FH, SEMICLASSICAL_TI
    }

    /**
     * Inner class for parameters
     */
    public static class VirialCO2SCParam extends ParameterBase {
        // don't change these
        public int[] nTypes = new int[]{2,1};
        public double temperature = 100;
        public long numSteps = 10000000;
        public double refFrac = -1;
        public double sigmaHSRef = -1;
        public Level level = Level.CLASSICAL;
        public boolean doHist = false;
        public boolean nonAdditive = false;
    }
}
