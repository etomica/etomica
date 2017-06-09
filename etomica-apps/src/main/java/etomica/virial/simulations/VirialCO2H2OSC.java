/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.api.ISpecies;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.models.co2.P2CO2H2OWheatley;
import etomica.models.co2.P2CO2Hellmann;
import etomica.models.water.P2WaterSzalewicz;
import etomica.models.water.P2WaterSzalewicz.Component;
import etomica.potential.*;
import etomica.potential.P2SemiclassicalAtomic.AtomInfo;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

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
            params.level = Level.SEMICLASSICAL_FH;
            params.temperature = 300;
            params.numSteps = 1000000;
            params.sigmaHSRef = 5;
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
        final Nonadditive nonAdditive = nPoints < 3 ? Nonadditive.NONE : params.nonAdditive;
        final boolean useSZ = nonAdditive != Nonadditive.NONE && params.useSZ && nTypes[1]==nPoints;

        final double HSB = Standard.BHS(nPoints, sigmaHSRef);

        System.out.println("Overlap sampling for CO2/H2O mixture at " + temperatureK + " K");
        System.out.println("nTypes: "+java.util.Arrays.toString(nTypes));
        if (nonAdditive != Nonadditive.NONE) {
            if (nonAdditive == Nonadditive.DISPERSION) {
                System.out.println("Including non-additive dispersion");
            }
            else {
                System.out.println("Including non-additive dispersion and induction");
            }
            if (useSZ) {
                System.out.println("    as prescribed by Gora...Szalewicz");
            }
        }
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
        AtomTypeAgentManager paramsManagerATM = null, paramsManagerInd = null;
        if (nonAdditive != Nonadditive.NONE) {
            if (useSZ && nPoints == nTypes[1]) {
                Component comp = nonAdditive == Nonadditive.FULL ? Component.NON_PAIR : Component.THREE_BODY;
                P2WaterSzalewicz p23cH2O = new P2WaterSzalewicz(space, 2);
                p23cH2O.setComponent(comp);
                IPotentialAtomic p23aH2O = level == Level.CLASSICAL ? null : (level == Level.SEMICLASSICAL_FH ? p23cH2O.makeSemiclassical(temperature) : new P2SemiclassicalAtomic(space, p23cH2O, temperature));
                PotentialMolecularMonatomic p23H2O = new PotentialMolecularMonatomic(space, level==Level.CLASSICAL ? p23cH2O : p23aH2O);

                P2WaterSzalewicz p3cH2O = new P2WaterSzalewicz(space, 3);
                p3cH2O.setComponent(comp);
                IPotentialAtomic p3aH2O = level == Level.CLASSICAL ? null : (level == Level.SEMICLASSICAL_FH ? p3cH2O.makeSemiclassical(temperature) : null);
                PotentialMolecularMonatomic p3H2O = new PotentialMolecularMonatomic(space, level==Level.CLASSICAL ? p3cH2O : p3aH2O);
                MayerFunctionMolecularThreeBody f3H2O = new MayerFunctionMolecularThreeBody(new PotentialNonAdditive(new IPotentialMolecular[]{p23H2O,p3H2O}));
                MayerFunctionNonAdditive[][][] allFNA = new MayerFunctionNonAdditive[2][2][2];
                allFNA[1][1][1] = f3H2O;
                targetCluster = new ClusterWheatleyMultibodyMix(nPoints, nTypes, allF, allFNA, 1e-12, true);
                ((ClusterWheatleyMultibodyMix)targetCluster).setDoCaching(false);
                targetCluster = new ClusterCoupledAtomFlipped(targetCluster, space, 20);
            }
            else {
                // CO2: alpha=2.913, E=13.7eV
                // H2O: alpha=1.444, E=12.6eV
                paramsManagerATM = new AtomTypeAgentManager(null);
                IPotentialAtomic p3 = null;
                p3 = new P3AxilrodTeller(space, paramsManagerATM);
                if (nonAdditive == Nonadditive.FULL) {
                    paramsManagerInd = new AtomTypeAgentManager(null);
                    p3 = new PotentialAtomicSum(new IPotentialAtomic[]{p3,new P3Induction(space, paramsManagerInd)});
                }

                MayerFunctionMolecularThreeBody f3 = new MayerFunctionMolecularThreeBody(new PotentialMolecularMonatomic(space, p3));
                MayerFunctionNonAdditive[][][] allFNA = new MayerFunctionNonAdditive[2][2][2];
                allFNA[0][0][0] = f3;
                allFNA[0][0][1] = f3;
                allFNA[0][1][0] = f3;
                allFNA[0][1][1] = f3;
                allFNA[1][0][0] = f3;
                allFNA[1][0][1] = f3;
                allFNA[1][1][0] = f3;
                allFNA[1][1][1] = f3;
                targetCluster = new ClusterWheatleyMultibodyMix(nPoints, nTypes, allF, allFNA, 1e-12, true);
                if (nonAdditive == Nonadditive.FULL && nTypes[1] > 0) {
                    ((ClusterWheatleyMultibodyMix)targetCluster).setDoCaching(false);
                    targetCluster = new ClusterCoupledAtomFlipped(targetCluster, space, 15);
                }
            }
        }

        
        targetCluster.setTemperature(temperature);

        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);

        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
 		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{speciesCO2,speciesH2O}, nTypes, temperature, refCluster, targetCluster);
        sim.init();
        int[] seeds = sim.getRandomSeeds();
        System.out.println("Random seeds: "+Arrays.toString(seeds));
        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        if (nonAdditive != Nonadditive.NONE && !useSZ) {
            double alphaCO2 = 2.913;
            double alphaH2O = 1.444;
            if (nonAdditive == Nonadditive.FULL) {
                Vector polCO2 = space.makeVector();
                double[] qCO2 = new double[7];
                Vector[] qSiteCO2 = new Vector[7];
                for (int i=0; i<7; i++) {
                    qCO2[i] = p2cCO2.getQ(i);
                    Vector r = space.makeVector();
                    r.setX(0, p2cCO2.getPos(i));
                    qSiteCO2[i] = r;
                }

                P3Induction.MyAgent agentCO2 = new P3Induction.MyAgent(new double[]{alphaCO2}, new Vector[]{polCO2}, qCO2, qSiteCO2);

                Vector polH2O = space.makeVector();
                double[] qH2O = P2WaterSzalewicz.getQ();
                Vector[] qSiteH2O = P2WaterSzalewicz.getSites(space);
                polH2O.E(qSiteH2O[0]);
                P3Induction.MyAgent agentH2O = new P3Induction.MyAgent(new double[]{alphaH2O}, new Vector[]{polH2O}, qH2O, qSiteH2O);

                paramsManagerInd.setAgent(speciesCO2.getLeafType(), agentCO2);
                paramsManagerInd.setAgent(speciesH2O.getLeafType(), agentH2O);
            }
            paramsManagerATM.setAgent(speciesCO2.getLeafType(), new P3AxilrodTeller.MyAgent(alphaCO2, ElectronVolt.UNIT.toSim(13.7)));
            paramsManagerATM.setAgent(speciesH2O.getLeafType(), new P3AxilrodTeller.MyAgent(alphaH2O, ElectronVolt.UNIT.toSim(12.6)));
        }

        if (level == Level.SEMICLASSICAL_TI && false) {
            if (true) throw new RuntimeException("implement me for anything other than CO2-CO2");
            final Vector[] rv = new Vector[4];
            for (int i=0; i<4; i++) {
                rv[i] = space.makeVector();
            }
            double om = Oxygen.INSTANCE.getMass();
            double bondLength = p2cCO2.posB[5] - p2cCO2.posB[1];
            rv[0].setX(0, om*bondLength*bondLength*0.25);
            rv[0].setX(1, om*bondLength*bondLength*0.25);
            ((P2SemiclassicalAtomic)p2aCO2).setAtomInfo(speciesCO2.getLeafType(), new AtomInfo() {
                public Vector[] getMomentAndAxes(IAtomOriented molecule) {
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
        
        if (nonAdditive != Nonadditive.NONE) {
            double r = 3;
            IAtomList atoms = sim.box[1].getLeafList();
            for (int i=1; i<nPoints; i++) {
                Vector pos = atoms.getAtom(i).getPosition();
                double theta = 2*i*Math.PI/nPoints;
                pos.setX(0, r*(1-Math.cos(theta)));
                pos.setX(1, r*Math.sin(theta));
            }
        }
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{40,40,40}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{40,40,40}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
//            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);

//            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
//            displayBox0.setColorScheme(colorScheme);
//            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
//            displayBox1.setColorScheme(colorScheme);
//            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 10);
                    sim.equilibrate(null, 20);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            simGraphic.makeAndDisplayFrame();
            return;
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
            refFileName = "refpref"+nPoints+"_"+(nonAdditive==Nonadditive.NONE?"2":("3"+(nonAdditive==Nonadditive.DISPERSION?"a":"ai")))+(useSZ?"_sz":"")+"_"+tempString;
            refFileName += level == Level.CLASSICAL ? "_c" : (level == Level.SEMICLASSICAL_FH ? "_fh" : "_ti");
            refFileName += "C";
        }

        sim.initRefPref(refFileName, steps/40);
        sim.equilibrate(refFileName, steps/20);
        
        System.out.println("equilibration finished");
        
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }

        if (params.doHist) {
            sim.setupTargetHistogram();
        }

        sim.getController().actionPerformed();
        
        if (params.doHist) {
            sim.printTargetHistogram();
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
    enum Nonadditive {
        NONE, DISPERSION, FULL
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
        public Nonadditive nonAdditive = Nonadditive.NONE;
        public boolean useSZ = false;
    }
}
