package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.histogram.HistogramSimple;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.models.water.P2WaterTIP4PHardCore;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.potential.PotentialMolecular;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.*;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * for water TIP4P
 * B2, dielectric constant
 * exp (-beta u(water) ) - exp(-beta u(DHS) )
 */
public class VirialB2WaterTIP4P_DHS_difference {

    public static void main(String[] args) {
        VirialParam params = new VirialParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.numSteps = 100000000;
            params.temperatureK = 500;
        }
        double temperature = Kelvin.UNIT.toSim(params.temperatureK);
        long steps = params.numSteps;
        final int nPoints = 2;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        System.out.println("TIP4P water dielectric B2 - dipolar HS at T=" + params.temperatureK + " Kelvin");
        System.out.println("temperature in sim:" + temperature);
        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        System.out.println("sigmaHSRef: " + sigmaHSRef);
        double coefficient = 4 * Math.PI / 9 / temperature;

        HSB[nPoints] = HSB[nPoints] * coefficient;
        System.out.println("B" + nPoints + "HS: " + HSB[nPoints]);
        double A_ = 600e3; // kcal A^12 / mol
        double C_ = 610; // kcal A^6 / mol
        double s6_ = A_ / C_;
        double sigmaLJ = Math.pow(s6_, 1.0 / 6.0);
        double sigmaHS = sigmaLJ;
        double epsilonLJ = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C_ / s6_ * 1000)) / 4.0;
        double OM = 0.15;
        double angle = 104.52 * Math.PI / 180.;
        double OH = 0.9572;
        double delta = OM + (OH * Math.cos(angle / 2) - OM) / 2.0;
        double chargeM = Electron.UNIT.toSim(-1.04);
        double chargeH = Electron.UNIT.toSim(0.52);
        double dipoleSeperation = OH * Math.cos(angle / 2) - OM;
        double mu = -chargeM * dipoleSeperation;//dipole moment in simulation unit
        double mu2 = mu * mu;
        System.out.println("charge on M:" + chargeM);
        System.out.println("dipole seperation: " + dipoleSeperation);
        System.out.println("mu:" + mu);
        System.out.println("seperation between LJ site and dipole:" + delta);
        System.out.println("sigmaHS:" + sigmaHS);
        System.out.println("sigmaLJ:" + sigmaLJ);
        System.out.println("epsilonLJ:" + epsilonLJ);
        double mu2Reduced = mu * mu / (sigmaLJ * sigmaLJ * sigmaLJ) / temperature;

        System.out.println("mu*:" + Math.sqrt(mu2Reduced));
        System.out.println("mu*^2:" + mu2Reduced);
        double x4 = 6 * mu2 / 48.0 / Math.pow(sigmaLJ, 12) / epsilonLJ;
        double x = Math.pow(x4, 1.0 / 4);
        double data = 2 * delta;
        System.out.println("data[0]: " + data);
        double hardCore = 0;
        for (int i = 1; i < 10; i++) {
            data = 2 * delta + Math.pow(data, 13.0 / 4) * x;
        }
        hardCore = data;
        System.out.println("hard core:" + hardCore);
        Space space = Space3D.getInstance();

        //HS ref system
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints > 3, eRef, true);
        refCluster.setTemperature(temperature);

        //Target system
        P2WaterTIP4PHardCore pTargetWater = new P2WaterTIP4PHardCore(space, hardCore, sigmaLJ, epsilonLJ, chargeM, chargeH);
        P2WaterDHS pTargetDHS = new P2WaterDHS(space, sigmaHS, mu);

        MuDeltaEBond muDeltaETarget = new MuDeltaEBond(space, pTargetWater, pTargetDHS);

        int nBondTypes = 1;
        ClusterBonds[] clusters = new ClusterBonds[0];
        int[][][] bondList = new int[nBondTypes][][];
        ClusterAbstract targetCluster = null;

        bondList[0] = new int[][]{{0, 1}};
        clusters = (ClusterBonds[]) Arrays.addObject(clusters, new ClusterBonds(nPoints, bondList, false));
        targetCluster = new ClusterSum(clusters, new double[]{1.0}, new MayerFunction[]{muDeltaETarget});
        //     targetCluster = new ClusterCoupledAtomFlipped(targetCluster, space);

        targetCluster.setTemperature(temperature);

        double refIntegral = HSB[nPoints];

        SpeciesWater4P species = new SpeciesWater4P(space);
        //simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, temperature, refCluster, targetCluster, false);

        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps + " steps (1000 blocks of " + steps / 1000 + ")");
        steps /= 1000;
        // displace the atoms to a certain distance
        for (int i = 0; i < 50 && sim.box[1].getSampleCluster().value(sim.box[1]) == 0; i++) {
            sim.mcMoveTranslate[1].doTrial();
            sim.mcMoveTranslate[1].acceptNotify();
            sim.box[1].trialNotify();
            sim.box[1].acceptNotify();
        }

        if (false) {
            double size = 20;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size, size, size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0 / size));
            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0 / size));
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            //set diameters
            DiameterHashByType diameter = new DiameterHashByType();
            diameter.setDiameter(species.getHydrogenType(), 1);
            diameter.setDiameter(species.getOxygenType(), 1);
            diameter.setDiameter(species.getMType(), 1);

            simGraphic.getDisplayBox(sim.box[0]).setDiameterHash(diameter);
            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
            ColorSchemeByType colorScheme = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
            colorScheme.setColor(species.getHydrogenType(), Color.red);
            colorScheme.setColor(species.getOxygenType(), Color.yellow);
            colorScheme.setColor(species.getMType(), Color.black);

            simGraphic.makeAndDisplayFrame();

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
            if (Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0) {
                throw new RuntimeException("Oops");
            }

            return;
        }
        long t1 = System.currentTimeMillis();
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref" + nPoints + "_" + temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps / 40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps / 20);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }

        System.out.println("equilibration finished");

        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int) steps);
        sim.ai.setMaxSteps(1000);
        sim.integratorOS.getMoveManager().setEquilibrating(false);

        for (int i = 0; i < 2; i++) {
            System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize() + " " + sim.mcMoveRotate[i].getStepSize());
        }
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        if (false) {
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {
                }

                public void integratorStepStarted(IntegratorEvent e) {
                }

                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount() * 10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount() + " steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    System.out.println("abs average: " + ratioAndError[0] * HSB[nPoints] + ", error: " + ratioAndError[1] * HSB[nPoints]);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }
        final HistogramSimple targHist = new HistogramSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple targPiHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();

        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {
            }

            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                double r2Min = Double.POSITIVE_INFINITY;
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i = 0; i < nPoints; i++) {
                    for (int j = i + 1; j < nPoints; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij < r2Min) r2Min = r2ij;
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }

                double v = finalTargetCluster.value(sim.box[1]);
                double r = Math.sqrt(r2Max);
                if (r > 1) {
                    r = Math.log(r);
                } else {
                    r -= 1;
                }
                targHist.addValue(r);
                targPiHist.addValue(r, Math.abs(v));
            }

            public void integratorInitialized(IntegratorEvent e) {
            }
        };
        if (false) {
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }
        sim.getController().actionPerformed();
        if (false) {
            double[] xValues = targHist.xValues();
            double[] h = targHist.getHistogram();
            double[] hPI = targPiHist.getHistogram();
            for (int i = 0; i < xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
                    if (xValues[i] > 0) {
                        System.out.println(Math.exp(xValues[i]) + " " + h[i] + "  " + hPI[i]);
                    } else {
                        System.out.println((xValues[i] + 1) + " " + h[i] + "  " + hPI[i]);
                    }
                }
            }
        }
        long t2 = System.currentTimeMillis();

        System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());
        sim.printResults(HSB[nPoints]);// s is the abs average
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    //inner class
    //for dielectric virial calculation
    // mu_i * mu_j * (ebond(P2TIP4Pwater, offset) - ebond(DHS))

    public static class MuDeltaEBond implements MayerFunction, java.io.Serializable {

        private P2WaterTIP4PHardCore potentialWater;
        private P2WaterDHS potentialWaterDHS;
        private Vector groupTranslationVector;
        private Vector O1O2;
        private Space space;

        public MuDeltaEBond(Space _space, P2WaterTIP4PHardCore potentialWater, P2WaterDHS potentialWaterDHS) {
            this.space = _space;
            this.potentialWater = potentialWater;
            this.potentialWaterDHS = potentialWaterDHS;
            groupTranslationVector = space.makeVector();
            O1O2 = space.makeVector();
        }

        public double f(IMoleculeList pair, double r2, double beta) {
            double u_water = potentialWater.energy(pair);
            double ebond_water = Math.exp(-beta * u_water);
            double u_waterDHS = potentialWaterDHS.energy(pair);
            Vector dipole1 = potentialWaterDHS.getDipole1();
            Vector dipole2 = potentialWaterDHS.getDipole2();
            double ebond_waterDHS = Math.exp(-beta * u_waterDHS);
            double ebond_difference = ebond_water - ebond_waterDHS;
            double dipole1_dot_dipole2 = dipole1.dot(dipole2);
            double cos12 = dipole1_dot_dipole2 / dipole1.squared();

//			  System.out.println("MuDeltaEBond class, u_water: "+u_water);
//              System.out.println("MuDeltaEBond class, u_waterDHS: "+u_waterDHS); 
//              System.out.println("MuDeltaEBond class, ebond_water: "+ebond_water);
//              System.out.println("MuDeltaEBond class, ebond_waterDHS: "+ebond_waterDHS);


            boolean debug = false;
            if (debug && r2 > 40000) {
                ///// ///////get O site position
                IMolecule water1 = pair.getMolecule(0);
                IMolecule water2 = pair.getMolecule(1);
                IAtomList atomList1 = water1.getChildList();
                IAtomList atomList2 = water2.getChildList();
                Vector O1 = atomList1.get(2).getPosition();
                Vector O2 = atomList2.get(2).getPosition();
                groupTranslationVector.Ev1Mv2(O2, O1);
                groupTranslationVector.normalize();
                for (int m = 0; m < 200; m++) {
                    double atomCount = atomList2.size();
                    for (int q = 0; q < atomCount; q++) {
                        IAtom atom = atomList2.get(q);
                        atom.getPosition().PE(groupTranslationVector);
                    }
                    O1O2.Ev1Mv2(O2, O1);
                    double distance = Math.sqrt(O1O2.squared());
                    double u_water_ = potentialWater.energy(pair);
                    double ebond_water_ = Math.exp(-beta * u_water_);
                    double u_waterDLJ_ = potentialWaterDHS.energy(pair);
                    double ebond_waterDLJ_ = Math.exp(-beta * u_waterDLJ_);
                    double ebond_difference_ = ebond_water_ - ebond_waterDLJ_;
                    System.out.println(distance + "  " + ebond_difference_);
                }
                System.exit(0);
            }

            return dipole1_dot_dipole2 * ebond_difference;
        }

        @Override
        public IPotential getPotential() {
            // TODO Auto-generated method stub
            return null;
        }

        @Override
        public void setBox(Box box) {
            // TODO Auto-generated method stub
        }

    }

    public static class P2WaterDHS extends PotentialMolecular {

        public double sigma, sigmaHS2;
        double mu;
        protected Boundary boundary;
        private Vector P1r, P2r;
        private Vector O1P1r, O2P2r;
        private Vector P1P2r;// vector between two dipoles
        private Vector mu1Normalized, mu2Normalized;
        private Vector dipole1, dipole2;

        public P2WaterDHS(Space space, double sigma, double mu) {
            super(2, space);
            this.sigma = sigma;
            this.mu = mu;
            sigmaHS2 = sigma * sigma;
            P1r = space.makeVector();
            P2r = space.makeVector();
            O1P1r = space.makeVector();
            O2P2r = space.makeVector();
            P1P2r = space.makeVector();
            mu1Normalized = space.makeVector();
            mu2Normalized = space.makeVector();
            dipole1 = space.makeVector();
            dipole2 = space.makeVector();
        }

        public void setBox(Box box) {
            boundary = box.getBoundary();
        }

        public double energy(IMoleculeList pair) {
            IMolecule water1 = pair.getMolecule(0);
            IMolecule water2 = pair.getMolecule(1);

            Vector O1r = (water1.getChildList().get(2)).getPosition();//H-H-O-M, so O is the third atom
            Vector O2r = (water2.getChildList().get(2)).getPosition();
            Vector H11r = (water1.getChildList().get(0)).getPosition();
            Vector H12r = (water1.getChildList().get(1)).getPosition();
            Vector H21r = (water2.getChildList().get(0)).getPosition();
            Vector H22r = (water2.getChildList().get(1)).getPosition();

            Vector M1r = water1.getChildList().get(3).getPosition();
            Vector M2r = water2.getChildList().get(3).getPosition();

            //get position of p1 and p2, the dipole position, also the HS center
            double offset = 0.36794113830914743;
            O1P1r.Ev1Mv2(M1r, O1r);
            O1P1r.normalize();
            mu1Normalized.E(O1P1r);
            O1P1r.TE(offset);
            //get P1
            P1r.Ev1Pv2(O1P1r, O1r);

            O2P2r.Ev1Mv2(M2r, O2r);
            O2P2r.normalize();
            mu2Normalized.E(O2P2r);
            O2P2r.TE(offset);
            //get P2
            P2r.Ev1Pv2(O2P2r, O2r);
            //get P1P2r, vector between dipoles, also between lj sites
            P1P2r.Ev1Mv2(P2r, P1r);
            // store vector between two dipoles in vDipoles
            boolean doTest = false;
            if (doTest) {
                Vector test = new Vector3D();
                System.out.println("====== test molecule 1 ====== ");
                test.Ev1Mv2(P1r, H11r);
                System.out.println("P1H11:" + Math.sqrt(test.squared()));
                test.Ev1Mv2(P1r, H12r);
                System.out.println("P1H12:" + Math.sqrt(test.squared()));
                test.Ev1Mv2(P1r, M1r);
                System.out.println("P1M:" + Math.sqrt(test.squared()));
                test.Ev1Mv2(P1r, O1r);
                System.out.println("P1O1:" + Math.sqrt(test.squared()));
                test.Ev1Mv2(O1r, M1r);
                System.out.println("O1M1:" + Math.sqrt(test.squared()));
                System.out.println("====== test molecule 2 ====== ");
                test.Ev1Mv2(P2r, H21r);
                System.out.println("P2H21:" + Math.sqrt(test.squared()));
                test.Ev1Mv2(P2r, H22r);
                System.out.println("P2H2:" + Math.sqrt(test.squared()));
                test.Ev1Mv2(P2r, M2r);
                System.out.println("P2M:" + Math.sqrt(test.squared()));
                test.Ev1Mv2(P2r, O2r);
                System.out.println("P2O2:" + Math.sqrt(test.squared()));
                System.out.println("mu1normalized:" + mu1Normalized.squared());
            }
            double r2 = P1P2r.squared();

            if (r2 < sigmaHS2) {
                return Double.POSITIVE_INFINITY;
            }
            //boolean printme = r2 > 400;
            P1P2r.normalize();
            double r = Math.sqrt(r2);//distance between two dipoles

            //use the other formula to get dipole-dipole potential
            double u_dipole = mu * mu * (mu1Normalized.dot(mu2Normalized) - 3.0 * mu1Normalized.dot(P1P2r) * P1P2r.dot(mu2Normalized)) / r / r2;
            //   if (printme) System.out.println("P2WaterDLJ,u_dipole:"+u_dipole);
            //  if (printme) System.out.println("P2WaterDLJ, uLJ:"+ulj);
            return u_dipole;
        }

        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

        public Vector getDipole1() {
            dipole1.E(mu1Normalized);
            dipole1.TE(mu);
            return dipole1;
        }

        public Vector getDipole2() {
            dipole2.E(mu2Normalized);
            dipole2.TE(mu);
            return dipole2;
        }

        public double getSigma() {
            return sigma;
        }

    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        public double temperatureK = 680;
        public long numSteps = 100000000;
        public double sigmaHSRef = 7.0;
        public double refFrac = -1;
    }
}
