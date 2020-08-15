package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.DiameterHashByType;
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
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.*;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;

/**
 * B3,C2, dielectric virial coefficient, TIP4P water
 * two f bonds, one bond of mu_i * mu_j
 * Dec 2014
 */
public class VirialB3C2TIP4Pwater {

    public static void main(String[] args) {
        VirialParam params = new VirialParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.numSteps = 10000000;
            params.temperatureK = 500;
        }

        double temperature = Kelvin.UNIT.toSim(params.temperatureK);
        long steps = params.numSteps;
        final int nPoints = 3;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        System.out.println("TIP4P water dielectric B3, C2 at T=" + params.temperatureK + " Kelvin");
        System.out.println("temperature in sim:" + temperature);
        final double[] HSB = new double[9];
        HSB[3] = Standard.B3HS(sigmaHSRef);
        double coefficient = 4.0 * Math.PI / 9.0 / temperature;
        HSB[nPoints] = HSB[nPoints] * coefficient;
        System.out.println("B" + nPoints + "HS: " + HSB[nPoints]);
        // TIP4P water parameters
        double A_ = 600e3; // kcal A^12 / mol
        double C_ = 610; // kcal A^6 / mol
        double s6_ = A_ / C_;
        double sigmaLJ = Math.pow(s6_, 1.0 / 6.0);
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
        double offsetRatio = delta / sigmaLJ;
        System.out.println("charge on M:" + chargeM);
        System.out.println("dipole seperation: " + dipoleSeperation);
        System.out.println("mu:" + mu);
        System.out.println("seperation between LJ site and dipole:" + delta);
        System.out.println("offset ratio is:" + offsetRatio);
        System.out.println("sigmaLJ:" + sigmaLJ);
        System.out.println("epsilonLJ:" + epsilonLJ);

        double mu2Reduced = mu * mu / (sigmaLJ * sigmaLJ * sigmaLJ) / temperature;
        System.out.println("mu*:" + Math.sqrt(mu2Reduced));
        System.out.println("mu*^2:" + mu2Reduced);
        //calculate hard core
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
        fPlusBetaEU fEUTarget = new fPlusBetaEU(space, pTargetWater);// bondlist[0]
        MuProduct muTarget = new MuProduct(pTargetWater);// bondlist[1]

        int nBondTypes = 2;//MuF bond and f bond
        ClusterBonds[] clusters = new ClusterBonds[0];
        int[][][] bondList = new int[nBondTypes][][];
        ClusterAbstract targetCluster = null;
        bondList[0] = new int[][]{{0, 2}, {0, 1}}; // f bond
        bondList[1] = new int[][]{{1, 2}}; // Muf  bond

        clusters = (ClusterBonds[]) Arrays.addObject(clusters, new ClusterBonds(nPoints, bondList, false));
        targetCluster = new ClusterSum(clusters, new double[]{1.0}, new MayerFunction[]{fEUTarget, muTarget});
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
            sim.box[0].getBoundary().setBoxSize(Vector.of(size, size, size));
            sim.box[1].getBoundary().setBoxSize(Vector.of(size, size, size));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0 / size));
            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0 / size));
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            //set diameters
            DiameterHashByType diameter = new DiameterHashByType();
            diameter.setDiameter(species.getAtomType(0), 1);
            diameter.setDiameter(species.getAtomType(1), 1);
            diameter.setDiameter(species.getAtomType(2), 1);

            simGraphic.getDisplayBox(sim.box[0]).setDiameterHash(diameter);
            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
            ColorSchemeByType colorScheme = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), Color.red);
            colorScheme.setColor(sim.getSpecies(0).getAtomType(1), Color.cyan);
            colorScheme.setColor(sim.getSpecies(0).getAtomType(2), Color.yellow);

            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20);
            sim.getController().addActivity(new ActivityIntegrate2(sim.integratorOS));
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
                    if ((sim.integratorOS.getStepCount() * 10) % sim.getController().getMaxSteps() != 0) return;
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
        if (false) {/////
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }
sim.getController().runActivityBlocking(new ActivityIntegrate2(sim.integratorOS), 1000);
        if (false) {////
            double[] xValues = targHist.xValues();
            double[] h = targHist.getHistogram();
            double[] hPI = targPiHist.getHistogram();
            for (int i = 0; i < xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
                    System.out.println(xValues[i] + " " + h[i] + "  " + hPI[i]);
                }
            }
        }
        long t2 = System.currentTimeMillis();
        System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());
        sim.printResults(HSB[nPoints]);
        System.out.println("time: " + (t2 - t1) / 1000.0);
    }

    // exp(-beta*u)-1 + exp(-beta*u`(lj))*beta*u(dipole),u` located at dipole center
    public static class fPlusBetaEU implements MayerFunction, java.io.Serializable {
        private P2WaterTIP4PHardCore potentialWater;
        private Space space;

        public fPlusBetaEU(Space space, P2WaterTIP4PHardCore potentialWater) {
            this.potentialWater = potentialWater;
            this.space = space;
        }

        public double f(IMoleculeList pair, double r2, double beta) {
            double u_water = potentialWater.energy(pair);
            double x = -beta * u_water;
            //get fbond
            double f = 0.0;
            if (Math.abs(x) < 0.01) {
                f = x + x * x / 2.0 + x * x * x / 6.0 + x * x * x * x / 24.0 + x * x * x * x * x / 120.0;
            } else {
                f = Math.exp(x) - 1.0;
            }
            // get exp(-beta u(LJ) ) at dipole center
            double uLJAtDipoleCenter = potentialWater.getULJAtDipoleCenter();

            double y = -beta * uLJAtDipoleCenter;
            double e_bond_LJAtDipleCenter = 0;

            if (Math.abs(y) < 0.01) {
                e_bond_LJAtDipleCenter = 1.0 + y + y * y / 2.0 + y * y * y / 6.0 + y * y * y * y / 24.0 + y * y * y * y * y / 120.0;
            } else {
                e_bond_LJAtDipleCenter = Math.exp(y);
            }
            //get u(dipole)
            double u_DipoleDipole = potentialWater.getU_DipoleDipole();
            double sum = f + beta * e_bond_LJAtDipleCenter * u_DipoleDipole;
//            System.out.println("uLJAtDipoleCenter:"+uLJAtDipoleCenter);
//            System.out.println("u_DipoleDipole:"+u_DipoleDipole);
//            System.out.println("sum:"+sum);          
            return sum;
        }

        /* (non-Javadoc)
         * @see etomica.virial.MayerFunction#getPotential()
         */
        public IPotential getPotential() {
            // TODO Auto-generated method stub
            return null;
        }

        public void setBox(Box newBox) {

        }
    }

    public static class MuProduct implements MayerFunction, java.io.Serializable {

        private P2WaterTIP4PHardCore potentialWater;

        public MuProduct(P2WaterTIP4PHardCore potentialWater) {
            this.potentialWater = potentialWater;
        }

        public double f(IMoleculeList pair, double r2, double beta) {
            potentialWater.energy(pair);
            Vector dipole1 = potentialWater.getDipole1();
            Vector dipole2 = potentialWater.getDipole2();
            double mu_Dot_mu = dipole1.dot(dipole2);
            return mu_Dot_mu;
        }

        public IPotential getPotential() {
            return null;
        }

        public void setBox(Box newBox) {
        }
    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        public double temperatureK = 350;
        public long numSteps = 100000000;
        public double sigmaHSRef = 7.0;
        public double refFrac = -1;
    }
}
