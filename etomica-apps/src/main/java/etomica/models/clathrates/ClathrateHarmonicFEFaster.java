/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.clathrates;

import etomica.action.BoxImposePbc;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationFileBinary;
import etomica.data.meter.MeterPotentialEnergyFasterer;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.models.water.ConfigurationFileTIP4P;
import etomica.models.water.P2WaterTIP4P;
import etomica.models.water.SpeciesWater4P;
import etomica.nbr.cell.PotentialMasterCellFasterer;
import etomica.normalmode.LatticeDynamics;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.potential.compute.PotentialComputeEwaldFourierLD;
import etomica.potential.ewald.P2Ewald1Real;
import etomica.potential.ewald.P2Ewald6Real;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Calorie;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

public class ClathrateHarmonicFEFaster extends Simulation {
    protected static double[] initialU;
    protected Box box;
    protected PotentialCompute potentialMaster, potentialMasterLD;
    protected SpeciesGeneral species;

    public ClathrateHarmonicFEFaster(Space space, int[] nC, double rCutRealES, double[] a0_sc, int numMolecule, String configFileName, boolean isIce, double kCut, boolean includeM) {
        super(space);
        species = SpeciesWater4P.create();
        addSpecies(species);
        box = this.makeBox(new BoundaryRectangularPeriodic(space, a0_sc));
        box.setNMolecules(species, numMolecule);
        double sigma, epsilon;
        if (isIce) {
            sigma = 3.1668;
            epsilon = Kelvin.UNIT.toSim(106.1);//TIP4P/Ice
//			sigma = 3.1589; epsilon = Kelvin.UNIT.toSim(93.2);//TIP4P/2005
        } else {//TIP4P
            double A = 600E3; // kcal A^12 / mol
            double C = 610.0; // kcal A^6 / mol
            double s6 = A / C;
            sigma = Math.pow(s6, 1.0 / 6.0);
            epsilon = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C / s6 * 1000)) / 4.0;
//            sigma = 3.154; epsilon = Kelvin.UNIT.toSim(78.0); //TIP4P
        }

        AtomType oType = species.getTypeByName("O");
        AtomType hType = species.getTypeByName("H");
        AtomType mType = species.getTypeByName("M");

        double chargeH = isIce ? 0.5897 : P2WaterTIP4P.qH;
        double chargeM = -2*chargeH;

        PotentialComputeEwaldFourier ewaldFourier = new PotentialComputeEwaldFourier(getSpeciesManager(), box);
        PotentialComputeEwaldFourier.EwaldParams params = ewaldFourier.getOptimalParams(3, 0);
        params.rCut = rCutRealES;
        params.kCut = kCut;
        params.alpha = 0.30472470011002206;

        ewaldFourier.setkCut(params.kCut);
        ewaldFourier.setCharge(mType, chargeM);
        ewaldFourier.setCharge(hType, chargeH);
        ewaldFourier.setAlpha(params.alpha);
        ewaldFourier.setR6Coefficient(oType, sigma, epsilon);
        ewaldFourier.setAlpha6(params.alpha);

        PotentialComputeEwaldFourierLD ewaldFourierLD = new PotentialComputeEwaldFourierLD(getSpeciesManager(), box);
        ewaldFourierLD.setkCut(params.kCut);
        ewaldFourierLD.setCharge(mType, chargeM);
        ewaldFourierLD.setCharge(hType, chargeH);
        ewaldFourierLD.setAlpha(params.alpha);
        ewaldFourierLD.setR6Coefficient(oType, sigma, epsilon);
        ewaldFourierLD.setAlpha6(params.alpha);

        PotentialMasterCellFasterer pm = new PotentialMasterCellFasterer(getSpeciesManager(), box, 3, BondingInfo.noBonding());

        TruncationFactory tf = new TruncationFactorySimple(space, params.rCut);
        P2SoftSphere p2OO12 = new P2SoftSphere(space, sigma, 4*epsilon, 12);
        P2Ewald6Real p2OO6 = new P2Ewald6Real(sigma, epsilon, sigma, epsilon, params.alpha);
        Potential2Soft p2OO = tf.make(p2OO12, p2OO6);
        Potential2Soft p2MM = tf.make(new P2Ewald1Real(chargeM*chargeM, params.alpha));
        Potential2Soft p2HH = tf.make(new P2Ewald1Real(chargeH*chargeH, params.alpha));
        Potential2Soft p2MH = tf.make(new P2Ewald1Real(chargeH*chargeM, params.alpha));
        pm.setPairPotential(oType, oType, p2OO);
        pm.setPairPotential(mType, mType, p2MM);
        pm.setPairPotential(hType, hType, p2HH);
        pm.setPairPotential(mType, hType, p2MH);

        PotentialMasterBonding pmBonding = ewaldFourier.makeIntramolecularCorrection();

        potentialMaster = new PotentialComputeAggregate(pm, pmBonding, ewaldFourier);
        potentialMasterLD = new PotentialComputeAggregate(pm, pmBonding, ewaldFourierLD);

        if (includeM) {
            ConfigurationFile config = new ConfigurationFile(configFileName);////to duplicate with M point!
            ConfigurationFileBinary.replicate(config, box, nC, space);
        } else {
            ConfigurationFileTIP4P config = new ConfigurationFileTIP4P(configFileName, space, isIce);//to duplicate w.o. M point!
            ConfigurationFileBinary.replicate(config, box, nC, space);
        }
        BoxImposePbc.imposePBC(box);

        potentialMaster.init();

        MeterPotentialEnergyFasterer meterPotentialEnergy = new MeterPotentialEnergyFasterer(potentialMaster);
        double E = Joule.UNIT.fromSim(meterPotentialEnergy.getDataAsScalar() / numMolecule) * 1.0E-3 * Constants.AVOGADRO;
        System.out.println(" E (kJ/mol) = " + E);

    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        final int[] nC = params.nC;
        int nBasis = params.nBasis;
        double kCut = params.kCut;
        String configFile = params.configFile;
        final double[] a0 = params.a0;
        boolean waveVectorMethod = params.waveVectorMethod;
        boolean isIce = params.isIce;
        boolean includeM = params.includeM;
        double rCutRealES = params.rCutRealES;

        final double[] a0_sc = new double[]{a0[0] * nC[0], a0[1] * nC[1], a0[2] * nC[2]};

        System.out.println(" nX = " + nC[0]);
        System.out.println(" rCut = " + rCutRealES);
        System.out.println(" kCut = " + kCut);

        int numMolecule = nBasis * nC[0] * nC[1] * nC[2];
        final Space space = Space3D.getInstance();
        final ClathrateHarmonicFEFaster sim = new ClathrateHarmonicFEFaster(space, nC, rCutRealES, a0_sc, numMolecule, configFile, isIce, kCut, includeM);

        if (false) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "string");
            final DisplayBox display = new DisplayBox(sim.getController(), sim.box);
            simGraphic.add(display);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getTypeByName("H"), Color.GREEN);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getTypeByName("O"), Color.RED);
            //Sabry
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("O"), 2.0);
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("H"), 1.0);

            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).setBackgroundColor(Color.white);
            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).setBoundaryFrameColor(Color.blue);
            simGraphic.makeAndDisplayFrame();
            if (!true) {
                return;
            }
        }

        Primitive primitive = new PrimitiveOrthorhombic(space, a0[0], a0[1], a0[2]);

        LatticeDynamics ld = new LatticeDynamics(sim.getSpeciesManager(), sim.box, primitive, nBasis);
        PotentialCallbackMoleculeHessian pcHessianWV = new PotentialCallbackMoleculeHessian(sim.getSpeciesManager(), sim.box, ld);
        pcHessianWV.reset();

        long t1 = System.nanoTime();
        sim.potentialMasterLD.computeAll(true, pcHessianWV);
        sim.potentialMaster.computeAll(true);
        pcHessianWV.intramolecularCorrection(sim.potentialMaster.getForces());
        long t2 = System.nanoTime();

        Tensor[][][][] matrix = ld.getMatrix();
        System.out.println("With WV  0 0  wv 0, real");
        System.out.println(matrix[0][0][0][0]);
        System.out.println("With WV  0 0  wv 1, real");
        System.out.println(matrix[1][0][0][0]);
        System.out.println("With WV  0 0  wv 1, imag");
        System.out.println(matrix[1][0][0][1]);

        System.out.println("time: "+(t2-t1)/1e9);
    }


    public static class SimParams extends ParameterBase {
        //    	public String configFile = "config_from_paper_HHO_shiftedL_2_sI";
        public String configFile = "finalPos";
        public double rCutRealES = 14;
        public double[] a0 = new double[]{12.03, 12.03, 12.03};//sI
        public int nBasis = 46;//sI
        //		public int nBasis = 136;//sII
//		public int nBasis = 68;//sH
        public double kCut = 2.6;
        //		public double[] a0 = new double[]{17.31, 17.31, 17.31};//sII
//		public double[] a0 = new double[]{12.21,21.15, 10.14};//sH
        public boolean waveVectorMethod = true;
        public boolean isIce = false;
        public boolean includeM = true;
        int nX = 2;
        public int[] nC = new int[]{nX, 1, 1};
    }
}
