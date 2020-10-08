/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.clathrates;

import etomica.action.MoleculeActionTranslateTo;
import etomica.action.WriteConfiguration;
import etomica.atom.*;
import etomica.box.Box;
import etomica.box.storage.Tokens;
import etomica.box.storage.VectorStorage;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationFileBinary;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.models.water.ConfigurationFileTIP4P;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.normalmode.MeterHarmonicEnergy;
import etomica.potential.*;
import etomica.potential.EwaldSummation.MyCharge;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

public class MinimizationTIP4P extends Simulation {
    private static final long serialVersionUID = 1L;
    protected static double[] initialU;
    protected static double selfELJ;
    protected Box box;
    protected PotentialMaster potentialMaster;
    protected SpeciesGeneral species;
    protected MeterHarmonicEnergy meterHarm;
    protected MeterPotentialEnergy meterPE;
    protected Potential2SoftSphericalLS potentialLJLS;
    protected EwaldSummation potentialES;

    public MinimizationTIP4P(Space space, double rCutLJ, double rCutRealES, double[] a0, int[] nC, int nBasis, boolean isIce, double kCut, String configFile, boolean includeM) {
        super(space);
        species = SpeciesWater4P.create();
        addSpecies(species);
        double[] a0_sc = new double[]{a0[0] * nC[0], a0[1] * nC[1], a0[2] * nC[2]};
        Boundary boundary = new BoundaryRectangularPeriodic(space, a0_sc);
        box = this.makeBox(boundary);
        box.setNMolecules(species, nBasis * nC[0] * nC[1] * nC[2]);
        ChargeAgentSourceRPM agentSource = new ChargeAgentSourceRPM(species, isIce);
        AtomLeafAgentManager<MyCharge> atomAgentManager = new AtomLeafAgentManager<MyCharge>(agentSource, box);
        double sigma, epsilon;
        if (isIce) {
            sigma = 3.1668;
            epsilon = Kelvin.UNIT.toSim(106.1);//TIP4P/Ice
        } else {//TIP4P
            double A = 600E3; // kcal A^12 / mol
            double C = 610.0; // kcal A^6 / mol
            double s6 = A / C;
            sigma = Math.pow(s6, 1.0 / 6.0);
            epsilon = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C / s6 * 1000)) / 4.0;
        }
        P2LennardJones potentialLJ = new P2LennardJones(space, sigma, epsilon);

        potentialLJLS = new Potential2SoftSphericalLS(space, rCutLJ, a0_sc, potentialLJ);

        potentialES = new EwaldSummation(box, atomAgentManager, space, kCut, rCutRealES);
        potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(potentialLJLS, new AtomType[]{species.getTypeByName("O"), species.getTypeByName("O")});
        potentialMaster.addPotential(potentialES, new AtomType[0]);
        potentialLJLS.setBox(box);

        if (includeM) {
            ConfigurationFile config = new ConfigurationFile(configFile);
            ConfigurationFileBinary.replicate(config, box, nC, space);
        } else {
            ConfigurationFileTIP4P config = new ConfigurationFileTIP4P(configFile, space, isIce);
            ConfigurationFileBinary.replicate(config, box, nC, space);
        }

//Wrap all MOLECULES (make O in the BOX)
        if (!true) {
            for (int i = 0; i < box.getMoleculeList().size(); i++) {
                IMolecule molecule = box.getMoleculeList().get(i);
                IAtomList childList = molecule.getChildList();
                Vector O = childList.get(2).getPosition();//O
                O.PE(boundary.centralImage(O));// to wrap all O inside the BOX; next steps will move Hs and M with O to keep the conformation.
                for (int iChild = 0; iChild < childList.size(); iChild++) {
                    if (iChild == 2) {//Ignore O
                        continue;
                    }
                    IAtom a = childList.get(iChild);
                    Vector ri = a.getPosition();
                    ri.ME(O);
                    boundary.nearestImage(ri);
                    ri.PE(O);
                }
            }
        }

        AtomPair selfAtomLJ = new AtomPair(box.getLeafList().get(2), box.getLeafList().get(2));
        selfELJ = potentialLJLS.energy(selfAtomLJ);// = 0 if nLJshells=0
        System.out.println("selfELJ = " + Joule.UNIT.fromSim(selfELJ) * 1.0E-3 * Constants.AVOGADRO);
        MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(potentialMaster, box);
        double E = Joule.UNIT.fromSim(meterPotentialEnergy.getDataAsScalar() / box.getMoleculeList().size() + selfELJ) * 1.0E-3 * Constants.AVOGADRO;
        System.out.println("E (kJ/mol)  = " + E);
        System.out.println("E (kCal/mol)  = " + E * 0.239005736);
        System.out.println("");
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        String configFileName = params.configFile;
        int[] nC = params.nC;
        int nBasis = params.nBasis;
        double[] a0 = params.a0;
        double rCutRealES = params.rCutRealES;
        double rCutLJ = params.rCutLJ;
        int nOuter = params.nOuter;
        int nInner = params.nInner;
        double kCut = params.kCut;
        boolean isIce = params.isIce;
        boolean isGraphics = params.isGraphics;
        boolean includeM = params.includeM;

        System.out.println(isIce ? "TIP4P/Ice" : "TIP4P");
        System.out.println(" nX = " + nC[0]);
        System.out.println(" nBasis = " + nBasis);
        System.out.println(" rCutLJ = " + rCutLJ + " Angs");
        System.out.println(" rCutRealES = " + rCutRealES + "Angs");
        System.out.println(" kCut = " + kCut + " Angs^-1");

        int ng = 3 * nBasis;
        int nd = 4 * nBasis;
        double[] d = new double[nd];
        double[] x0 = new double[ng];
        final MinimizationTIP4P sim = new MinimizationTIP4P(Space3D.getInstance(3), rCutLJ, rCutRealES, a0, nC, nBasis, isIce, kCut, configFileName, includeM);
        final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
        VectorStorage forces = sim.box.getAtomStorage(Tokens.FORCES);
        PotentialCalculationForceSum pcForce = new PotentialCalculationForceSum(forces);

        MoleculeActionTranslateTo translator = new MoleculeActionTranslateTo(sim.space);
        MoleculePositionGeometricCenter pos = new MoleculePositionGeometricCenter(sim.space);

        translator.setAtomPositionDefinition(pos);

        Vector p = sim.space.makeVector();
        double step1 = 0;
        Vector[] orient0 = new Vector[nBasis];
        Vector[] orientf = new Vector[nBasis];
        Vector[] torques = new Vector[nBasis];
        Vector[] axes = new Vector[nBasis];

        for (int i = 0; i < nBasis; i++) {
            axes[i] = sim.space.makeVector();
            torques[i] = sim.space.makeVector();
            orient0[i] = sim.space.makeVector();
            orientf[i] = sim.space.makeVector();
        }
        Vector f = sim.space.makeVector();
        Vector f4 = sim.space.makeVector();
        Vector dr = sim.space.makeVector();

        for (int outer = 0; outer < nOuter; outer++) {
            System.out.println();
            System.out.println("****************************** " + outer + " ******************************");

            double totalD = 0;
            double step = 0;
            double radianFac = 1;
            RotationTensor3D rTensor = new RotationTensor3D();

            for (int iter = 0; iter < nInner; iter++) {
                double[] g = new double[ng];
                double t = 0;
                IteratorDirective id = new IteratorDirective();
                id.includeLrc = false;
                sim.potentialMaster.calculate(sim.box, id, pcForce);
                latticeEnergy = Joule.UNIT.fromSim(meterPotentialEnergy.getDataAsScalar() / nBasis + selfELJ);
                latticeEnergy *= (1E-3 * Constants.AVOGADRO);
                System.out.println("latticeEnergy/N (kJ/mol) = " + latticeEnergy);
                System.out.println("latticeEnergy/N (kCal/mol) = " + latticeEnergy * 0.239005736);
                for (int i = 0; i < nBasis; i++) {
                    torques[i].E(0);
                    IMolecule iMol = sim.box.getMoleculeList().get(i);
                    IAtomList atoms = iMol.getChildList();
                    Vector pO = atoms.get(SpeciesWater4P.indexO).getPosition();
                    for (int j = 0; j < atoms.size(); j++) {
                        IAtom atomj = atoms.get(j);
                        Vector fj = forces.get(atomj);
                        f.PE(fj);
                        if (j != SpeciesWater4P.indexO) {
                            dr.Ev1Mv2(atomj.getPosition(), pO);
                            sim.box.getBoundary().nearestImage(dr);
                            dr.XE(fj);
                            torques[i].PE(dr);
                        }
                    }//End j
//                    pcForce.reset();

                    for (int j = 0; j < 3; j++) {
                        g[i * 3 + j] = -f.getX(j);
                        t += g[i * 3 + j] * g[i * 3 + j];
                    }
                    f.E(0);

                    if (iter == 0) {
                        double t2 = torques[i].squared();
                        double sqrtT = Math.sqrt(t2);
                        t += t2 / (radianFac * radianFac);
                        d[ng + i] = sqrtT / radianFac;
                        axes[i].E(torques[i]);
                        axes[i].TE(1 / sqrtT);
                        if (torques[i].isZero()) {
                            axes[i].E(new double[]{1.0, 0, 0});
                        }
                    }
                }//End i
                pcForce.reset();
                if (iter == 0) {
                    t = Math.sqrt(t);
                    for (int j = 0; j < ng; j++) {
                        d[j] = -g[j] / t;
                        totalD += g[j] * d[j];
                    }
                    for (int j = ng; j < nd; j++) {
                        d[j] /= t;
                        totalD -= d[j] * torques[j - ng].dot(axes[j - ng]) / radianFac;
                    }
                    System.out.println("totalD " + totalD);
                    if (step1 == 0) {
                        step = 0.0000001;
                    } else {
                        step = step1 / 5;
                    }
                    System.out.println(iter + ">>>>step " + step);
                } else {
                    double newTotalD = 0;
                    for (int j = 0; j < ng; j++) {
                        newTotalD += g[j] * d[j];
                    }
                    for (int j = ng; j < nd; j++) {
                        newTotalD -= d[j] * torques[j - ng].dot(axes[j - ng]) / radianFac;
                    }

                    double oldStep = step;
                    if (newTotalD > totalD || totalD > 0.0) {
                        step = -newTotalD * step / (newTotalD - totalD);
                        if (iter == 1) {
                            step1 = oldStep + step;
                        }
                    } else {

                        step = step * 2;
                        if (iter == 1) {
                            step1 = step * 4;
                        }
                    }
                    totalD = newTotalD;
                    System.out.println("totalD " + newTotalD);
                    System.out.println(iter + ">>>>step " + step);
                }//end iter

                for (int i = 0; i < nBasis; i++) {
                    IMolecule iMol = sim.box.getMoleculeList().get(i);
                    p.E(pos.position(iMol));
                    for (int j = 0; j < 3; j++) {
                        if (x0[i * 3 + j] == 0) {
                            x0[i * 3 + j] = p.getX(j);
                        }
                        p.setX(j, p.getX(j) + step * d[i * 3 + j]);
                    }
                    translator.setDestination(p);
                    translator.actionPerformed(iMol);
                    if (orient0[i].isZero()) {
                        orient0[i].E(iMol.getChildList().get(0).getPosition());
                        orient0[i].ME(pos.position(iMol));
                        orient0[i].normalize();
                    }

                    rTensor.setRotationAxis(axes[i], step * d[ng + i] / radianFac);
                    doTransform(iMol, sim.box.getBoundary(), rTensor);
                }
            }//iter
        }//outer

        double[] xf = new double[ng];
        double disp = 0.0;
        double angleDisp = 0.0;
        for (int i = 0; i < nBasis; i++) {
            IMolecule iMol = sim.box.getMoleculeList().get(i);
            p.E(pos.position(iMol));
            for (int j = 0; j < 3; j++) {
                xf[i * 3 + j] = p.getX(j);
                double dx = xf[i * 3 + j] - x0[i * 3 + j];
                disp += dx * dx;
            }

            orientf[i].E(iMol.getChildList().get(0).getPosition());
            orientf[i].ME(pos.position(iMol));
            orientf[i].normalize();
            angleDisp += orientf[i].Mv1Squared(orient0[i]);
        }
        disp = Math.sqrt(disp) / 3.0 / nBasis;
        angleDisp = Math.sqrt(angleDisp);
        System.out.println("\ndisp " + disp + "  angleDisp " + angleDisp);
        double newLatticeEnergy = Joule.UNIT.fromSim(meterPotentialEnergy.getDataAsScalar());
        newLatticeEnergy *= (1E-3 * Constants.AVOGADRO / nBasis);

        System.out.println("Old Lattice Energy (per molecule): " + latticeEnergy);
        System.out.println("New Lattice Energy (per molecule): " + newLatticeEnergy);

        Vector fSum = sim.space.makeVector();
        Vector r = sim.space.makeVector();
        Vector T = sim.space.makeVector();

        IteratorDirective id = new IteratorDirective();
        id.includeLrc = false;
        sim.potentialMaster.calculate(sim.box, id, pcForce);

        for (int I = 0; I < sim.box.getMoleculeList().size(); I++) {
            IMolecule iMol = sim.box.getMoleculeList().get(I);
            fSum.E(0);
            T.E(0);
            for (int j = 0; j < iMol.getChildList().size(); j++) {
                IAtom jAtom = iMol.getChildList().get(j);
                Vector fj = forces.get(jAtom);
                fSum.PE(fj);
                r.Ev1Mv2(jAtom.getPosition(), iMol.getChildList().get(2).getPosition());
                sim.box.getBoundary().nearestImage(r);
                r.XE(fj);
                T.PE(r);
            }
            System.out.println("fSum(" + I + ")" + fSum);
            System.out.println("Torque(" + I + ")" + T);
        }
        WriteConfiguration writeConfig = new WriteConfiguration(sim.space);
        writeConfig.setBox(sim.box);
        writeConfig.setDoApplyPBC(false);//false ... ok

        writeConfig.setFileName("finalPos.pos");
        writeConfig.actionPerformed();

        if (isGraphics) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "string");
            final DisplayBox display = new DisplayBox(sim, sim.box);
            simGraphic.add(display);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getTypeByName("H"), Color.GREEN);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getTypeByName("O"), Color.RED);
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("O"), 2.0);
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("H"), 1.0);

            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).setBackgroundColor(Color.white);
            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).setBoundaryFrameColor(Color.blue);
            simGraphic.makeAndDisplayFrame();
            if (!true) {
                return;
            }
        }

    }//end main

    protected static void doTransform(IMolecule molecule, Boundary boundary, Tensor rotationTensor) {
        IAtomList childList = molecule.getChildList();
        Vector O = childList.get(2).getPosition();//O
        O.PE(boundary.centralImage(O));// to wrap all O inside the BOX; next steps will move Hs and M with O to keep the conformation.
        for (int iChild = 0; iChild < childList.size(); iChild++) {
            if (iChild == 2) {//Ignore O
                continue;
            }
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(O);
            boundary.nearestImage(r);
            rotationTensor.transform(r);
            r.PE(O);
        }
    }

    // ********************* charges & epsilon_r ********************************//
    public static class ChargeAgentSourceRPM implements AtomLeafAgentManager.AgentSource<MyCharge> {
        private final Map<AtomType, MyCharge> myCharge;

        public ChargeAgentSourceRPM(SpeciesGeneral species, boolean isIce) {
            myCharge = new HashMap<>();
            double chargeH;
            if (isIce) {
                chargeH = Electron.UNIT.toSim(0.5897); //TIP4P/Ice
            } else {
                chargeH = Electron.UNIT.toSim(0.52); // TIP4P
            }
            myCharge.put(species.getTypeByName("H"), new MyCharge(chargeH));
            myCharge.put(species.getTypeByName("O"), new MyCharge(0));
            myCharge.put(species.getTypeByName("M"), new MyCharge(-2.0 * chargeH));
        }

        // *********************** set half(even # of particles ) as +ion, the other half -ion ***********************
        public MyCharge makeAgent(IAtom a, Box agentBox) {
            return myCharge.get(a.getType());
        }

        public void releaseAgent(MyCharge agent, IAtom atom, Box agentBox) {
            // Do nothing
        }
    }

    public static class SimParams extends ParameterBase {
        public String configFile = "config_sI";
        public int nBasis = 46;//sI
        public double[] a0 = new double[]{12.03, 12.03, 12.03};//sI
        public double rCutLJ = 1.0;
        public double rCutRealES = 1.724496;
        public double kCut = 0.1103499;
        public int nOuter = 0;
        public int nInner = 3;
        public boolean isIce = false, isGraphics = !false;
        public boolean includeM = false;
        int nX = 1;
        public int[] nC = new int[]{nX, nX, nX};
    }
}
