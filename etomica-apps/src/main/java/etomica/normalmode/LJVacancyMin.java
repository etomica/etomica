/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;


public class LJVacancyMin extends Simulation {

    public final CoordinateDefinitionLeaf coordinateDefinition;
    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public PotentialMasterList potentialMaster;
    public Potential2SoftSpherical potential;
    public SpeciesSpheresMono species;
    public LJVacancyMin(Space _space, int numAtoms, double density, double rc, boolean ss) {
        super(_space);

        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, space);

        // TARGET
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);

        primitive = new PrimitiveCubic(space, n * L);

        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        potential = ss ? new P2SoftSphere(space, 1.0, 4.0, 12) : new P2LennardJones(space, 1.0, 1.0);
        potential = new P2SoftSphericalTruncated(space, potential, rc);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        potentialMaster.lrcMaster().setEnabled(false);

        int cellRange = 7;
        potentialMaster.setRange(rc);
        potentialMaster.setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMaster.getNeighborManager(box).reset();
        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }

        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimLJHTTI.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        if (args.length == 0) {
            params.numAtoms = 500;
            params.density = 1;
            params.rc = 3;
            params.ss = false;
            params.verbose = true;
        }
        else {
            ParseArgs.doParseArgs(params, args);
        }
        boolean ss = params.ss;
        double density = params.density;
        final int numAtoms = params.numAtoms;
        double rc = params.rc;
        boolean verbose = params.verbose;

        System.out.println("Minimizing "+(ss?"soft-sphere":"Lennard-Jones")+" with a vacancy");
        System.out.println(numAtoms+" atoms at density "+density);
        System.out.println("rc: "+rc);

        //instantiate simulation
        final LJVacancyMin sim = new LJVacancyMin(Space.getInstance(3), numAtoms, density, rc*Math.pow(density, -1.0/3.0), ss);
        //start simulation
        IteratorDirective all = new IteratorDirective();
        PotentialCalculationEnergySum pcEnergy = new PotentialCalculationEnergySum();
        MeterPressure meterP = new MeterPressure(sim.space);
        meterP.setBox(sim.box);
        meterP.setIncludeLrc(false);
        meterP.setPotentialMaster(sim.potentialMaster);
        meterP.setTemperature(0);
        double pLat = meterP.getDataAsScalar();
        System.out.println("pLat: "+pLat);
        pcEnergy.zeroSum();
        sim.potentialMaster.calculate(sim.box, all, pcEnergy);
        double uLat = pcEnergy.getSum();
        System.out.println("uLat: "+uLat/numAtoms);

        sim.box.removeMolecule(sim.box.getMoleculeList().get(0));

        pcEnergy.zeroSum();
        sim.potentialMaster.calculate(sim.box, all, pcEnergy);
        double u0 = pcEnergy.getSum();

        PotentialCalculationForceSum pc = new PotentialCalculationForceSum();
        AtomLeafAgentManager<Vector> forceManager = new AtomLeafAgentManager<>(a -> sim.space.makeVector(), sim.box);
        pc.setAgentManager(forceManager);

        double[] d = new double[3*(numAtoms-1)];
        double[] dir = new double[d.length];
        double step0 = 1e-5;
        double uLast = Double.POSITIVE_INFINITY;
        for (int i=0; i<200; i++) {
            pcEnergy.zeroSum();
            sim.potentialMaster.calculate(sim.box, all, pcEnergy);
            double u = pcEnergy.getSum();
            if (u > uLast) break;
            uLast = u;
            // compute direction of minimization
            pc.reset();
            sim.potentialMaster.calculate(sim.box, all, pc);
            double dSum = 0;
            for (int j=0; j<numAtoms-1; j++) {
                IAtom jAtom = sim.box.getLeafList().get(j);
                Vector fj = forceManager.getAgent(jAtom);
                Vector pj = jAtom.getPosition();
                for (int k=0; k<pj.getD(); k++) {
                    d[3*j+k] = fj.getX(k);
                    dSum += d[3*j+k]*d[3*j+k];
                }
            }
            double dNorm0 = Math.sqrt(dSum);
            if (dNorm0 < 1e-8) {
                if (verbose) System.out.println("==> "+dNorm0);
                break;
            }
            if (verbose) System.out.print(String.format("%4d  %25.15e  %10.4e", i, u-u0, dNorm0));
            for (int l=0; l<d.length; l++) {
                dir[l] = d[l]/dNorm0;
            }
            // take a small step
            for (int j=0; j<numAtoms-1; j++) {
                IAtom jAtom = sim.box.getLeafList().get(j);
                Vector pj = jAtom.getPosition();
                for (int k=0; k<pj.getD(); k++) {
                    pj.setX(k, pj.getX(k)+step0*dir[3*j+k]);
                }
            }

            pc.reset();
            sim.potentialMaster.calculate(sim.box, all, pc);
            dSum = 0;
            for (int j=0; j<numAtoms-1; j++) {
                IAtom jAtom = sim.box.getLeafList().get(j);
                Vector fj = forceManager.getAgent(jAtom);
                Vector pj = jAtom.getPosition();
                for (int k=0; k<pj.getD(); k++) {
                    d[3*j+k] = fj.getX(k);
                    dSum += d[3*j+k]*d[3*j+k];
                }
            }
            double dNorm1 = Math.sqrt(dSum);
            if (verbose) System.out.print(String.format(" =>  %10.4e", dNorm1));

            // take a second step; try to jump to 0
            // dNorm2 = dNorm1 + (dNorm1-dNorm0)*step1/step0 = 0
            double step1 = -step0*dNorm1/(dNorm1-dNorm0);
            if (verbose) System.out.print(String.format("  %10.4e\n", step1));
            for (int j=0; j<numAtoms-1; j++) {
                IAtom jAtom = sim.box.getLeafList().get(j);
                Vector pj = jAtom.getPosition();
                for (int k=0; k<pj.getD(); k++) {
                    pj.setX(k, pj.getX(k)+step1*dir[3*j+k]);
                }
            }
            step0 = step1/100;
        }
        pcEnergy.zeroSum();
        sim.potentialMaster.calculate(sim.box, all, pcEnergy);
        double u1 = pcEnergy.getSum();
        System.out.println("deltaU: "+(u1-uLat));
        double pLat1 = meterP.getDataAsScalar();
        System.out.println("deltaP: "+(pLat1-pLat));
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numAtoms = 256;
        public double density = 1.28;
        public double rc = 2.5;
        public boolean ss = false;
        public boolean verbose = false;
    }
}
