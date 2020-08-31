/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import Jama.Matrix;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.integrator.IntegratorMD;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;

/**
 * Simulation class of hard spheres in 1D or 3D that calculates the dq/dx
 * Jacobian.
 */
public class SimCalcJ extends Simulation {

    public SimCalcJ(Space _space, int numAtoms) {
        super(_space);

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        Basis basis;
        int[] nCells;
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1);
            bdry = new BoundaryRectangularPeriodic(space, numAtoms);
            basis = new BasisMonatomic(space);
            nCells = new int[]{numAtoms};
        } else if (space.D() == 2) {
            primitive = new PrimitiveCubic(space, 1);
            int n = (int) Math.round(Math.pow(numAtoms / 2, 1.0 / 2.0));
            nCells = new int[]{n, n};
            bdry = new BoundaryDeformableLattice(primitive, nCells);
            basis = new BasisOrthorhombicHexagonal();
        } else {
            primitive = new PrimitiveCubic(space, 1);
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            nCells = new int[]{n, n, n};
            bdry = new BoundaryDeformableLattice(primitive, nCells);
            basis = new BasisCubicFcc();
        }
        box = this.makeBox(bdry);
        box.setNMolecules(species, numAtoms);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 108;

        // parse arguments
        if (args.length > 0) {
            nA = Integer.parseInt(args[0]);
        }

        System.out.println("Calculating "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere Jacobian for N="+nA);

        // construct simulation
        SimCalcJ sim = new SimCalcJ(Space.getInstance(D), nA);

        // set up normal-mode meter
        WaveVectorFactory waveVectorFactory;
        if (D == 1) {
            waveVectorFactory = new WaveVectorFactory1D();
        } else if (D == 2) {
            waveVectorFactory = new WaveVectorFactory2D(sim.primitive, sim.space);
        } else {
            waveVectorFactory = new WaveVectorFactorySimple(sim.primitive, sim.space);
        }
        CalcJacobian meterJacobian = new CalcJacobian();
        meterJacobian.setWaveVectorFactory(waveVectorFactory);
        meterJacobian.setCoordinateDefinition(sim.coordinateDefinition);

        double[][] jac = meterJacobian.getJacobian();
        
        if (jac.length < 100) {
            for (int i=0; i<jac.length; i++) {
                for (int j=0; j<jac[0].length; j++) {
                    System.out.print(jac[i][j]+"  ");
                }
                System.out.println();
            }
        }
        
        Matrix m = new Matrix(jac);
        double d = Math.abs(m.det());
        System.out.println("dx/dq det = "+1/d+" (log2 det = "+Math.log(1.0/d)/Math.log(2)+")");
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMD integrator;

    public Box box;
    public Boundary bdry;
    public Primitive primitive;
    public CoordinateDefinition coordinateDefinition;
}
