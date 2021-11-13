/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
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
import etomica.species.SpeciesGeneral;

/**
 * Applet to illustrate HTTP method
 * 
 * @author Tai Boon Tan
 */
public class HTTPSoftSphereSim extends Simulation {

    public IntegratorMC integrator;
    protected PotentialMasterList potentialMaster;
    protected double latticeEnergy;
    protected SpeciesGeneral species;
    protected P1ConstraintNbr p1Constraint;
    protected CoordinateDefinitionLeaf coordinateDefinition;
    protected double rc;

    public HTTPSoftSphereSim(Space _space, int numAtoms, double temperature) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        double density = 1.1964;
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        Boundary boundary = new BoundaryRectangularPeriodic(space, n * L);
        Box box = new Box(boundary, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        double nbrDistance = L / Math.sqrt(2);
        Primitive primitive = new PrimitiveCubic(space, n * L);
        rc = 0.495 * n * L;

        potentialMaster = new PotentialMasterList(getSpeciesManager(), box, 7, rc, BondingInfo.noBonding());
        potentialMaster.setDoDownNbrs(true);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        integrator.getEventManager().removeListener(potentialMaster);
        MCMoveAtomCoupled atomMove = new MCMoveAtomCoupled(potentialMaster, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        integrator.getMoveManager().addMCMove(atomMove);

        int[] nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, nCells);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, 12);
        if (potentialMaster instanceof PotentialMasterList) {
            potential = new P2SoftSphericalTruncated(space, potential, rc);

        } else {
            potential = new P2SoftSphericalTruncatedShifted(space, potential, rc);
        }
        AtomType sphereType = species.getLeafType();
        potentialMaster.setPairPotential(sphereType, sphereType, potential);

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

        p1Constraint = new P1ConstraintNbr(space, nbrDistance);
        p1Constraint.initBox(box);
        p1Constraint.setBox(box);
        atomMove.setConstraint(p1Constraint);

        potentialMaster.init();
        latticeEnergy = potentialMaster.computeAll(false);

        getController().addActivity(new ActivityIntegrate(integrator));
    }
    
}
