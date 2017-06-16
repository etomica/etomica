/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
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
import etomica.species.SpeciesSpheresMono;

/**
 * Applet to illustrate HTTP method
 * 
 * @author Tai Boon Tan
 */
public class HTTPSoftSphereSim extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    protected PotentialMaster potentialMaster;
    protected double latticeEnergy;
    protected SpeciesSpheresMono species;
    protected P1ConstraintNbr p1Constraint;
    protected CoordinateDefinitionLeaf coordinateDefinition;
    public HTTPSoftSphereSim(Space _space, int numAtoms,  double temperature) {
        super(_space);

        double density = 1.1964;

        potentialMaster = new PotentialMasterList(this, space);

        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        Box box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        MCMoveAtomCoupled atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        double L = Math.pow(4.0/density, 1.0/3.0);
        double nbrDistance = L / Math.sqrt(2);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        Primitive primitive = new PrimitiveCubic(space, n*L);
        double rc = 0.495*n*L;

        int[] nCells = new int[]{n,n,n};
        Boundary boundary = new BoundaryRectangularPeriodic(space, n * L);
        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, nCells);

        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, 12);
     	if(potentialMaster instanceof PotentialMasterList){
			potential = new P2SoftSphericalTruncated(space, potential, rc);

        } else {
			potential = new P2SoftSphericalTruncatedShifted(space, potential, rc);

        }
        atomMove.setPotential(potential);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

        p1Constraint = new P1ConstraintNbr(space, nbrDistance, this);
        p1Constraint.initBox(box);
        atomMove.setConstraint(p1Constraint);

        potentialMaster.lrcMaster().setEnabled(false);

        integrator.setBox(box);

		if (potentialMaster instanceof PotentialMasterList) {
            int cellRange = 7;
            ((PotentialMasterList)potentialMaster).setRange(rc);
            ((PotentialMasterList)potentialMaster).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMaster).getNeighborManager(box).reset();
            int potentialCells = ((PotentialMasterList)potentialMaster).getNbrCellManager(box).getLattice().getSize()[0];
            if (potentialCells < cellRange*2+1) {
                throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
            }
		}

        latticeEnergy = meterPE.getDataAsScalar();

		ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        if (potentialMaster instanceof PotentialMasterList) {
            // extend potential range, so that atoms that move outside the truncation range will still interact
            // atoms that move in will not interact since they won't be neighbors
            ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
        }
    }
    
}
