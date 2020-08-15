/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.action.BoxInflate;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.Argon;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.MoleculeArrayList;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

/**
 * Atomic simulation for Droplet module.
 *
 * @author Andrew Schultz
 */
public class DropletAtomic extends Simulation {

    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono species;
    public final Box box;
    public final IntegratorVelocityVerlet integrator;

    public final PotentialMasterList potentialMaster;
    public final P2LennardJones p2LJ;
    public final P2SoftSphericalTruncatedShifted p2LJt;
    public final P1Smash p1Smash;
    protected int nNominalAtoms;
    protected double dropRadius;
    protected double xDropAxis;
    protected double density;
    protected double sigma;

    public DropletAtomic() {
        super(Space3D.getInstance());
        //species
        species = new SpeciesSpheresMono(space, Argon.INSTANCE);
        species.setIsDynamic(true);
        addSpecies(species);

        //construct box
        box = this.makeBox(new BoundaryRectangularPeriodic(space));
        double pRange = 3;
        sigma = 3.35;
        nNominalAtoms = 32000;
        dropRadius = 0.4;
        xDropAxis = 1;
        density = 0.6;

        potentialMaster = new PotentialMasterList(this, sigma * pRange * 1.5, space);

        //controller and integrator
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setTimeStep(0.005);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(5000);
        getController().addActivity(new ActivityIntegrate(integrator));
        integrator.setTemperature(Kelvin.UNIT.toSim(118));

        //potentials
        AtomType leafType = species.getLeafType();

        p2LJ = new P2LennardJones(space);
        p2LJ.setEpsilon(Kelvin.UNIT.toSim(118));
        p2LJ.setSigma(sigma);
        p2LJt = new P2SoftSphericalTruncatedForceShifted(space, p2LJ, sigma * pRange);
        potentialMaster.addPotential(p2LJt, new AtomType[]{leafType, leafType});

        p1Smash = new P1Smash(space);
        p1Smash.setG(4);
        potentialMaster.addPotential(p1Smash, new AtomType[]{leafType});


        makeDropShape();

        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        DropletAtomic sim = new DropletAtomic();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), Long.MAX_VALUE);
    }//end of main
    
    public void makeDropShape() {
        box.setNMolecules(species, nNominalAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density/(sigma*sigma*sigma));
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        IAtomList leafList = box.getLeafList();
        Vector v = space.makeVector();
        Vector dim = box.getBoundary().getBoxSize();
        double dropRadiusSq = 0.25*dropRadius*dropRadius*dim.getX(0)*dim.getX(0);
        int ambientCount = 0;
        MoleculeArrayList outerMolecules = new MoleculeArrayList();
        for (int i = 0; i<leafList.size(); i++) {
            v.E(leafList.get(i).getPosition());
            v.setX(0, v.getX(0)/xDropAxis);
            if (v.squared() > dropRadiusSq) {
                ambientCount++;
                if (ambientCount == 20) {
                    ambientCount = 0;
                }
                else {
                    outerMolecules.add(leafList.get(i).getParentGroup());
                }
            }
        }
        for (int i = 0; i<outerMolecules.size(); i++) {
            box.removeMolecule(outerMolecules.get(i));
        }
    }
}
