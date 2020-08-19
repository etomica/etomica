/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.interfacial;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardBond;
import etomica.potential.P2HardSphere;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialGroup;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresHetero;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation for interfacial tension module.  Simulation itself is just a
 * simple LJ system.
 *
 * @author Andrew Schultz
 */
public class InterfacialSW extends Simulation {

    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono species;
    public final SpeciesSpheresHetero surfactant;
    public final Box box;
    public final IntegratorHard integrator;

    public final AtomType leafType, headType, tailType;
    public final P2SquareWell p2Head, p2HeadHead;
    public final P2HardSphere p2TailTail, p2Tail, p2HeadTail;
    public final P2HardBond p2Bond;

    public InterfacialSW(Space _space) {
        super(_space);

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        surfactant = new SpeciesSpheresHetero(this, space, 2);
        surfactant.setIsDynamic(true);
        surfactant.setChildCount(new int[]{1, 1});
        surfactant.setTotalChildren(2);
        ((ConformationLinear) surfactant.getConformation()).setBondLength(0.9);
        addSpecies(surfactant);

        double pRange = 2.0;
        PotentialMasterList potentialMaster = new PotentialMasterList(this, pRange, space);

        int N = 643;  //number of atoms

        //controller and integrator
        box = this.makeBox();
        integrator = new IntegratorHard(this, potentialMaster, box);
        if (space.D() == 2) {
            integrator.setTemperature(0.4);
            N = 300;
        }
        integrator.setIsothermal(true);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        getController().addActivity(new ActivityIntegrate(integrator));
        integrator.setTimeStep(0.01);

        //species and potentials
        leafType = species.getLeafType();
        headType = surfactant.getAtomType(0); // head likes the monatomic species
        tailType = surfactant.getAtomType(1);

        //instantiate several potentials for selection in combo-box
        P2SquareWell p2SW = new P2SquareWell(space, 1.0, 1.5, 1.0, true);
        potentialMaster.addPotential(p2SW, new AtomType[]{leafType, leafType});
        p2Head = new P2SquareWell(space, 1.0, 1.5, 1.0, true);
        potentialMaster.addPotential(p2Head, new AtomType[]{leafType, headType});
        p2HeadHead = new P2SquareWell(space, 1.0, 1.5, 1.0, true);
        potentialMaster.addPotential(p2HeadHead, new AtomType[]{headType, headType});

        p2TailTail = new P2HardSphere(space, 1.0, true);
        potentialMaster.addPotential(p2TailTail, new AtomType[]{tailType, tailType});
        p2Tail = new P2HardSphere(space, 1.0, true);
        potentialMaster.addPotential(p2Tail, new AtomType[]{leafType, tailType});
        p2HeadTail = new P2HardSphere(space, 1.0, true);
        potentialMaster.addPotential(p2HeadTail, new AtomType[]{headType, tailType});

        p2Bond = new P2HardBond(space, 0.8, 0.2, true);
        PotentialGroup p1Surfactant = potentialMaster.makePotentialGroup(1);
        p1Surfactant.addPotential(p2Bond, ApiBuilder.makeAdjacentPairIterator());
        potentialMaster.addPotential(p1Surfactant, new ISpecies[]{surfactant});

        //construct box
        Vector dim = space.makeVector();
        if (space.D() == 2) {
            dim.E(new double[]{30, 15});
        } else {
            dim.E(new double[]{12, 10, 10});
        }
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        if (space.D() == 2) {
            new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        } else {
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        }
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
    }
    
    public static void main(String[] args) {
        Space space = Space2D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
            
        InterfacialSW sim = new InterfacialSW(space);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }//end of main
}
