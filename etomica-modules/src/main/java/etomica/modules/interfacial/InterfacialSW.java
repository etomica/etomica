/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.interfacial;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManagerHard;
import etomica.potential.P2HardGeneric;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

import java.util.ArrayList;
import java.util.List;

/**
 * Simulation for interfacial tension module.  Simulation itself is just a
 * simple LJ system.
 *
 * @author Andrew Schultz
 */
public class InterfacialSW extends Simulation {

    public final SpeciesGeneral species, speciesGhost;
    public final SpeciesGeneral surfactant;
    public final Box box;
    public final IntegratorHard integrator;

    public final AtomType leafType, headType, tailType, ghostType;
    public final P2HardGeneric p2Head, p2HeadHead;
    public final P2HardGeneric p2TailTail, p2Tail, p2HeadTail;
    public final P2HardGeneric p2Bond;
    public final P2HardGeneric p2Ghost, p2GhostHead, p2GhostTail;
    public final ConfigurationPerturbed config;

    public InterfacialSW(Space _space) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        surfactant = new SpeciesBuilder(space)
                .addCount(AtomType.simpleFromSim(this), 1)
                .addCount(AtomType.simpleFromSim(this), 1)
                .withConformation(new ConformationLinear(space, 0.9))
                .setDynamic(true)
                .build();
        addSpecies(surfactant);

        speciesGhost = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesGhost);

        double pRange = 2.0;
        box = this.makeBox();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        NeighborListManagerHard neighborManager = new NeighborListManagerHard(getSpeciesManager(), box, 2, pRange, pmBonding.getBondingInfo());
//        NeighborManagerSimpleHard neighborManager = new NeighborManagerSimpleHard(box);
        PotentialComputePair potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        //species and potentials
        leafType = species.getLeafType();
        headType = surfactant.getAtomType(0); // head likes the monatomic species
        tailType = surfactant.getAtomType(1);
        ghostType = speciesGhost.getAtomType(0);

        //instantiate several potentials for selection in combo-box
        P2HardGeneric p2SW = new P2HardGeneric(new double[]{1.0, 1.5}, new double[]{Double.POSITIVE_INFINITY, -1.0}, false);
        potentialMaster.setPairPotential(leafType, leafType, p2SW);
        p2Head = new P2HardGeneric(new double[]{1.0, 1.5}, new double[]{Double.POSITIVE_INFINITY, -1.0}, false);
        potentialMaster.setPairPotential(leafType, headType, p2Head);
        p2HeadHead = new P2HardGeneric(new double[]{1.0, 1.5}, new double[]{Double.POSITIVE_INFINITY, -1.0}, false);
        potentialMaster.setPairPotential(headType, headType, p2HeadHead);

        p2TailTail = new P2HardGeneric(new double[]{1.0}, new double[]{Double.POSITIVE_INFINITY}, true);
        potentialMaster.setPairPotential(tailType, tailType, p2TailTail);
        p2Tail = new P2HardGeneric(new double[]{1.0}, new double[]{Double.POSITIVE_INFINITY}, true);
        potentialMaster.setPairPotential(leafType, tailType, p2Tail);
        p2HeadTail = new P2HardGeneric(new double[]{1.0}, new double[]{Double.POSITIVE_INFINITY}, true);
        potentialMaster.setPairPotential(headType, tailType, p2HeadTail);

        p2Ghost = new P2HardGhost(true);
        potentialMaster.setPairPotential(ghostType, leafType, p2Ghost);
        p2GhostHead = new P2HardGhost(true);
        potentialMaster.setPairPotential(ghostType, headType, p2GhostHead);
        p2GhostTail = new P2HardGhost(false);
        potentialMaster.setPairPotential(ghostType, tailType, p2GhostTail);

        p2Bond = new P2HardGeneric(new double[]{0.6, 1.0, 2}, new double[]{Double.POSITIVE_INFINITY, 0, Double.POSITIVE_INFINITY}, true) {
            @Override
            protected double collisionTime(Vector r12, Vector v12, int collisionState, double[] cd2) {
                double bij = r12.dot(v12);

                if (bij > 0.0) {
                    // moving apart
                    double r2 = r12.squared();
                    if (fixOverlap && r2 > cd2[1]) {
                        // overlapped collide now
                        double v2 = v12.squared();
                        return 0.001 * Math.sqrt(r2 / v2);
                    }
                }
                return super.collisionTime(r12, v12, collisionState, cd2);
            }

            public int getState(IAtom atom1, IAtom atom2, Vector r12) {
                if (fixOverlap) return 1;
                return super.getState(atom1, atom2, r12);
            }

            public double u(double r2) {
                if (fixOverlap) return 0;
                return super.u(r2);
            }
        };

        List<int[]> bonds = new ArrayList<>();
        bonds.add(new int[]{0, 1});
        pmBonding.setBondingPotentialPair(surfactant, p2Bond, bonds);

        int N = 643;  //number of atoms

        //controller and integrator
        integrator = new IntegratorHard(potentialMaster.getPairPotentials(), neighborManager, random,
                0.01, space.D() == 2 ? 0.4 : 1, box, getSpeciesManager(), pmBonding.getBondingInfo()) {
            public void doThermostat() {
                thermostatting = true;
                thermostatCount = thermostatInterval;

                if (initialized) {
                    IAtomList atomList = box.getLeafList();
                    int atomCount = atomList.size();
                    if (atomCount > 0) {
                        int index = random.nextInt(atomList.size());
                        IAtomKinetic a = (IAtomKinetic) atomList.get(index);
                        double m = a.getType().getMass();
                        if (m == Double.POSITIVE_INFINITY) return;
                        currentKineticEnergy -= 0.5 * m * a.getVelocity().squared();
                        randomizeMomentum(a);
                        if (a.getType() == ghostType) {
                            a.getVelocity().setX(0, 0);
                        }
                        currentKineticEnergy += 0.5 * m * a.getVelocity().squared();
                    }
                }
                thermostatting = false;
            }
        };
        if (space.D() == 2) {
            N = 300;
        }
        integrator.setIsothermal(true);
        integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        getController().addActivity(new ActivityIntegrate(integrator));

        //construct box
        Vector dim = space.makeVector();
        if (space.D() == 2) {
            dim.E(new double[]{30, 15});
        } else {
            dim.E(new double[]{12, 10, 10});
        }
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        ConfigurationLattice configLattice;
        if (space.D() == 2) {
            configLattice = new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space);
        } else {
            configLattice = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        }
        config = new ConfigurationPerturbed(configLattice, random, 0.01);
        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        Space space = Space2D.getInstance();
        if (args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        }

        InterfacialSW sim = new InterfacialSW(space);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }//end of main

    private static class P2HardGhost extends P2HardGeneric {
        public P2HardGhost(boolean sqw) {
            super(sqw ? new double[]{1, 1.5} : new double[]{1}, sqw ? new double[]{0, 0} : new double[]{0});
        }

        @Override
        public double collisionTime(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector v12, int collisionState, double falseTime) {
            return Double.POSITIVE_INFINITY;
        }
    }
}
