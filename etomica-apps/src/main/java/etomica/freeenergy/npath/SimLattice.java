/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.action.BoxInflate;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.atom.AtomSetSinglet;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimLattice extends Simulation {

    public SpeciesSpheresMono species;
    public IBox box;
    public P1ImageHarmonic p1ImageHarmonic;

    public SimLattice(int numAtoms, double temperature, double density, double w, int offsetDim) {
        super(Space3D.getInstance());
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        IVectorMutable l = space.makeVector();
        l.E(10);
        for (int i=0; i<=offsetDim; i++) {
            l.setX(i,20);
        }
        box.getBoundary().setBoxSize(l);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        IVectorMutable offset = space.makeVector();
        offset.setX(offsetDim, box.getBoundary().getBoxSize().getX(offsetDim)*0.5);
        p1ImageHarmonic = new P1ImageHarmonic(space, offset, w);

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
    }
    
    public static void main(String[] args) {

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = false;
            params.numAtoms = 1000;
            params.density = 1.3;
            params.T = 2.00093081574;
            params.w = 0.1353352832366127;
        }

        final int numAtoms = params.numAtoms;
        final double temperature = params.T;
        final double density = params.density;
        boolean graphics = params.graphics;
        double w = params.w;
        int offsetDim = params.offsetDim;

        if (!graphics) {
            System.out.println("Running lattice MC with N="+numAtoms+" at rho="+density+" T="+temperature);
            System.out.println("w: "+w);
        }

        final SimLattice sim = new SimLattice(numAtoms, temperature, density, w, offsetDim);
        sim.p1ImageHarmonic.setBox(sim.box);
        IAtomList atoms = sim.box.getLeafList();
        IAtom atom = atoms.getAtom(numAtoms/2);
        IVectorMutable p = atom.getPosition();
        AtomSetSinglet singlet = new AtomSetSinglet(atom);
        double u = 2*sim.p1ImageHarmonic.energy(singlet);
        double exp = Math.exp(-u/temperature);
        double uSum = u*exp;
        double wSum = exp;
        for (int i=1; i<numAtoms; i++) {
            if (i==numAtoms/2) continue;
            p.E(atoms.getAtom(i).getPosition());
            u = 2*sim.p1ImageHarmonic.energy(singlet);
            exp = Math.exp(-u/temperature);
            uSum += u*exp;
            wSum += exp;
        }
        double avg = uSum/wSum;

        System.out.println("spring energy: "+avg);
    }

    public static class LjMC3DParams extends ParameterBase {
        public int numAtoms = 500;
        public double T = 2.0;
        public double density = 0.3;
        public boolean graphics = false;
        public double w = 1;
        public int offsetDim = 0;
    }

}
