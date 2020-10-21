/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.action.BoxInflate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimLattice extends Simulation {

    public SpeciesGeneral species;
    public Box box;

    public SimLattice(int numAtoms, double temperature, double density, double w, int offsetDim) {
        super(Space3D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, numAtoms);
        Vector l = space.makeVector();
        l.E(10);
        for (int i = 0; i <= offsetDim; i++) {
            l.setX(i, 20);
        }
        box.getBoundary().setBoxSize(l);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
    }

    public static final double pi32 = Math.pow(Math.PI, 1.5);

    public static double q(double alpha, double betaw, double R2) {
        return Math.exp(-alpha*betaw*R2/(alpha+betaw))*pi32/Math.pow(alpha+betaw,1.5);
    }

    public static double uIntOverQ(double alpha, double betaw, double R2) {
        double abw = alpha+betaw;
        return (3*abw+2*alpha*alpha*R2)/(2*abw*abw);
    }
    
    public static void main(String[] args) {

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = false;
            params.numAtoms = 1000;
            params.density = 1.3;
            params.T = 3.39406146937929;
            params.w = Math.exp(7);
            params.alpha = 503;
        }

        final int numAtoms = params.numAtoms;
        final double temperature = params.T;
        final double density = params.density;
        boolean graphics = params.graphics;
        double w = params.w;
        int offsetDim = params.offsetDim;
        double alpha = params.alpha;
        boolean doq2 = params.doq2;

        if (!graphics) {
            System.out.println("Running lattice MC with N="+numAtoms+" at rho="+density+" T="+temperature);
            System.out.println("w: "+w);
            System.out.println("alpha: "+alpha);
        }

        final SimLattice sim = new SimLattice(numAtoms, temperature, density, w, offsetDim);
        Boundary boundary = sim.box.getBoundary();
        Vector offset = sim.space.makeVector();
        offset.setX(offsetDim, sim.box.getBoundary().getBoxSize().getX(offsetDim)*0.5);
        IAtomList atoms = sim.box.getLeafList();
        IAtom atom0 = atoms.get(0);
        IAtom atom1 = atoms.get(numAtoms/2);
        Vector p0 = atom0.getPosition();
        Vector p1 = atom1.getPosition();

        Vector dr = sim.space.makeVector();
        double betaw = w/temperature;

        double uSum = 0;
        double wSum = 0;
        for (int i=1; i<numAtoms; i++) {
            dr.Ev1Mv2(atoms.get(i).getPosition(),p0);
            dr.ME(offset);
            boundary.nearestImage(dr);
            double R2 = dr.squared();
            double q = q(alpha, betaw, R2);
            double u = uIntOverQ(alpha, betaw, R2);
//            System.out.println(Math.sqrt(R2)+" "+u+" "+q);
            if (doq2) {
                uSum += 2 * u * q * q;
                wSum += q * q;
            }
            else {
                uSum += u * q;
                wSum += q;
            }
        }
        double avg = uSum/wSum;
        if (doq2) avg /= 2;

        System.out.println("spring energy: "+w*avg);
    }

    public static class LjMC3DParams extends ParameterBase {
        public int numAtoms = 500;
        public double T = 2.0;
        public double density = 0.3;
        public boolean graphics = false;
        public double w = 1;
        public int offsetDim = 0;
        public double alpha = Double.POSITIVE_INFINITY;
        public boolean doq2 = false;
    }

}
