/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;import etomica.molecule.IMolecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;import etomica.species.SpeciesManager;

/**
 * Potential that includes virtual sites whose position can be determined
 * from concrete IAtom positions using SiteReconstructors.
 */
public class PotentialMoleculePairImplicit extends PotentialMoleculePair {


    protected final SiteReconstructor reconstructor1, reconstructor2;
    protected final Vector dr;
    public static boolean debug;
    public SiteSet sites1;
    public SiteSet sites2;

    public PotentialMoleculePairImplicit(
            Space space, SpeciesManager sm,
            SiteReconstructor reconstructor1, SiteReconstructor reconstructor2) {
        super(space, sm);
        this.dr = space.makeVector();
        this.reconstructor1 = reconstructor1;
        this.reconstructor2 = reconstructor2;
        sites1 = reconstructor1.makeSites();
        sites2 = reconstructor2.makeSites();
    }
    public PotentialMoleculePairImplicit(
            Space space, IPotential2[][] atomPotentials,
            SiteReconstructor reconstructor1, SiteReconstructor reconstructor2) {
        super(space, atomPotentials);
        this.dr = space.makeVector();
        this.reconstructor1 = reconstructor1;
        this.reconstructor2 = reconstructor2;
        sites1 = reconstructor1.makeSites();
        sites2 = reconstructor2.makeSites();
    }

    public double energy(IMolecule molecule1, IMolecule molecule2) {
        reconstructor1.reconstructSites(molecule1, sites1);
        reconstructor2.reconstructSites(molecule2, sites2);


        double u = 0.0;
        for (int i = 0; i < sites1.nSites; i++) {
            IPotential2[] iPotentials = atomPotentials[sites1.types[i].getIndex()];
            for (int j = 0; j < sites2.nSites; j++) {
                IPotential2 p2 = iPotentials[sites2.types[j].getIndex()];
                if (p2 == null) continue;
                dr.Ev1Mv2(sites2.pos[j], sites1.pos[i]);
                if (debug) System.out.println(i + " " + j + " " + dr.squared() + " " +p2.u(dr.squared()));
                u += p2.u(dr.squared());
            }
        }
        if (debug) System.exit(0);
        return u;
    }

    public static class SiteSet {
        public final Vector[] pos;
        public final AtomType[] types;
        public int nSites;

        public SiteSet(AtomType[] types) {
            nSites = types.length;
            pos = new Vector[nSites];
            this.types = types;
            for (int i = 0; i < nSites; i++) {
                pos[i] = new Vector3D();
            }
        }
    }

    public interface SiteReconstructor {

        SiteSet makeSites();

        void reconstructSites(IMolecule molecule, SiteSet sites);
    }
    public static class SiteReconstructorNull implements SiteReconstructor{
        public final ISpecies species;
        public SiteReconstructorNull(ISpecies species){
            this.species = species;
        }
        public PotentialMoleculePairImplicit.SiteSet makeSites() {
            AtomType[] types = new AtomType[species.getLeafAtomCount()];
            for (int i = 0; i < types.length; i++){
                types[i] = species.getAtomType(i);
            }
            return new SiteSet(types);
        }


        public void reconstructSites(IMolecule molecule, SiteSet sites){
            for (IAtom a: molecule.getChildList()){
                sites.pos[a.getIndex()].E(a.getPosition());
            }
        }
    }
}