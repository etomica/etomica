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

    public PotentialMoleculePairImplicit(
            Space space, SpeciesManager sm,
            SiteReconstructor reconstructor1, SiteReconstructor reconstructor2) {
        super(space, sm);
        this.dr = space.makeVector();
        this.reconstructor1 = reconstructor1;
        this.reconstructor2 = reconstructor2;
    }
    public PotentialMoleculePairImplicit(
            Space space, IPotential2[][] atomPotentials,
            SiteReconstructor reconstructor1, SiteReconstructor reconstructor2) {
        super(space, atomPotentials);
        this.dr = space.makeVector();
        this.reconstructor1 = reconstructor1;
        this.reconstructor2 = reconstructor2;
    }

    public double energy(IMolecule molecule1, IMolecule molecule2) {
        reconstructor1.reconstructSites(molecule1);
        reconstructor2.reconstructSites(molecule2);

        SiteSet sites1 = reconstructor1.getSites();
        SiteSet sites2 = reconstructor2.getSites();

        double u = 0.0;
        for (int i = 0; i < sites1.nSites; i++) {
            IPotential2[] iPotentials = atomPotentials[sites1.types[i].getIndex()];
            for (int j = 0; j < sites2.nSites; j++) {
                IPotential2 p2 = iPotentials[sites2.types[j].getIndex()];
                if (p2 == null) continue;
                dr.Ev1Mv2(sites2.pos[j], sites1.pos[i]);
                u += p2.u(dr.squared());
            }
        }
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

        SiteSet getSites();

        void reconstructSites(IMolecule molecule);
    }
    public static class SiteReconstructorNull implements SiteReconstructor{
        public SiteSet sites;
        public PotentialMoleculePairImplicit.SiteSet getSites() {
            return sites;
        }

        public void reconstructSites(IMolecule molecule){
            if(sites == null) {
                AtomType[] types = new AtomType[molecule.getChildList().size()];
                for (IAtom a: molecule.getChildList()){
                    types[a.getIndex()] = a.getType();
                }
                sites = new SiteSet(types);

            }
            for (IAtom a: molecule.getChildList()){
                sites.pos[a.getIndex()].E(a.getPosition());
            }
        }
    }
}