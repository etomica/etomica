/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.storage.VectorStorage;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListenerMD;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;

public class FixedWall implements IntegratorListenerMD {

    private final VectorStorage forces;
    protected final Box box;
    protected final ISpecies species;
    
    public FixedWall(Space space, Box box, VectorStorage forces, ISpecies species) {
        this.box = box;
        this.forces = forces;
        this.species = species;
    }
    
    public void integratorInitialized(IntegratorEvent e) {
        IMoleculeList molecules = box.getMoleculeList(species);
        double zTotMomentum = 0;
        double totMass = 0;
        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            for (int j = 0; j<atoms.size(); j++) {
                IAtomKinetic jAtom = (IAtomKinetic)atoms.get(j);
                double m = jAtom.getType().getMass();
                zTotMomentum += m*jAtom.getVelocity().getX(2);
                totMass += m;
            }
        }

        double vz = zTotMomentum/totMass;
        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            for (int j = 0; j<atoms.size(); j++) {
                IAtomKinetic jAtom = (IAtomKinetic)atoms.get(j);
                Vector v = jAtom.getVelocity();
                v.E(0);
                v.setX(2, vz);
            }
        }
    }

    public void integratorStepStarted(IntegratorEvent e) {}

    public void integratorStepFinished(IntegratorEvent e) {}
    
    public void integratorForcePrecomputed(IntegratorEvent e) {}
    
    public void integratorForceComputed(IntegratorEvent e) {
        IMoleculeList molecules = box.getMoleculeList(species);
        double fz = 0;
        double totMass = 0;
        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            for (int j = 0; j<atoms.size(); j++) {
                IAtom jAtom = atoms.get(j);
                fz += forces.get(jAtom).getX(2);
                totMass += jAtom.getType().getMass();
            }
        }

        
        fz /= totMass;
        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            for (int j = 0; j<atoms.size(); j++) {
                IAtom jAtom = atoms.get(j);
                Vector jf = forces.get(jAtom);
                jf.E(0);
                jf.setX(2, fz*jAtom.getType().getMass());
            }
        }
    }

}
