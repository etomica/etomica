package etomica.interfacial;

import etomica.api.*;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.space.Vector;
import etomica.space.Space;

public class FixedWall implements IIntegratorListenerMD {

    protected final AtomLeafAgentManager<MyAgent> agentManager;
    protected final Box box;
    protected final ISpecies species;
    
    public FixedWall(Space space, Box box, AtomLeafAgentManager<MyAgent> agentManager, ISpecies species) {
        this.box = box;
        this.agentManager = agentManager;
        this.species = species;
    }
    
    public void integratorInitialized(IIntegratorEvent e) {
        IMoleculeList molecules = box.getMoleculeList(species);
        double zTotMomentum = 0;
        double totMass = 0;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList atoms = molecules.getMolecule(i).getChildList();
            for (int j=0; j<atoms.getAtomCount(); j++) {
                IAtomKinetic jAtom = (IAtomKinetic)atoms.getAtom(j);
                double m = jAtom.getType().getMass();
                zTotMomentum += m*jAtom.getVelocity().getX(2);
                totMass += m;
            }
        }

        double vz = zTotMomentum/totMass;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList atoms = molecules.getMolecule(i).getChildList();
            for (int j=0; j<atoms.getAtomCount(); j++) {
                IAtomKinetic jAtom = (IAtomKinetic)atoms.getAtom(j);
                Vector v = jAtom.getVelocity();
                v.E(0);
                v.setX(2, vz);
            }
        }
    }

    public void integratorStepStarted(IIntegratorEvent e) {}

    public void integratorStepFinished(IIntegratorEvent e) {}
    
    public void integratorForcePrecomputed(IIntegratorEvent e) {}
    
    public void integratorForceComputed(IIntegratorEvent e) {
        IMoleculeList molecules = box.getMoleculeList(species);
        double fz = 0;
        double totMass = 0;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList atoms = molecules.getMolecule(i).getChildList();
            for (int j=0; j<atoms.getAtomCount(); j++) {
                IAtom jAtom = atoms.getAtom(j);
                fz += agentManager.getAgent(jAtom).force.getX(2);
                totMass += jAtom.getType().getMass();
            }
        }

        
        fz /= totMass;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList atoms = molecules.getMolecule(i).getChildList();
            for (int j=0; j<atoms.getAtomCount(); j++) {
                IAtom jAtom = atoms.getAtom(j);
                Vector jf = agentManager.getAgent(jAtom).force;
                jf.E(0);
                jf.setX(2, fz*jAtom.getType().getMass());
            }
        }
    }

}
