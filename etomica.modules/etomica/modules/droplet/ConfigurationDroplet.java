package etomica.modules.droplet;

import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.config.Configuration;
import etomica.space.IVectorRandom;

public class ConfigurationDroplet implements Configuration {

    public ConfigurationDroplet(IRandom random) {
        this.random = random;
    }
    
    public void initializeCoordinates(IBox box) {
        IAtomList leafList = box.getLeafList();
        int numAtoms = leafList.getAtomCount();
        for (int i=0; i<numAtoms; i++) {
            IVectorMutable r = ((IAtomPositioned)leafList.getAtom(i)).getPosition();
            ((IVectorRandom)r).setRandomInSphere(random);
            r.TE(2);
        }
    }

    protected final IRandom random;
}
