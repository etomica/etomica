package etomica.modules.droplet;

import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.config.Configuration;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class ConfigurationDroplet implements Configuration {

    public ConfigurationDroplet(IRandom random, ISpace space) {
        this.random = random;
        axis = space.makeVector();
        center = space.makeVector();
    }
    
    public void initializeCoordinates(IBox box) {
        
        double factor = (1+deformation) / (1-deformation);
        axis.E(Math.pow(factor, -1.0/3.0));
        axis.setX(2,1.0/(axis.getX(0)*axis.getX(0)));
        
        IAtomList leafList = box.getLeafList();
        int numAtoms = leafList.getAtomCount();
        for (int i=0; i<numAtoms; i++) {
            IVectorMutable r = ((IAtomPositioned)leafList.getAtom(i)).getPosition();
            ((IVectorRandom)r).setRandomInSphere(random);
            r.TE(axis);
        }

        center.E(0);
        for (int i=0; i<leafList.getAtomCount(); i++) {
            center.PE(((IAtomPositioned)leafList.getAtom(i)).getPosition());
        }
        center.TE(1.0/leafList.getAtomCount());

        for (int i=0; i<leafList.getAtomCount(); i++) {
            ((IAtomPositioned)leafList.getAtom(i)).getPosition().ME(center);
        }
    }
    
    public void setDeformation(double newDeformation) {
        deformation = newDeformation;
    }
    
    public double getDeformation() {
        return deformation;
    }
    
    protected double deformation;
    protected final IRandom random;
    protected final IVectorMutable axis, center;
}
