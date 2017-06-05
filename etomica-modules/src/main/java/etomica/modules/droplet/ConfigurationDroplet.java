/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.config.Configuration;
import etomica.space.Space;

public class ConfigurationDroplet implements Configuration {

    public ConfigurationDroplet(IRandom random, Space space) {
        this.random = random;
        axis = space.makeVector();
        center = space.makeVector();
    }
    
    public void initializeCoordinates(Box box) {
        
        double factor = (1+deformation) / (1-deformation);
        axis.E(Math.pow(factor, -1.0/3.0));
        axis.setX(2,1.0/(axis.getX(0)*axis.getX(0)));
        
        IAtomList leafList = box.getLeafList();
        int numAtoms = leafList.getAtomCount();
        for (int i=0; i<numAtoms; i++) {
            Vector r = leafList.getAtom(i).getPosition();
            r.setRandomInSphere(random);
            r.TE(axis);
        }

        center.E(0);
        for (int i=0; i<leafList.getAtomCount(); i++) {
            center.PE(leafList.getAtom(i).getPosition());
        }
        center.TE(1.0/leafList.getAtomCount());

        for (int i=0; i<leafList.getAtomCount(); i++) {
            leafList.getAtom(i).getPosition().ME(center);
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
    protected final Vector axis, center;
}
