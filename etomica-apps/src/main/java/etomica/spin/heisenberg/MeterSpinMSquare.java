/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.data.IEtomicaDataSource;
import etomica.space.ISpace;
import etomica.space3d.Vector3D;
import etomica.units.Undefined;


/**
 * returns dipole square
 *
 * @author  Weisong Lin
 *
 */
public class MeterSpinMSquare extends DataSourceScalar implements IEtomicaDataSource {

    /**
     * 
     */
	
    public MeterSpinMSquare(ISpace space,IBox box,double dipoleMagnitude) {
        super("Spin",Undefined.DIMENSION);
        sum = space.makeVector();
        this.box = box;
        this.dipoleMagnitude = dipoleMagnitude;
    }

    public double getDataAsScalar() {
        sum.E(0.0);

//        int count = 0;
//        iterator.setBox(box);
//        iterator.reset();
//        for (IAtomOriented atom = (IAtomOriented) iterator.nextAtom(); atom != null;
//             atom = (IAtomOriented) iterator.nextAtom()) {
//            sum.PE(atom.getOrientation().getDirection());
//            count++;
//        }

        // test for new way to get MSquare
        if (box == null) throw new IllegalStateException("no box");
        IAtomList leafList = box.getLeafList();
        int nM = leafList.getAtomCount();
        for (int i = 0;i < nM; i++){
            IAtomOriented atom = (IAtomOriented) leafList.getAtom(i);
            sum.PE(atom.getOrientation().getDirection());
        }//i loop

        return sum.squared()*dipoleMagnitude*dipoleMagnitude;
    }
    
    
    

    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(IBox box) {
        this.box = box;
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private final IVectorMutable sum;
    private final double dipoleMagnitude;
}
