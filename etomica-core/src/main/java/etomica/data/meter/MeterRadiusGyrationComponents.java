/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataTensor;
import etomica.data.types.DataVector;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Length;

/**
 * Meter for tabulation of the radius of gyration of a set of chain molecules. 
 * 
 * @author David Kofke
 */
public class MeterRadiusGyrationComponents implements IDataSource {
    private final Box box;
    private final Vector cm, realPos;
    private final Vector dr;
    private final DataVector data;
    private final Tensor t;
    private final DataVector.DataInfoVector datainfo;
    private final DataTag tag;
    public MeterRadiusGyrationComponents(Box box) {

        this.box = box;
        cm = box.getSpace().makeVector();
        realPos = box.getSpace().makeVector();
        dr = box.getSpace().makeVector();
        data= new DataVector(box.getSpace());
        t=box.getSpace().makeTensor();
        datainfo=new DataVector.DataInfoVector("Radius Of Gyration",Length.DIMENSION,box.getSpace());
        tag=new DataTag();
        datainfo.addTag(tag);
    }



    @Override
    public IData getData() {
        if (box == null)
            throw new IllegalStateException(
                    "must call setBox before using meter");
        t.E(0);
        Boundary boundary = box.getBoundary();
        int nLeafAtomsTot = 0;
        double r2Tot = 0.0;
        for (IMolecule molecule : box.getMoleculeList()) {
            // loop over molecules
            IAtomList childList = molecule.getChildList();
            if (childList.size() < 2) {
                // a monatomic molecule
                continue;
            }

            // find center of mass
            //do the first iterate explicitly, assume there is at least
            // one leaf atom
            IAtom firstAtom = childList.get(0);
            int nLeafAtoms = 1;
            realPos.E(firstAtom.getPosition());
            cm.E(realPos);
            Vector prevPosition = firstAtom.getPosition();
            for (int iChild = 1; iChild < childList.size(); iChild++) {
                IAtom a = childList.get(iChild);
                nLeafAtoms++;
                Vector position = a.getPosition();
                dr.Ev1Mv2(position, prevPosition);
                //molecule might be wrapped around the box.  calculate
                //the real difference in position
                boundary.nearestImage(dr);
                //realPos is now the effective position of a
                realPos.PE(dr);
                cm.PE(realPos);
                prevPosition = position;
            }
            cm.TE(1.0 / nLeafAtoms);
            // calculate Rg^2 for this chain
            double r2 = 0.0;
            prevPosition = firstAtom.getPosition();
            realPos.E(firstAtom.getPosition());
            for (int iChild = 0; iChild < childList.size(); iChild++) {
                IAtom a = childList.get(iChild);
                Vector position = a.getPosition();
                dr.Ev1Mv2(position, prevPosition);
                //molecule might be wrapped around the box.  calculate
                //the real difference in position
                boundary.nearestImage(dr);
                //realPos is now the effective position of a
                realPos.PE(dr);
                dr.Ev1Mv2(realPos, cm);// = realPos.M(cm);
               t.PEv1v2(dr,dr);
                prevPosition = position;
            }
            r2Tot += r2;
            nLeafAtomsTot += nLeafAtoms;
        }
       t.TE( 1.0 /nLeafAtomsTot);
        Matrix m = new Matrix(t.toArray(),t.D());
        EigenvalueDecomposition eigen= new EigenvalueDecomposition(m);
       data.x.E(eigen.getRealEigenvalues());

        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return datainfo;
    }
}
