package etomica.data.meter;


import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataVector;
import etomica.molecule.IMolecule;
import etomica.molecule.iterator.MoleculeIteratorAllMolecules;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Area;
import etomica.units.dimensions.Length;

/**
 * Meter for tabulation of the gyration tensor of a set of chain molecules.
 *
 * @author Arpit Bansal
 */

public class MeterGyrationTensor implements IDataSource {
    public MeterGyrationTensor(Space space){
        data = new DataDoubleArray(4);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("Gyration Tensor", Area.DIMENSION, new int[]{4});
        tag = new DataTag();
        iterator = new MoleculeIteratorAllMolecules();
        cm = space.makeVector();
        realPos = space.makeVector();
        dr = space.makeVector();
    }


    /**
     * Mutator method for the iterator that generates the atom pairs used to
     * tabulate the ROG. By setting this iterator the
     * meter can be configured to compute pair distribution for any set of atom
     * pairs. At construction the default is an instance of ApiLeafAtoms, which
     * generates pairs from all leaf atoms in the box.
     *
     * @param iter
     */
    public void setIterator(MoleculeIteratorAllMolecules iter) {
        iterator = iter;
    }

    /**
     * Accessor method for the iterator that generates the atom pairs used to
     * tabulate the ROG
     *
     * @return
     */
    public MoleculeIteratorAllMolecules getIterator() {
        return iterator;
    }

    public IData getData() {
        if (box == null)
            throw new IllegalStateException(
                    "must call setBox before using meter");
        Boundary boundary = box.getBoundary();
        iterator.setBox(box);
        iterator.reset();
        double[] x = data.getData();
        int nLeafAtomsTot = 0;
        double r2Tot = 0.0;
        for (IMolecule molecule = iterator.nextMolecule(); molecule != null;
             molecule = iterator.nextMolecule()) {
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
            // calculate Gyration Tensor for this chain
            double [][] GT = new double[box.getSpace().getD()][box.getSpace().getD()];
            for (int iChild = 0; iChild < childList.size(); iChild++) {
                IAtom a = childList.get(iChild);
                Vector position = a.getPosition();
                dr.Ev1Mv2(position, cm);
                boundary.nearestImage(dr);
                for(int i = 0; i < box.getSpace().getD(); i++){
                    for(int j = 0; j < box.getSpace().getD(); j++){
                        GT[i][j] += dr.getX(i)*dr.getX(j)/nLeafAtoms;
                    }
                }
            }
            Matrix A = new Matrix(GT);
            EigenvalueDecomposition EVD = new EigenvalueDecomposition(A);
            Matrix EVT = EVD.getD();
            double lamdbaX2 = EVT.get(0,0);
            double lamdbaY2 = EVT.get(1,1);
            double lamdbaZ2 = EVT.get(2,2);
            /*
            Shape descriptors from: https://en.wikipedia.org/wiki/Gyration_tensor
            x[0]: Asphericity
            x[1]: Acylindricity
            x[2]: Anisotropy
             */
            x[0] = lamdbaX2 + lamdbaY2 + lamdbaZ2;
            x[1] = (lamdbaZ2 - (lamdbaX2 + lamdbaY2)/2)/x[0];
            x[2] = (lamdbaY2 - lamdbaX2)/x[0];
            x[3] = (3.0/2.0)*((lamdbaX2*lamdbaX2 + lamdbaY2*lamdbaY2 + lamdbaZ2*lamdbaZ2)/(x[0]*x[0])) - (1.0/2.0);
        }
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    /**
     * @param box
     *            The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
    }

    private static final long serialVersionUID = 1L;
    private Box box;
    private MoleculeIteratorAllMolecules iterator;
    private final Vector cm, realPos;
    private final Vector dr;
    protected DataDoubleArray data;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;

    protected DataTag tag;
}
