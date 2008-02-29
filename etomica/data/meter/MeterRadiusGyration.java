package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.DataSourceScalar;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.units.Length;

/**
 * Meter for tabulation of the radius of gyration of a set of chain molecules. 
 * 
 * @author David Kofke
 */
public class MeterRadiusGyration extends DataSourceScalar {

    public MeterRadiusGyration(Space space) {
        super("Radius of Gyration", Length.DIMENSION);
        iterator = new AtomIteratorAllMolecules();
        cm = space.makeVector();
        realPos = space.makeVector();
        dr = space.makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Calculates radius of gyration");
        return info;
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
    public void setIterator(AtomIteratorAllMolecules iter) {
        iterator = iter;
    }

    /**
     * Accessor method for the iterator that generates the atom pairs used to
     * tabulate the ROG
     * 
     * @return
     */
    public AtomIteratorAllMolecules getIterator() {
        return iterator;
    }

    public double getDataAsScalar() {
        if (box == null)
            throw new IllegalStateException(
                    "must call setBox before using meter");
        NearestImageTransformer nearestImageTransformer = box.getBoundary();
        iterator.setBox(box);
        iterator.reset();
        int nLeafAtomsTot = 0;
        double r2Tot = 0.0;
        for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            // loop over molecules
            AtomSet childList = ((IMolecule)iterator.nextAtom()).getChildList();
            if (childList.getAtomCount() < 2) {
                // a monatomic molecule
                continue;
            }

            // find center of mass
            //do the first iterate explicitly, assume there is at least
            // one leaf atom
            IAtomPositioned firstAtom = (IAtomPositioned)childList.getAtom(0);
            int nLeafAtoms = 1;
            realPos.E(firstAtom.getPosition());
            cm.E(realPos);
            IVector prevPosition = firstAtom.getPosition();
            for (int iChild = 1; iChild < childList.getAtomCount(); iChild++) {
                IAtomPositioned a = (IAtomPositioned)childList.getAtom(iChild);
                nLeafAtoms++;
                IVector position = a.getPosition();
                dr.Ev1Mv2(position, prevPosition);
                //molecule might be wrapped around the box.  calculate
                //the real difference in position
                nearestImageTransformer.nearestImage(dr);
                //realPos is now the effective position of a
                realPos.PE(dr);
                cm.PE(realPos);
                prevPosition = position;
            }
            cm.TE(1.0 / nLeafAtoms);
            // calculate Rg^2 for this chain
            double r2 = 0.0;
            realPos.E(firstAtom.getPosition());
            for (int iChild = 1; iChild < childList.getAtomCount(); iChild++) {
                IAtomPositioned a = (IAtomPositioned)childList.getAtom(iChild);
                IVector position = a.getPosition();
                dr.Ev1Mv2(position, prevPosition);
                //molecule might be wrapped around the box.  calculate
                //the real difference in position
                nearestImageTransformer.nearestImage(dr);
                //realPos is now the effective position of a
                realPos.PE(dr);
                dr.Ev1Mv2(realPos, cm);// = realPos.M(cm);
                r2 += dr.squared();
                prevPosition = position;
            }
            r2Tot += r2;
            nLeafAtomsTot += nLeafAtoms;
        }
        return r2Tot / nLeafAtomsTot;
    }

    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }

    /**
     * @param box
     *            The box to set.
     */
    public void setBox(IBox box) {
        this.box = box;
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
    private AtomIteratorAllMolecules iterator;
    private final IVector cm, realPos;
    private final IVector dr;

}