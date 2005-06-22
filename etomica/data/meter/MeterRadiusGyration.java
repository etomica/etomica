package etomica.data.meter;

import etomica.Atom;
import etomica.DataInfo;
import etomica.EtomicaInfo;
import etomica.Meter;
import etomica.Phase;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.data.DataSourceScalar;
import etomica.space.CoordinatePair;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).
 * 
 * @author David Kofke
 */
public class MeterRadiusGyration extends DataSourceScalar implements Meter {

    /**
     * Creates meter with default to compute pair correlation for all leaf atoms
     * in a phase.
     * 
     * @param parent
     */
    public MeterRadiusGyration(Space space) {
        super(new DataInfo("Radius of Gyration", Dimension.LENGTH));
        iterator = new AtomIteratorAllMolecules();
        cPair = space.makeCoordinatePair();
        cm = space.makeVector();
        realPos = space.makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Calculates radius of gyration");
        return info;
    }

    /**
     * Mutator method for the iterator that generates the atom pairs used to
     * tabulate the radial distribution function. By setting this iterator the
     * meter can be configured to compute pair distribution for any set of atom
     * pairs. At construction the default is an instance of ApiLeafAtoms, which
     * generates pairs from all leaf atoms in the phase.
     * 
     * @param iter
     */
    public void setIterator(AtomIteratorAllMolecules iter) {
        iterator = iter;
    }

    /**
     * Accessor method for the iterator that generates the atom pairs used to
     * tabulate the radial distribution function.
     * 
     * @return
     */
    public AtomIteratorAllMolecules getIterator() {
        return iterator;
    }

    /**
     * Computes RDF for the current configuration of the given phase.
     */
    public double getDataAsScalar() {
        if (phase == null)
            throw new IllegalStateException(
                    "must call setPhase before using meter");
        cPair.setNearestImageTransformer(phase.boundary());
        iterator.setPhase(phase);
        iterator.reset();
        AtomIteratorTree leafIterator = new AtomIteratorTree();
        int nLeafAtomsTot = 0;
        double r2Tot = 0.0;
        while (iterator.hasNext()) {
            // loop over molecules
            leafIterator.setRoot(iterator.nextAtom());
            leafIterator.reset();
            if (!leafIterator.hasNext()) {
                continue;
            }
            // find center of mass
            Atom prevAtom = leafIterator.nextAtom();
            cm.E(prevAtom.coord.position());
            int nLeafAtoms = 1;
            realPos.E(prevAtom.coord.position());
            while (leafIterator.hasNext()) {
                nLeafAtoms++;
                Atom a = leafIterator.nextAtom();
                cPair.reset(prevAtom.coord, a.coord);
                Vector dr = cPair.dr();
                realPos.PE(dr);
                cm.PE(realPos);
                prevAtom = a;
            }
            cm.TE(1.0 / nLeafAtoms);
            // calculate Rg^2 for this chain
            double r2 = 0.0;
            leafIterator.reset();
            prevAtom = leafIterator.nextAtom();
            realPos.E(prevAtom.coord.position());
            while (leafIterator.hasNext()) {
                Atom a = leafIterator.nextAtom();
                cPair.reset(prevAtom.coord, a.coord);
                Vector dr = cPair.dr();
                realPos.PE(dr);
                dr = realPos.M(cm);
                r2 += dr.squared();
                prevAtom = a;
            }
            r2Tot += r2;
            nLeafAtomsTot += nLeafAtoms;
        }
        return r2Tot / nLeafAtomsTot;
    }

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }

    /**
     * @param phase
     *            The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    private AtomIteratorAllMolecules iterator;
    private final CoordinatePair cPair;
    private final Vector cm, realPos;

}