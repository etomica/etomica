package etomica.data.meter;
import etomica.Atom;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Space;
import etomica.Space.CoordinatePair;
import etomica.Space.Vector;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.units.Dimension;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).
 *
 * @author David Kofke
 */
public class MeterRadiusGyration extends MeterScalar implements EtomicaElement {
	
	/**
	 * Creates meter with default to compute pair correlation for all
	 * leaf atoms in a phase.
	 * @param parent
	 */
    public MeterRadiusGyration(Space space) {
	    super();
	    setLabel("Rg");
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
     * Mutator method for the iterator that generates the atom pairs used to tabulate
     * the radial distribution function.  By setting this iterator the meter can be
     * configured to compute pair distribution for any set of atom pairs.  At construction
     * the default is an instance of ApiLeafAtoms, which generates pairs from all leaf
     * atoms in the phase.
     * @param iter
     */
    public void setIterator(AtomIteratorAllMolecules iter) {iterator = iter;}
    
    /**
     * Accessor method for the iterator that generates the atom pairs used to
     * tabulate the radial distribution function.
     * @return
     */
    public AtomIteratorAllMolecules getIterator() {return iterator;}
    
    /**
     * Returns Dimension.NULL, indicating that the measured quantity is dimensionless.
     */
    public Dimension getDimension() {return Dimension.LENGTH;}

	/**
	 * Computes RDF for the current configuration of the given phase.
	 */
	public double getDataAsScalar(Phase aPhase) {
        cPair.setNearestImageTransformer(aPhase.boundary());
		iterator.setPhase(aPhase);
	    iterator.reset();
        AtomIteratorTree leafIterator = new AtomIteratorTree();
        int nLeafAtomsTot = 0;
        double r2Tot = 0.0;
	    while(iterator.hasNext()) {
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
                cPair.reset(prevAtom.coord,a.coord);
                Space.Vector dr = cPair.dr();
                realPos.PE(dr);
                cm.PE(realPos);
                prevAtom = a;
            }
            cm.DE(nLeafAtoms);
            // calculate Rg^2 for this chain
            double r2 = 0.0;
            leafIterator.reset();
            prevAtom = leafIterator.nextAtom();
            realPos.E(prevAtom.coord.position());
            while (leafIterator.hasNext()) {
                Atom a = leafIterator.nextAtom();
                cPair.reset(prevAtom.coord,a.coord);
                Space.Vector dr = cPair.dr();
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
	
    private AtomIteratorAllMolecules iterator;
    private final Space.CoordinatePair cPair;
    private final Space.Vector cm, realPos;
	
}
    