package etomica.potential;

import etomica.api.IAtomType;
import etomica.api.IPotential;
import etomica.space.Space;

/**
 * 2-body soft truncated Potential class for use between two monatomic
 * molecules.  The 0-body LRC potential created is that for the wrapped leaf
 * potential, although the AtomTypes passed are molecule types.  This should be
 * fine since the LRC potential only cares about how many atoms of the given
 * type(s) there are.
 *
 * @author Andrew Schultz
 */
public class P2SoftMoleculeMonatomicTruncated extends P2SoftMoleculeMonatomic
        implements PotentialTruncated {

    public P2SoftMoleculeMonatomicTruncated(Space space, IPotential potential) {
        super(space, potential);
    }

    public Potential0Lrc makeLrcPotential(IAtomType[] types) {
        return ((PotentialTruncated)wrappedPotential).makeLrcPotential(types);
    }

    private static final long serialVersionUID = 1L;
}
