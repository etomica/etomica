package etomica.potential;

import etomica.AtomType;


/**
 * Interface for a potential that is artificially truncated, and thus which can
 * provide a zero-body potential that (approximately) corrects for the
 * truncation. The PotentialMaster gets a LRC potential from a potential
 * implementing this interface when the potential is given to the
 * PotentialMaster's setSpecies method.  PotentialGroup also recognizes
 * this interface and adds a LRC potential to the PotentialMaster when a
 * PotentialTruncated is added to it via the type-specifying addPotential method.
 * 
 * @see PotentialMaster
 * @see PotentialGroup
 */

/*
 * History
 * Created on Mar 28, 2005
 */
public interface PotentialTruncated {

    /**
     * Returns a class that calculates the long-range contribution to the potential
     * that becomes neglected by the truncation.  May return null.
     */
    public Potential0Lrc makeLrcPotential(AtomType[] types);
}
