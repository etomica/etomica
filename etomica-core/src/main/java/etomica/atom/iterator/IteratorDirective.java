/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.util.EnumeratedType;

/**
 * Encapsulation of a set of instructions that an AtomsetIterator
 * uses to select the atoms it presents on iteration.
 *
 * @author David Kofke and Andrew Schultz
 */

public final class IteratorDirective implements java.io.Serializable {
    
    public IteratorDirective() {
        this(Direction.UP);
    }
    public IteratorDirective(Direction direction) {
        setDirection(direction);
    }
    public IteratorDirective(Direction direction, IAtom atom) {
        this(direction);
        targetAtom = atom;
        targetMolecule = null;
    }
    public IteratorDirective(Direction direction, IMolecule mole) {
        this(direction);
        targetAtom = null;
        targetMolecule = mole;
    }

    /**
     * Puts directive in default state of no atoms specified, up direction, no
     * potential criteria applied, no LRC included.
     */
    public void clear() {
        setDirection(Direction.UP);
        targetAtom = null;
        targetMolecule = null;
        potentialCriteriaHead = null;
        includeLrc = false;
    }
    
    public final void setDirection(Direction direction) {
        this.direction = direction;
    }

    public final Direction direction() {return direction;}
    
    public final boolean excludes(IPotential p) {
        for(PotentialCriterion crit=potentialCriteriaHead; crit!=null; crit=crit.nextCriterion()) {
            if(crit.excludes(p)) return true;
        }
        return false;
    }
    
    /**
     * Adds criterion to set of criteria for potentials.  Criteria
     * are considered in the reverse order of their addition (last-added is considered first).
     * There is no way to remove a criterion.
     */
    public final void addCriterion(PotentialCriterion newCriterion) {
        if(newCriterion == null) return;
        newCriterion.setNextCriterion(potentialCriteriaHead);
        potentialCriteriaHead = newCriterion;
    }

    public void setTargetAtom(IAtom atom) {
        targetAtom = atom;
        targetMolecule = null;
    }

    public void setTargetMolecule(IMolecule mole) {
        targetMolecule = mole;
        targetAtom = null;
    }

    public IAtom getTargetAtom() {return targetAtom;}
    public IMolecule getTargetMolecule() {return targetMolecule;}

    private static final long serialVersionUID = 1L;
    /**
     * Flag indicating whether long-range correction contributions should
     * be included in calculation.  Default is <b>true</b>.
     */
    public boolean includeLrc = true;
    private Direction direction;
    PotentialCriterion potentialCriteriaHead;
    private IAtom targetAtom = null;
    private IMolecule targetMolecule = null;
    
    //IteratorDirective.Direction
    public static final class Direction extends EnumeratedType {
        
        private static final long serialVersionUID = 1L;

        private Direction(String label) {
            super(label);
        }

        public static final Direction UP = new Direction("Up");
        public static final Direction DOWN = new Direction("Down");

        public static EnumeratedType[] choices() {return new Direction[] {UP,DOWN,null};}

        /**
         * Required to guarantee singleton when deserializing.
         * @return the singleton INSTANCE
         */
        private Object readResolve() {
            return this.toString().equals("Up") ? Direction.UP : Direction.DOWN;
        }
        
    }//end of Direction
    
    /**
     * Class used to define a criterion that must be satisfied by a potential
     * in order for its atoms to be approved for iteration by an iterator directive.
     * Multiple criteria are ordered into a linked list by the iterator directive.
     * This is made cloneable to support IteratorDirective.copy functionality.
     */
    public static abstract class PotentialCriterion implements Cloneable, java.io.Serializable {
        /**
         * Definition of criterion.  If this method returns true, the potential's atoms
         * are excluded from iteration.
         * Note that any subclasses should be sure to override clone method if more than
         * a shallow copy is appropriate.
         */
        public abstract boolean excludes(IPotential pot);
        
        //Linked-list constructs
        private PotentialCriterion nextCriterion;
        public void setNextCriterion(PotentialCriterion next) {nextCriterion = next;}
        public PotentialCriterion nextCriterion() {return nextCriterion;}
        
        public Object clone() {
            Object obj = null;
            try {
                obj = super.clone();
            } catch(CloneNotSupportedException e) {e.printStackTrace();}
            return obj;
        }
    }//end of PotentialCriterion
    
    /**
     * Sets flag indicating if lrc potentials (long-range correction) should be
     * included.
     * @return boolean
     */
    public boolean isIncludeLrc() {
        return includeLrc;
    }

    /**
     * Sets flag indicating if lrc potentials (long-range correction) should be
     * included.
     * @param includeLrc The flag value to set
     * @return this IteratorDirective, for in-line use of the method.
     */
    public IteratorDirective setIncludeLrc(boolean includeLrc) {
        this.includeLrc = includeLrc;
        return this;
    }
}    
