/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.molecule.IMolecule;

/**
 * Encapsulation of a set of instructions that an AtomsetIterator
 * uses to select the atoms it presents on iteration.
 *
 * @author David Kofke and Andrew Schultz
 */

public final class IteratorDirective {

    /**
     * Flag indicating whether long-range correction contributions should
     * be included in calculation.  Default is <b>true</b>.
     */
    public boolean includeLrc = true;
    private Direction direction;
    private IAtom targetAtom = null;
    private IMolecule targetMolecule = null;

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
        includeLrc = false;
    }

    public void setDirection(Direction direction) {
        this.direction = direction;
    }

    public Direction direction() {
        return direction;
    }

    public IAtom getTargetAtom() {
        return targetAtom;
    }

    public void setTargetAtom(IAtom atom) {
        targetAtom = atom;
        targetMolecule = null;
    }

    public IMolecule getTargetMolecule() {
        return targetMolecule;
    }

    public void setTargetMolecule(IMolecule mole) {
        targetMolecule = mole;
        targetAtom = null;
    }

    /**
     * Sets flag indicating if lrc potentials (long-range correction) should be
     * included.
     *
     * @return boolean
     */
    public boolean isIncludeLrc() {
        return includeLrc;
    }

    /**
     * Sets flag indicating if lrc potentials (long-range correction) should be
     * included.
     *
     * @param includeLrc The flag value to set
     * @return this IteratorDirective, for in-line use of the method.
     */
    public IteratorDirective setIncludeLrc(boolean includeLrc) {
        this.includeLrc = includeLrc;
        return this;
    }

    //IteratorDirective.Direction
    public enum Direction {
        UP,
        DOWN
    }
}    
