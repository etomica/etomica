package etomica.atom.iterator;

import etomica.action.AtomActionAdapter;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.atom.AtomFilter;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.phase.Phase;

/**
 * Wraps an AtomIterator and filters its iterates so that only those meeting
 * specified criteria are returned. Constructor is private, and new instance is
 * created using the static makeIterator method, which attempts to return an
 * iterator that implements interface of wrapped iterator.
 */
public class AtomIteratorFiltered implements AtomIterator, java.io.Serializable {

    /**
     * Returns the iterates of the given iterator that meet the critertia of the
     * given filter. Filter is final and cannot be changed after construction.
     * 
     * @param iterator
     *            the wrapped iterator
     * @param filter
     *            iterator returns only those atoms for which this filter's
     *            accept method returns true
     */
    //private so that only way to construct is via makeIterator method
    protected AtomIteratorFiltered(AtomIterator iterator, AtomFilter filter) {
        this.iterator = iterator;
        this.filter = filter;
        actionWrapper = new ActionWrapper(filter);
    }

    /**
     * Wraps given iterator with a filtering iterator that attempts to implement
     * the same interfaces as the given iterator. Given iterator can implement
     * any combination of:
     * <ul>
     * <li>AtomIteratorAtomDependent
     * <li>AtomIteratorPhaseDependent
     * <li>AtomsetIteratorDirectable
     * <li>AtomsetIteratorTargetable
     * <li>AtomsetIteratorBasisDependent (which extends AtomsetIteratorTargetable)
     * </ul>
     * and the returned iterator will implement the same combination of
     * interfaces. Iterator may instead implement AtomIteratorListDependent, but
     * not in combination with any of the others (if it does, an exception is
     * thrown). If given iterator implements none of these interfaces, then
     * returned iterator is simply an instance of AtomIteratorFiltered.
     * 
     * @param iterator
     *            wrapped iterator
     * @param filter
     *            AtomFilter such that only atoms for which its accept method
     *            returns true will be given by the iterator. Filter is final
     *            and cannot be changed after construction.
     * @return filtering iterator that implements the interfaces of the wrapped
     *         iterator.
     * @throws IllegalArgumentException
     *             if the given iterator implements AtomIteratorListDependent
     *             and any of the interfaces: AtomIteratorAtomDependent,
     *             AtomIteratorPhaseDependent, AtomIteratorDirectable, and/or
     *             AtomIteratorTargetable
     */
    public static AtomIteratorFiltered makeIterator(AtomIterator iterator,
            AtomFilter filter) {
        int key = 0;
        if (iterator instanceof AtomIteratorAtomDependent)
            key += (1 << 0);//1
        if (iterator instanceof AtomIteratorPhaseDependent)
            key += (1 << 1);//2
        if (iterator instanceof AtomsetIteratorDirectable)
            key += (1 << 2);//4
        if (iterator instanceof AtomsetIteratorTargetable)
            key += (1 << 3);//8
        boolean basis = iterator instanceof AtomsetIteratorBasisDependent;

        switch (key) {
        case 0:
            return new AtomIteratorFiltered(iterator, filter);
        case 1:
            return new AIFAtom(iterator, filter);
        case 2:
            return new AIFPhase(iterator, filter);
        case 3:
            return new AIFAtomPhase(iterator, filter);
        case 4:
            return new AIFDirectable(iterator, filter);
        case 5:
            return new AIFAtomDirectable(iterator, filter);
        case 6:
            return new AIFPhaseDirectable(iterator, filter);
        case 7:
            return new AIFAtomPhaseDirectable(iterator, filter);
        case 8:
            return basis ? new AIFBasis(iterator, filter) : new AIFTarget(iterator, filter);
        case 9:
            return basis ? new AIFAtomBasis(iterator, filter)
                    : new AIFAtomTarget(iterator, filter);
        case 10:
            return basis ? new AIFPhaseBasis(iterator, filter)
                    : new AIFPhaseTarget(iterator, filter);
        case 11:
            return basis ? new AIFAtomPhaseBasis(iterator, filter)
                    : new AIFAtomPhaseTarget(iterator, filter);
        case 12:
            return basis ? new AIFBasisDirectable(iterator, filter)
                    : new AIFTargetDirectable(iterator, filter);
        case 13:
            return basis ? new AIFAtomBasisDirectable(iterator, filter)
                    : new AIFAtomTargetDirectable(iterator, filter);
        case 14:
            return basis ? new AIFPhaseBasisDirectable(iterator, filter)
                    : new AIFPhaseTargetDirectable(iterator, filter);
        case 15:
            return basis ? new AIFAtomPhaseBasisDirectable(iterator, filter)
                    : new AIFAtomPhaseTargetDirectable(iterator, filter);
        default:
            throw new IllegalArgumentException(
                    " not ready for that kind of iterator");
        }
    }

    /**
     * Indicates whether iterator has another iterate to return.
     */
    public boolean hasNext() {
        return next != null;
    }

    /**
     * Puts iterator in state ready for iteration.
     */
    public void reset() {
        iterator.reset();
        nextAtom();
    }

    /**
     * Sets iterator so that hasNext returns false.
     */
    public void unset() {
        iterator.unset();
        next = null;
    }

    /**
     * Returns the next atom from the iterator that meets the filter's criteria.
     */
    public IAtom nextAtom() {
        IAtom nextAtom = next;
        next = null;
        while (iterator.hasNext() && next == null) {
            next = iterator.nextAtom();
            if (!filter.accept(next))
                next = null;
        }
        return nextAtom;
    }

    /**
     * Same as nextAtom.
     */
    public AtomSet next() {
        return nextAtom();
    }

    /**
     * Returns next atom without advancing the iterator.
     */
    public AtomSet peek() {
        return next;
    }

    /**
     * Performs the given action on all atoms from iterator that meet the filter
     * criteria.
     */
    public void allAtoms(AtomsetAction action) {
        actionWrapper.action = action;
        iterator.allAtoms(actionWrapper);
    }

    /**
     * Returns the number of iterates given by the iterator that meet the
     * criteria of the filter.
     */
    public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
    }

    /**
     * Returns 1, indicating that this is an atom iterator.
     */
    public final int nBody() {
        return 1;
    }

    /**
     * Returns the iterator wrapped in this filter.
     */
    public AtomIterator getWrappedIterator() {
        return iterator;
    }

    private static final long serialVersionUID = 2L;
    protected final AtomIterator iterator;
    private final AtomFilter filter;
    private IAtom next;
    private final ActionWrapper actionWrapper;

    /**
     * Defines a new action that wraps an action such that action is performed
     * only on the atoms meeting the filter's criteria.  Used by allAtoms method.
     */
    private static class ActionWrapper extends AtomActionAdapter {

        AtomsetAction action;
        private final AtomFilter myFilter;

        public ActionWrapper(AtomFilter filter) {
            myFilter = filter;
        }

        public void actionPerformed(IAtom atom) {
            if (myFilter.accept(atom))
                action.actionPerformed(atom);
        }
    }
    
    private static class AIFAtom extends AtomIteratorFiltered implements AtomIteratorAtomDependent {
        AIFAtom(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setAtom(IAtom atom) {((AtomIteratorAtomDependent)iterator).setAtom(atom);}
    }
    private static class AIFPhase extends AtomIteratorFiltered implements AtomIteratorPhaseDependent {
        AIFPhase(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setPhase(Phase phase) {((AtomIteratorPhaseDependent)iterator).setPhase(phase);}
    }
    private static class AIFAtomPhase extends AtomIteratorFiltered implements AtomIteratorAtomDependent, AtomIteratorPhaseDependent {
        AIFAtomPhase(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setAtom(IAtom atom) {((AtomIteratorAtomDependent)iterator).setAtom(atom);}
        public void setPhase(Phase phase) {((AtomIteratorPhaseDependent)iterator).setPhase(phase);}
    }
    private static class AIFDirectable extends AtomIteratorFiltered implements AtomsetIteratorDirectable {
        AIFDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setDirection(Direction direction) {((AtomsetIteratorDirectable)iterator).setDirection(direction);}
    }
    private static class AIFAtomDirectable extends AtomIteratorFiltered implements AtomIteratorAtomDependent, AtomsetIteratorDirectable {
        AIFAtomDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setAtom(IAtom atom) {((AtomIteratorAtomDependent)iterator).setAtom(atom);}
        public void setDirection(Direction direction) {((AtomsetIteratorDirectable)iterator).setDirection(direction);}
    }
    private static class AIFPhaseDirectable extends AtomIteratorFiltered implements AtomIteratorPhaseDependent, AtomsetIteratorDirectable {
        AIFPhaseDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setPhase(Phase phase) {((AtomIteratorPhaseDependent)iterator).setPhase(phase);}
        public void setDirection(Direction direction) {((AtomsetIteratorDirectable)iterator).setDirection(direction);}
    }
    private static class AIFAtomPhaseDirectable extends AtomIteratorFiltered implements AtomIteratorAtomDependent, AtomIteratorPhaseDependent, AtomsetIteratorDirectable {
        AIFAtomPhaseDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setAtom(IAtom atom) {((AtomIteratorAtomDependent)iterator).setAtom(atom);}
        public void setPhase(Phase phase) {((AtomIteratorPhaseDependent)iterator).setPhase(phase);}
        public void setDirection(Direction direction) {((AtomsetIteratorDirectable)iterator).setDirection(direction);}
    }
    private static class AIFTarget extends AtomIteratorFiltered implements AtomsetIteratorTargetable {
        AIFTarget(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setTarget(IAtom atom) {((AtomsetIteratorTargetable)iterator).setTarget(atom);}
    }
    private static class AIFBasis extends AIFTarget implements AtomsetIteratorBasisDependent {
        AIFBasis(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public int basisSize() {return ((AtomsetIteratorBasisDependent)iterator).basisSize();}
        public void setBasis(AtomSet atoms) {((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);}
        public boolean haveTarget(IAtom atom) {return ((AtomsetIteratorBasisDependent)iterator).haveTarget(atom);}
    }
    private static class AIFAtomTarget extends AtomIteratorFiltered implements AtomIteratorAtomDependent, AtomsetIteratorTargetable {
        AIFAtomTarget(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setAtom(IAtom atom) {((AtomIteratorAtomDependent)iterator).setAtom(atom);}
        public void setTarget(IAtom atom) {((AtomsetIteratorTargetable)iterator).setTarget(atom);}
    }
    private static class AIFAtomBasis extends AIFAtomTarget implements AtomsetIteratorBasisDependent {
        AIFAtomBasis(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public int basisSize() {return ((AtomsetIteratorBasisDependent)iterator).basisSize();}
        public void setBasis(AtomSet atoms) {((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);}
        public boolean haveTarget(IAtom atom) {return ((AtomsetIteratorBasisDependent)iterator).haveTarget(atom);}
    }
    private static class AIFPhaseTarget extends AtomIteratorFiltered implements AtomIteratorPhaseDependent, AtomsetIteratorTargetable {
        AIFPhaseTarget(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setPhase(Phase phase) {((AtomIteratorPhaseDependent)iterator).setPhase(phase);}
        public void setTarget(IAtom atom) {((AtomsetIteratorTargetable)iterator).setTarget(atom);}
    }
    private static class AIFPhaseBasis extends AIFPhaseTarget implements AtomsetIteratorBasisDependent {
        AIFPhaseBasis(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public int basisSize() {return ((AtomsetIteratorBasisDependent)iterator).basisSize();}
        public void setBasis(AtomSet atoms) {((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);}
        public boolean haveTarget(IAtom atom) {return ((AtomsetIteratorBasisDependent)iterator).haveTarget(atom);}
    }
    private static class AIFAtomPhaseTarget extends AtomIteratorFiltered implements AtomIteratorAtomDependent, AtomIteratorPhaseDependent, AtomsetIteratorTargetable {
        AIFAtomPhaseTarget(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setAtom(IAtom atom) {((AtomIteratorAtomDependent)iterator).setAtom(atom);}
        public void setPhase(Phase phase) {((AtomIteratorPhaseDependent)iterator).setPhase(phase);}
        public void setTarget(IAtom atom) {((AtomsetIteratorTargetable)iterator).setTarget(atom);}
    }
    private static class AIFAtomPhaseBasis extends AIFAtomPhaseTarget implements AtomsetIteratorBasisDependent {
        AIFAtomPhaseBasis(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public int basisSize() {return ((AtomsetIteratorBasisDependent)iterator).basisSize();}
        public void setBasis(AtomSet atoms) {((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);}
        public boolean haveTarget(IAtom atom) {return ((AtomsetIteratorBasisDependent)iterator).haveTarget(atom);}
    }
    private static class AIFTargetDirectable extends AtomIteratorFiltered implements AtomsetIteratorTargetable, AtomsetIteratorDirectable {
        AIFTargetDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setTarget(IAtom atom) {((AtomsetIteratorTargetable)iterator).setTarget(atom);}
        public void setDirection(Direction direction) {((AtomsetIteratorDirectable)iterator).setDirection(direction);}
    }
    private static class AIFBasisDirectable extends AIFTargetDirectable implements AtomsetIteratorBasisDependent {
        AIFBasisDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public int basisSize() {return ((AtomsetIteratorBasisDependent)iterator).basisSize();}
        public void setBasis(AtomSet atoms) {((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);}
        public boolean haveTarget(IAtom atom) {return ((AtomsetIteratorBasisDependent)iterator).haveTarget(atom);}
    }
    private static class AIFAtomTargetDirectable extends AtomIteratorFiltered implements AtomIteratorAtomDependent, AtomsetIteratorTargetable, AtomsetIteratorDirectable {
        AIFAtomTargetDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setAtom(IAtom atom) {((AtomIteratorAtomDependent)iterator).setAtom(atom);}
        public void setTarget(IAtom atom) {((AtomsetIteratorTargetable)iterator).setTarget(atom);}
        public void setDirection(Direction direction) {((AtomsetIteratorDirectable)iterator).setDirection(direction);}
    }
    private static class AIFAtomBasisDirectable extends AIFAtomTargetDirectable implements AtomsetIteratorBasisDependent {
        AIFAtomBasisDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public int basisSize() {return ((AtomsetIteratorBasisDependent)iterator).basisSize();}
        public void setBasis(AtomSet atoms) {((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);}
        public boolean haveTarget(IAtom atom) {return ((AtomsetIteratorBasisDependent)iterator).haveTarget(atom);}
    }
    private static class AIFPhaseTargetDirectable extends AtomIteratorFiltered implements AtomIteratorPhaseDependent, AtomsetIteratorTargetable, AtomsetIteratorDirectable {
        AIFPhaseTargetDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setPhase(Phase phase) {((AtomIteratorPhaseDependent)iterator).setPhase(phase);}
        public void setTarget(IAtom atom) {((AtomsetIteratorTargetable)iterator).setTarget(atom);}
        public void setDirection(Direction direction) {((AtomsetIteratorDirectable)iterator).setDirection(direction);}
    }
    private static class AIFPhaseBasisDirectable extends AIFPhaseTargetDirectable implements AtomsetIteratorBasisDependent {
        AIFPhaseBasisDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public int basisSize() {return ((AtomsetIteratorBasisDependent)iterator).basisSize();}
        public void setBasis(AtomSet atoms) {((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);}
        public boolean haveTarget(IAtom atom) {return ((AtomsetIteratorBasisDependent)iterator).haveTarget(atom);}
    }
    private static class AIFAtomPhaseTargetDirectable extends AtomIteratorFiltered implements AtomIteratorAtomDependent, AtomIteratorPhaseDependent, AtomsetIteratorTargetable, AtomsetIteratorDirectable {
        AIFAtomPhaseTargetDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public void setAtom(IAtom atom) {((AtomIteratorAtomDependent)iterator).setAtom(atom);}
        public void setPhase(Phase phase) {((AtomIteratorPhaseDependent)iterator).setPhase(phase);}
        public void setTarget(IAtom atom) {((AtomsetIteratorTargetable)iterator).setTarget(atom);}
        public void setDirection(Direction direction) {((AtomsetIteratorDirectable)iterator).setDirection(direction);}
    }
    private static class AIFAtomPhaseBasisDirectable extends AIFAtomPhaseTargetDirectable implements AtomsetIteratorBasisDependent {
        AIFAtomPhaseBasisDirectable(AtomIterator iterator, AtomFilter filter) {super(iterator, filter);}
        public int basisSize() {return ((AtomsetIteratorBasisDependent)iterator).basisSize();}
        public void setBasis(AtomSet atoms) {((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);}
        public boolean haveTarget(IAtom atom) {return ((AtomsetIteratorBasisDependent)iterator).haveTarget(atom);}
    }

}
