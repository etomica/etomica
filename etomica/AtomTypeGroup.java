/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica;


/**
 * Type for an AtomGroup atom.
 */
public class AtomTypeGroup extends AtomType {
    public /*final*/ boolean childrenAreGroups;
    public AtomTypeGroup(AtomFactory creator) {//presently unable to determine other fields when group is constructed, so don't make them final but set after group is built
        super(creator);
    }
    public AtomTypeGroup(AtomFactory creator,
                  boolean childrenAreGroups) {
        super(creator);
        this.childrenAreGroups = childrenAreGroups;
    }
}