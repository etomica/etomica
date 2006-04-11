/*
 * History
 * Created on Nov 27, 2004 by kofke
 */
package etomica.nbr;


/**
 * Specifies that all atoms pairs are to be considered neighbors, but still 
 * claims to be range-dependent.  This is useful for Cell-based neighbors since
 * the cell listing provides the range depedence.  Should
 * not be used for species in which atoms are being added/removed by integrator.
 */
public class CriterionRangedAll extends CriterionAll {

    /**
     * Returns true.
     */
    public boolean isRangeDependent() {
        return true;
    }
}
