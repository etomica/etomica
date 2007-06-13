package etomica.nbr;

import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;

/**
 * A PotentialGroup class to work with PotentialMasterHybrid.  Holds a reference
 * to cell, list and hybrid PotentialMasters so that all of them can be notified
 *
 * @author Andrew Schultz
 */
public class PotentialGroupHybrid extends PotentialGroupNbr {

    public PotentialGroupHybrid(int nBody, Space space) {
        super(nBody, space);
    }

    public void setPotentialMaster(PotentialMaster newPotentialMaster) {
        // we'll get notification of each PotentialMaster instance (cell, list, 
        // hybrid), but we only care about the hyrbid one.  It will notify the
        // others.
        if (newPotentialMaster instanceof PotentialMasterHybrid) {
            // we got what we wanted.  
            super.setPotentialMaster(newPotentialMaster);
        }
    }
    
    private static final long serialVersionUID = 1L;
    protected PotentialMasterCell potentialMasterCell;
    protected PotentialMasterList potentialMasterList;
}
