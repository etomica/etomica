package etomica.surfacetension;

import etomica.action.IAction;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.exception.ConfigurationOverlapException;
import etomica.nbr.cell.PotentialMasterCell;

public class ActionRecenter implements IAction {
    
    private final LJMC sim;

    public ActionRecenter(LJMC sim) {
        this.sim = sim;
    }

    public void actionPerformed() {

        double L = sim.box.getBoundary().getBoxSize().getX(0);

        // calculate structure factor and phase angle for lowest-frequency
        // concentration wave (delta rho (x)).
        
        IAtomList leafAtoms = sim.box.getLeafList();
        int nTot = leafAtoms.getAtomCount();
        double sumCos = 0, sumSin = 0;
        double q = 2*Math.PI/L;
        for (int i=0; i<nTot; i++) {
            Vector pos = leafAtoms.getAtom(i).getPosition();
            double qx = q*pos.getX(0);
            double sinx = Math.sin(qx);
            double cosx = Math.cos(qx);
            sumCos += cosx;
            sumSin += sinx;
        }
        double amplitude = (2*Math.sqrt((sumCos*sumCos+sumSin*sumSin))/nTot);
        // concentration wave amplitude must be large enough to correspond
        // to a substantial segregation (phase separation).
        if (amplitude < 0.75) {
//          System.out.println("not centering "+amplitude);
            return;
        }
//      System.out.println("centering "+amplitude);
        // phase angle = atan(sumCos/sumSin), where delta rho = 0
        double center = -(Math.atan2(sumCos, sumSin)/q-0.25*L);
        if (center > 1) {
            center = 1;
        }
        else if (center < -1) {
            center = -1;
        }
        for (int i=0; i<nTot; i++) {
            Vector pos = leafAtoms.getAtom(i).getPosition();
            pos.setX(0, pos.getX(0) - center);
        }
        ((PotentialMasterCell)sim.integrator.getPotentialMaster()).getNbrCellManager(sim.box).assignCellAll();
        try {
            sim.integrator.reset();
        }
        catch (ConfigurationOverlapException e) {
            // we can cause overlap by increasing tail diameter
            // if so, that might not have been resolved yet
        }
    }
}
