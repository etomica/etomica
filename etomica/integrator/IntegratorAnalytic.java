package etomica.integrator;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.exception.ConfigurationOverlapException;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.util.IRandom;

/**
 * Integrator that generates atom trajectories from an analytic formula.
 * Takes an IntegratorAnalytic.AtomAction instance and performs action on all atoms in each
 * time step; intention is for the action to set the position of the atom for the
 * current time (but it could do anything).
 *
 * @author David Kofke
 * @author Nancy Cribbin
 */
public class IntegratorAnalytic extends IntegratorMD {
    
    private static final long serialVersionUID = 1L;
    private AtomAction action;
    
    public IntegratorAnalytic(Simulation sim) {
        this(sim.getPotentialMaster(),sim.getRandom(),sim.getDefaults().timeStep);
    }
    
    public IntegratorAnalytic(PotentialMaster potentialMaster, IRandom random,
                              double timeStep) {
        super(potentialMaster,random,timeStep,0);
    }
    
    public void doStepInternal() {
        if(action == null) return;
        elapsedTime += getTimeStep();
        action.setTime(elapsedTime);
        AtomArrayList leafList = phase.getSpeciesMaster().getLeafList();
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            AtomLeaf a = (AtomLeaf)leafList.get(iLeaf);
            action.actionPerformed(a);
        }
    }
    
    public void reset() throws ConfigurationOverlapException {
        elapsedTime = 0.0;
        super.reset();
    }
    
    public void setAction(AtomAction action) {this.action = action;}
    
    public AtomAction getAction() {return action;}
    
    private double elapsedTime = 0.0;
    
    /**
     * Extends AtomAction class to add a method to set the time.
     */
    public static abstract class AtomAction extends etomica.action.AtomActionAdapter {
        protected double time;
        public void setTime(double t) {time = t;}
    }
    
 }//end of IntegratorAnalytic
 
