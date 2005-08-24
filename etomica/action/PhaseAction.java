package etomica.action;

import etomica.Phase;

 /**
  * Elementary action performed on a phase.
  */
public interface PhaseAction extends Action {

    public void setPhase(Phase phase);
    
    public Phase getPhase();

}