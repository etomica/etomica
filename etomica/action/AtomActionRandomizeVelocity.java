package etomica.action;

import etomica.Atom;
import etomica.space.ICoordinateKinetic;
import etomica.space.Vector;
import etomica.statmech.MaxwellBoltzmann;


/**
 * Action that sets the velocity vector of a given atom to a randomly
 * chosen value sampled from a Boltzmann distribution.
 * 
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jan 27, 2005 by kofke
 */
public class AtomActionRandomizeVelocity extends AtomActionAdapter {

    public AtomActionRandomizeVelocity() {
        setLabel("Randomize velocity");
    }

    /* (non-Javadoc)
     * @see etomica.action.AtomAction#actionPerformed(etomica.Atom)
     */
    public void actionPerformed(Atom a) {
        Vector velocity = ((ICoordinateKinetic)a.coord).velocity();
        double mass = a.type.getMass();
        int D = velocity.D();
        for(int i=0; i<D; i++) {
            velocity.setX(i,MaxwellBoltzmann.randomMomentumComponent(temperature,mass)/mass);
        }
    }
    
    /**
     * @return Returns the Boltzmann-distribution temperature used to sample
     * the new velocity.
     */
    public double getTemperature() {
        return temperature;
    }
    /**
     * @param temperature The temperature to set.
     */
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }
    private double temperature;

}
