/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.thermointegration;

import etomica.action.ResetAccumulators;
import etomica.data.DataPump;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.Integrator;

import java.util.ArrayList;


/**
 * Abstract class suitable as a parent class for thermodynamic integration
 * integrators.  Subclasses should implement advanceParameters to advance
 * the thermodynamic parameter (or parameters) of interest.  Data collection
 * can be performed as interval actions.
 *
 * @author Andrew Schultz
 */
public abstract class IntegratorThermodynamic extends Integrator {

    /**
     * Constructs a thermodynamic integrator that operates with the given
     * sub integrator, which is responsible for exploring configruations at any
     * given state point.
     */
    public IntegratorThermodynamic(Integrator subIntegrator) {
        super();
        dataPumps = new ArrayList<DataPump>();
        this.subIntegrator = subIntegrator;
        setNumSubsteps(1000);
        setNumEquilibrationSubsteps(10);
    }

    /**
     * Resets this integrator and passes on the reset to the subIntegrator.
     * Thermodynamic integration will start back at the initial parameter(s).
     * An data collected will be reset
     */
    public void reset() throws ConfigurationOverlapException {
        ConfigurationOverlapException overlapException = null;
        try {
            subIntegrator.reset();
        }
        catch (ConfigurationOverlapException e) {
            if (overlapException == null) {
                overlapException = e;
            }
        }
        if (overlapException != null) {
            throw overlapException;
        }
        
        iSubstep = 0;
        
        ResetAccumulators resetAccumulators = new ResetAccumulators();
        for (int i=0; i<dataPumps.size(); i++) {
            resetAccumulators.actionPerformed(dataPumps.get(i));
        }
        
        super.reset();
    }

    /**
     * Returns the subIntegrator, responsible for actually walking through
     * configurational space.
     */
    public Integrator getSubIntegrator() {
        return subIntegrator;
    }

    /**
     * Performs a step in the subIntegrator.  If numSubsteps have been
     * performed since the last advance of the thermodynamic parameter(s),
     * the parameters are advanced again.
     */
    protected void doStepInternal() {
        iSubstep++;
        subIntegrator.doStep();
        if (iSubstep > numEquilibrationSubsteps) {
            for (int i=0; i<dataPumps.size(); i++) {
                dataPumps.get(i).actionPerformed();
            }
        }
        if (iSubstep == numSubsteps) {
            advanceParameters();
            iSubstep = 0;
        }
    }
    
    /**
     * Method to advance the thermodynamic parameter(s) of interest as part
     * of the thermodynamic integration.  The subclass implements this method.
     * The subclass is also responsible for pulling out whatever data it wants
     * (from the DataPumps, or elsewhere) and storing that.
     */
    protected abstract void advanceParameters();

    /**
     * Sets the number of steps performed at each thermodynamic state point.
     */
    public void setNumSubsteps(int newNumSubsteps) {
        numSubsteps = newNumSubsteps;
        if (iSubstep >= numSubsteps) {
            iSubstep = numSubsteps - 1;
        }
    }

    /**
     * Returns the number of steps performed at each thermodynamic state point.
     */
    public int getNumSubsteps() {
        return numSubsteps;
    }
    
    /**
     * Sets the number of steps equilibration steps performed at each
     * thermodynamic state point before data collection begins.
     */
    public void setNumEquilibrationSubsteps(int newNumEquilibrationSubsteps) {
        numEquilibrationSubsteps = newNumEquilibrationSubsteps;
    }

    /**
     * Returns the number of steps equilibration steps performed at each
     * thermodynamic state point before data collection begins.
     */
    public int getNumEquilibrationSubsteps() {
        return numEquilibrationSubsteps;
    }
    
    public ArrayList<DataPump> getDataPumps() {
        return dataPumps;
    }
    
    private static final long serialVersionUID = 1L;
    protected Integrator subIntegrator;
    protected int numSubsteps;
    protected int iSubstep;
    protected int numEquilibrationSubsteps;
    protected ArrayList<DataPump> dataPumps;
}
