/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Standard Monte Carlo atom-displacement trial move.  Two atoms are moved at a
 * time in such a way that the geometric center of the system is not changed.
 * 
 * for sampling the Bennet's Overlap
 *
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MCMoveAtomCoupledBennet extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtomList;
    protected final MeterPotentialEnergy energyMeter;
    protected final MeterHarmonicEnergy harmonicEnergyMeter;
    protected final Vector translationVector;
    protected IAtom atom0, atom1;
    protected double uOld, uNew;
    protected double gamma_Old, gamma_New;
    protected AtomSource atomSource;
    protected final IRandom random;
    protected final AtomPair pair;
    private double refPref;
    public NormalModes normalModes;
    public CoordinateDefinition coordinateDefinition;
    public double temperature;
    public double latticeEnergy;



	public MCMoveAtomCoupledBennet(PotentialMaster potentialMaster, IRandom random,
                                   CoordinateDefinition coordinateDef, NormalModes normalMode, double refpref, Space _space) {
        super(potentialMaster);
        this.random = random;
        this.refPref = refpref;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setCoordinateDefinition(coordinateDef);
        setNormalModes(normalMode);
        harmonicEnergyMeter = new MeterHarmonicEnergy(coordinateDefinition, normalModes);
        
        translationVector = _space.makeVector();
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.1);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        affectedAtomList = new AtomArrayList(2);
        affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtomList);
        pair = new AtomPair();
    }
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom0 = atomSource.getAtom();
        atom1 = atomSource.getAtom();
        if (atom0 == null || atom1 == null || atom0 == atom1) return false;
       
        uOld = energyMeter.getDataAsScalar() - latticeEnergy;
        double exp_uOld = Math.exp(-uOld /temperature);
        double exp_uOldHarmonic = Math.exp(-harmonicEnergyMeter.getDataAsScalar() /temperature);
        gamma_Old = (exp_uOld*exp_uOldHarmonic) / (exp_uOld + refPref*exp_uOldHarmonic);
        
//        System.out.println("exp_uOld: "+ exp_uOld);
//        System.out.println("exp_uOldHarmonic: "+exp_uOldHarmonic);
//        System.out.println("gamma_Old: "+ gamma_Old);
        
        pair.atom0 = atom0;
        pair.atom1 = atom1;
        
        if( 0 > gamma_Old ) {
            throw new ConfigurationOverlapException(box);
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom0.getPosition().PE(translationVector);
        atom1.getPosition().ME(translationVector);

        uNew = energyMeter.getDataAsScalar() - latticeEnergy;
        double exp_uNew = Math.exp(-uNew /temperature);
        double exp_uNewHarmonic = Math.exp(-harmonicEnergyMeter.getDataAsScalar()/ temperature);
        gamma_New = (exp_uNew*exp_uNewHarmonic) / (exp_uNew + refPref*exp_uNewHarmonic);
        
//        System.out.println("\nexp_uNew: "+ exp_uNew);
//        System.out.println("uexp_NewHarmonic: "+exp_uNewHarmonic);
//        System.out.println("gamma_New: "+ gamma_New+ "\n");
        
        return true;
    }//end of doTrial
    
    
    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double getA() {return gamma_New / gamma_Old;}
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double getB() {
        return 0;
    }
    
    public double energyChange() {
    	return uNew - uOld;
    	}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
        atom0.getPosition().ME(translationVector);
        atom1.getPosition().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomList.clear();
        affectedAtomList.add(atom0);
        affectedAtomList.add(atom1);
        return affectedAtomIterator;
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        atomSource.setBox(p);
    }
    
    /**
     * @return Returns the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
    /**
     * @param atomSource The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }

	public NormalModes getNormalModes() {
		return normalModes;
	}

	public void setNormalModes(NormalModes normalModes) {
		this.normalModes = normalModes;
	}

	public CoordinateDefinition getCoordinateDefinition() {
		return coordinateDefinition;
	}

	public void setCoordinateDefinition(CoordinateDefinition coordinateDefinition) {
		this.coordinateDefinition = coordinateDefinition;
	}
    public double getTemperature() {
		return temperature;
	}

	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}
	

	public double getLatticeEnergy() {
		return latticeEnergy;
	}

	public void setLatticeEnergy(double latticeEnergy) {
		this.latticeEnergy = latticeEnergy;
	}

}
