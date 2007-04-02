package etomica.threaded;

import etomica.potential.PotentialCalculation;

public interface PotentialCalculationThreaded {

	public PotentialCalculation[] getPotentialCalculations();

	public void writeData();
}