package etomica.threaded;

import etomica.potential.PotentialCalculation;

public interface IPotentialCalculationThreaded {

	public PotentialCalculation[] getPotentialCalculations();

	public void writeData();
}