package etomica.virial;

import etomica.AtomSet;
import etomica.Default;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential0;
import etomica.PotentialCalculation;
import etomica.PotentialCalculationEnergySum;
import etomica.SimulationElement;

/**
 * @author David Kofke
 *
 * Pair potential given according to the Mayer bonds in a cluster integral.
 */
public class P2Cluster extends Potential0 {

	protected Cluster cluster;
	protected PairSet pairs;
	protected double g = 0;
	/**
	 * Constructor for P2Cluster.
	 * @param parent
	 */
	public P2Cluster(SimulationElement parent, PairSet pairs) {
		super(parent);
		this.pairs = pairs;
		setTemperature(Default.TEMPERATURE);
	}
	
	public void setCluster(Cluster cluster) {
		this.cluster = cluster;
	}
	
	public double pi() {return cluster.value(pairs, beta);}

	public void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
	   if(!enabled) return;
	   switch(id.atomCount()) {
	   		case 0: g = cluster.value(pairs, beta); break;
	   		case 1: g = cluster.value(id.atom1(), pairs, beta); break;
	   		default: throw new RuntimeException();
	   }
	   ((PotentialCalculationEnergySum)pc).set(this).actionPerformed((Phase)null);
	}//end of calculate

	public double energy(Phase phase) {
		return (g==1.0) ? 0.0 : -temperature*Math.log(g);
	}

	protected double temperature, beta;
	/**
	 * Returns the temperature.
	 * @return double
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets the temperature.
	 * @param temperature The temperature to set
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
		beta = 1.0/temperature;
	}

}
