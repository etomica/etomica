package etomica;
import etomica.units.Dimension;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).
 *
 * @author David Kofke
 */
public class MeterRDF extends MeterFunction implements EtomicaElement {
	
	/**
	 * Creates meter with default to compute pair correlation for all
	 * leaf atoms in a phase.
	 * @param parent
	 */
    public MeterRDF(SimulationElement parent) {
	    super(parent, new DataSourceUniform());
	    xDataSourceUniform = (DataSourceUniform)xDataSource;
	    xDataSourceUniform.setDimension(Dimension.LENGTH);
	    xDataSourceUniform.setLabel("r");
	    setLabel("rdf");
	    iterator = new ApiLeafAtoms();
	    cPair = parent.space.makeCoordinatePair();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabulates radial distribution function");
        return info;
    }
    
    /**
     * Mutator method for the iterator that generates the atom pairs used to tabulate
     * the radial distribution function.  By setting this iterator the meter can be
     * configured to compute pair distribution for any set of atom pairs.  At construction
     * the default is an instance of ApiLeafAtoms, which generates pairs from all leaf
     * atoms in the phase.
     * @param iter
     */
    public void setIterator(AtomsetIteratorPhaseDependent iter) {iterator = iter;}
    
    /**
     * Accessor method for the iterator that generates the atom pairs used to
     * tabulate the radial distribution function.
     * @return
     */
    public AtomsetIteratorPhaseDependent getIterator() {return iterator;}
    
    /**
     * Returns Dimension.NULL, indicating that the measured quantity is dimensionless.
     */
    public Dimension getDimension() {return Dimension.NULL;}

	/**
	 * Computes RDF for the current configuration of the given phase.
	 */
	public double[] getDataAsArray(Phase phase) {
        cPair.setBoundary(phase.boundary());
	    for(int i=0; i<nDataPerPhase; i++) {phaseData[i] = 0.0;}  //zero histogram
	    double xMax = xDataSourceUniform.getXMax();
	    double xMaxSquared = xMax*xMax;
	    int count = 0;
		iterator.setPhase(phase);
	    iterator.reset();
	    while(iterator.hasNext()) {                 //iterate over all pairs
	    	Atom[] pair = iterator.next();
	    	cPair.reset(pair[0].coord, pair[1].coord);
	    	double r2 = cPair.r2();       //compute pair separation
	        if(r2 < xMaxSquared) {
	            int index = xDataSourceUniform.getIndex(Math.sqrt(r2));  //determine histogram index
	            phaseData[index]++;                        //add once for each atom
	        }
	        //TODO consider if this (count) is being used correctly
	        count++;
	    }
//	    int n = phase.atomCount();             //compute normalization: divide by
	    double norm = count/phase.volume();    //n, and density*(volume of shell)
	    double[] x = xDataSourceUniform.getData();
	    Space space = phase.simulation().space();
	    double dx2 = 0.5*(xMax - xDataSourceUniform.getXMin())/(double)nDataPerPhase;
	    for(int i=0; i<nDataPerPhase; i++) {
	        double vShell = space.sphereVolume(x[i]+dx2)-space.sphereVolume(x[i]-dx2);
	    	phaseData[i] /= (norm*vShell);
	    }
	    return phaseData;
	}
	
    private AtomsetIteratorPhaseDependent iterator;
    private final Space.CoordinatePair cPair;
    private final DataSourceUniform xDataSourceUniform;//local copy of xDataSource to avoid repeated casts
	
	/**
	 * main method to demonstrate and test use of class.
	 */
/*	 public static void main(String[] args) {
	    
	    etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
	    Simulation.instance = sim;
	    
	    MeterRDF meter = new MeterRDF(sim);
	    etomica.graphics.DisplayPlot plot = new etomica.graphics.DisplayPlot(sim);
	    etomica.graphics.DisplayTableFunction table = new etomica.graphics.DisplayTableFunction(sim);
	    ApiIntragroupAA iterator = new ApiIntragroupAA(sim);
	    sim.elementCoordinator.go();
	    SpeciesAgent agent = sim.phase.getAgent(sim.species);
	    meter.setIterator(iterator);
		iterator.setBasis(agent, agent);
	    
	    etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
	 }//end of main
	*/   
}
    