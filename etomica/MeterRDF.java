package etomica;
import etomica.units.Dimension;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF)
 * Should work ok in ensembles with fluctuating volume and particle numbers
 * Not suited for tabulation of RDFs of mixtures or heteroatomic molecules.
 *
 * @author David Kofke
 */
public class MeterRDF extends MeterFunction implements EtomicaElement {
    private AtomPairIterator iterator;
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private double[] vShell;
    protected double delr;
    
    public MeterRDF(SimulationElement parent) {
	    super(parent, new DataSourceUniform());
	    ((DataSourceUniform)xDataSource).setDimension(Dimension.LENGTH);
	    ((DataSourceUniform)xDataSource).setLabel("r");
	    setLabel("rdf");
	    iterator = AtomPairIterator.NULL;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabulates radial distribution function");
        return info;
    }
    
    public void setIterator(AtomPairIterator iter) {iterator = iter;}
    public AtomPairIterator getIterator() {return iterator;}
    
    public IteratorDirective getDirective() {return iteratorDirective;}

    public Dimension getDimension() {return Dimension.NULL;}

	/**
	 * Computes RDF for the current configuration
	 *    For future development: It may be possible to extend to particular atom pairs by changing iterator and using a different normalization
	 */
	public double[] getDataAsArray(Phase phase) {
		//TODO set up the iterator and phase
		iterator.reset(iteratorDirective);          //prepare iterator of atom pairs
	    for(int i=0; i<nDataPerPhase; i++) {phaseData[i] = 0.0;}  //zero histogram
	    double xMax = ((DataSourceUniform)xDataSource).getXMax();
	    double xMaxSquared = xMax*xMax;
	    while(iterator.hasNext()) {                 //iterate over all pairs
	        double r2 = iterator.next().r2();       //compute pair separation
	        if(r2 < xMaxSquared) {
	            int index = ((DataSourceUniform)xDataSource).getIndex(Math.sqrt(r2));  //determine histogram index
	            phaseData[index]++;                        //add once for each atom
	        }
	    }
	    int n = phase.atomCount();                  //compute normalization: divide by
	    double norm = (n*(n-1)/2)/phase.volume();           //n, and density*(volume of shell)
	    double[] x = xDataSource.getData();
	    Space space = phase.simulation().space();
	    double xMin = ((DataSourceUniform)xDataSource).getXMin();
	    double dx2 = 0.5*(xMax - xMin)/(double)nDataPerPhase;
	    for(int i=0; i<nDataPerPhase; i++) {
	        double vShell = space.sphereVolume(x[i]+dx2)-space.sphereVolume(x[i]-dx2);
	    	phaseData[i] /= (norm*vShell);
	    }
	    return phaseData;
	}
	
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
    