package etomica.data.meter;
import etomica.AtomPair;
import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.EtomicaInfo;
import etomica.Meter;
import etomica.Phase;
import etomica.Space;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;
import etomica.data.DataSourceUniform;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.CoordinatePair;
import etomica.units.Dimension;
import etomica.utility.NameMaker;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).
 *
 * @author David Kofke
 */
public class MeterRDF implements DataSource, Meter {
	
	/**
	 * Creates meter with default to compute pair correlation for all
	 * leaf atoms in a phase.
	 * @param parent
	 */
    public MeterRDF(Space space) {
	    this.space = space;
        data = new DataFunction(new DataInfo("RDF",Dimension.FRACTION),new DataInfo("r",Dimension.LENGTH));
        xDataSource = new DataSourceUniform();
        xDataSource.setTypeMax(DataSourceUniform.HALF_STEP);
        xDataSource.setTypeMin(DataSourceUniform.HALF_STEP);
        DataDoubleArray xData = (DataDoubleArray)xDataSource.getData();
        data.setLength(xData.getLength());
        data.getTData().E(xData);
	    iterator = new ApiLeafAtoms();
	    cPair = space.makeCoordinatePair();
        setName(NameMaker.makeName(this.getClass()));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabulates radial distribution function");
        return info;
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
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
	public Data getData() {
        cPair.setNearestImageTransformer(phase.boundary());
        if (data.getLength() != xDataSource.getNValues()) {
            data.setLength(xDataSource.getNValues());
        }
        final double[] y = data.getData();
	    for(int i=0; i<y.length; i++) {y[i] = 0.0;}  //zero histogram
	    double xMax = xDataSource.getXMax();
	    double xMaxSquared = xMax*xMax;
	    int count = 0;
		iterator.setPhase(phase);
	    iterator.reset();
	    while(iterator.hasNext()) {                 //iterate over all pairs
	    	AtomPair pair = (AtomPair)iterator.next();
	    	cPair.reset(pair.atom0.coord, pair.atom1.coord);
	    	double r2 = cPair.r2();       //compute pair separation
	        if(r2 < xMaxSquared) {
	            int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
	            y[index]++;                        //add once for each atom
	        }
	        //TODO consider if this (count) is being used correctly
	        count++;
	    }
//	    int n = phase.atomCount();             //compute normalization: divide by
	    double norm = count/phase.volume();    //n, and density*(volume of shell)
	    double[] x = ((DataDoubleArray)xDataSource.getData()).getData();
	    double dx2 = 0.5*(xMax - xDataSource.getXMin())/data.getLength();
	    for(int i=0;i<data.getLength(); i++) {
	        double vShell = space.sphereVolume(x[i]+dx2)-space.sphereVolume(x[i]-dx2);
	    	y[i] /= (norm*vShell);
	    }
	    return data;
	}
	
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private Phase phase;
    private final Space space;
    private final DataFunction data;
    private AtomsetIteratorPhaseDependent iterator;
    private final CoordinatePair cPair;
    private final DataSourceUniform xDataSource;
    private String name;
	
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
    