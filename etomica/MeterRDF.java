package etomica;
import java.util.Observer;
import java.util.Observable;
import etomica.units.Dimension;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF)
 * Should work ok in ensembles with fluctuating volume and particle numbers
 * Not suited for tabulation of RDFs of mixtures or heteroatomic molecules.
 *
 * @author David Kofke
 */
public class MeterRDF extends MeterFunction implements EtomicaElement {
    
    private AtomPair.Iterator iterator;
    private double[] vShell;
    private double delr;
    private DataSourceUniform xSource;
    
    public MeterRDF() {
        this(Simulation.instance);
    }
    public MeterRDF(Simulation sim) {
	    super(sim);
	    xSource = new DataSourceUniform();
	    xSource.setLabel("r");
	    setXLabel("r");
	    setLabel("rdf");
	    setActive(true);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabulates radial distribution function");
        return info;
    }

    public Dimension getDimension() {return Dimension.NULL;}
    public Dimension getXDimension() {return Dimension.LENGTH;}
    
    /**
     * This meter needs iterators to do its measurements, so this method overrides the no-op method of AbstractMeter 
     * It obtains the necessary iterators from the phase.
     */
	protected void setPhaseIteratorFactory(IteratorFactory factory) {
	    iterator = factory.makeAtomPairIteratorAll();
	}
	/**
	 * Returns true to flag that this meter uses an iterator.
	 */
	protected boolean usesPhaseIteratorFactory() {return true;}

	protected void setPhaseBoundary(Space.Boundary b) {
//	    if(phase == null || phase.boundary() == null) return;
	    setX(0.0, 0.5*b.dimensions().component(0),nPoints);
	}
	/**
	 * Returns true to flag that this meter uses the phase boundary.
	 */
	protected boolean usesPhaseBoundary() {return true;}

    /**
     * Sets phase, evaluates abscissa values to maximum of half system edge-length, and constructs atom iterator
     */
	public void setPhase(Phase p) {
	    super.setPhase(p);
	    setX(xMin, xMax, nPoints);
	}
	
	/**
	 * Computes RDF for the current configuration
	 *    For future development: It may be possible to extend to particular atom pairs by changing iterator and using a different normalization
	 */
	public double[] currentValue() {
	    iterator.reset();                           //prepare iterator of atom pairs
	    for(int i=0; i<nPoints; i++) {y[i] = 0.0;}  //zero histogram
	    while(iterator.hasNext()) {                 //iterate over all pairs
	        double r = Math.sqrt(iterator.next().r2()); //compute pair separation
	        if(r < xMax) {
	            int index = (int)(r/delr);          //determine histogram index
	            y[index]+=2;                        //add once for each atom
	        }
	    }
	    int n = phase.atomCount();                  //compute normalization: divide by
	    double norm = n*n/phase.volume();           //n, and density*(volume of shell)
	    for(int i=0; i<nPoints; i++) {y[i] /= (norm*vShell[i]);}
	    return y;
	}
	
	/**
	 * Sets the radial values at which RDF is tabulated
	 */
    public void setX(double min, double max, int n) {
	    super.setX(min, max, n);
	    if(phase == null) return;
	    //Compute normalization constants for RDF, including shell volumes for ideal-gas particle numbers
	    vShell = new double[n];         
	    Space space = phase.parentSimulation().space();
	    double dx2 = 0.5*(xMax - xMin)/(double)nPoints;
	    for(int i=0; i<nPoints-1; i++) {
	        vShell[i] = space.sphereVolume(x[i]+dx2)-space.sphereVolume(x[i]-dx2);
	    }
	    delr = xMax/(double)(nPoints-1);
	}
	
	/**
	 * main method to demonstrate and test use of class.
	 */
	 public static void main(String[] args) {
	    
	    javax.swing.JFrame frame = new javax.swing.JFrame();
        frame.setSize(600,350);
	    
	    etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
	    Simulation.instance = sim;
	    
	    MeterRDF meter = new MeterRDF(sim);
	    meter.setActive(true);
	    DisplayPlot plot = new DisplayPlot(sim);
	    DisplayTableFunction table = new DisplayTableFunction(sim);
	    
	    sim.elementCoordinator.go();
	    
	    frame.getContentPane().add(sim.panel());
	    frame.pack();
	    frame.show();
        frame.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
	 }//end of main
	    
}
    