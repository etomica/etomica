package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataSourceIndependent;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.phase.Phase;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Length;
import etomica.units.Null;
import etomica.util.NameMaker;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).
 *
 * @author David Kofke
 */
public class MeterRDF implements DataSource, Meter, DataSourceIndependent, java.io.Serializable {
	
	/**
	 * Creates meter with default to compute pair correlation for all
	 * leaf atoms in a phase.
	 * @param parent
	 */
    public MeterRDF(Space space) {
	    this.space = space;

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        
        rData = (DataDoubleArray)xDataSource.getData();
        data = new DataFunction(new int[] {rData.getLength()});
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);

	    iterator = new ApiLeafAtoms();
        setName(NameMaker.makeName(this.getClass()));
        dr = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabulates radial distribution function");
        return info;
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Mutator method for the iterator that generates the atom pairs used to tabulate
     * the radial distribution function.  By setting this iterator the meter can be
     * configured to compute pair distribution for any set of atom pairs.  At construction
     * the default is an instance of ApiLeafAtoms, which generates pairs from all leaf
     * atoms in the phase.
     */
    public void setIterator(AtomsetIteratorPhaseDependent iter) {
        iterator = iter;
    }
    
    /**
     * Accessor method for the iterator that generates the atom pairs used to
     * tabulate the radial distribution function.
     */
    public AtomsetIteratorPhaseDependent getIterator() {
        return iterator;
    }
    
	/**
	 * Computes RDF for the current configuration of the given phase.
	 */
	public Data getData() {
        boolean needUpdate = false;
        if (rData != xDataSource.getData()) {
            rData = (DataDoubleArray)xDataSource.getData();
            needUpdate = true;
        }
        if (data.getLength() != rData.getLength()) {
            needUpdate = true;
        }
        if (needUpdate) {
            data = new DataFunction(new int[] {rData.getLength()});
            dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);
            dataInfo.addTag(tag);
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
            dr.Ev1Mv2(((AtomLeaf)pair.atom1).getCoord().getPosition(),((AtomLeaf)pair.atom0).getCoord().getPosition());
            nearestImageTransformer.nearestImage(dr);
	    	double r2 = dr.squared();       //compute pair separation
	        if(r2 < xMaxSquared) {
	            int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
	            y[index]++;                        //add once for each atom
	        }
	        //TODO consider if this (count) is being used correctly
	        count++;
	    }
//	    int n = phase.atomCount();             //compute normalization: divide by
	    double norm = count/phase.volume();    //n, and density*(volume of shell)
	    double[] r = rData.getData();
	    double dx2 = 0.5*(xMax - xDataSource.getXMin())/r.length;
	    for(int i=0;i<r.length; i++) {
	        double vShell = space.sphereVolume(r[i]+dx2)-space.sphereVolume(r[i]-dx2);
	        y[i] /= (norm*vShell);
	    }
	    return data;
	}
    
    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }
	
    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray)xDataSource.getData();
    }
    
    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)xDataSource.getDataInfo();
    }
    
    public int getIndependentArrayDimension() {
        return 1;
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
        nearestImageTransformer = phase.getBoundary();
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    private static final long serialVersionUID = 1L;
    private Phase phase;
    private final Space space;
    protected DataFunction data;
    private DataInfo dataInfo;
    private DataDoubleArray rData;
    private AtomsetIteratorPhaseDependent iterator;
    private final Vector dr;
    private NearestImageTransformer nearestImageTransformer;
    private final DataSourceUniform xDataSource;
    private String name;
    protected final DataTag tag;
}
