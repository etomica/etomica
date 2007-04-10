package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.action.Action;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataSourceIndependent;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.units.Length;
import etomica.units.Null;
import etomica.util.NameMaker;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).  The
 * meter takes data via actionPerformed and returns the average RDF via
 * getData.
 *
 * @author David Kofke
 */
public class MeterRDF implements Action, DataSource, DataSourceIndependent, java.io.Serializable {
	
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
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);

	    iterator = new ApiLeafAtoms();
        setName(NameMaker.makeName(this.getClass()));
        dr = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public String getLabel() {
        throw new RuntimeException("You deserve this");
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabulates radial distribution function");
        return info;
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        rData = (DataDoubleArray)xDataSource.getData();
        xMax = xDataSource.getXMax();
        data = new DataFunction(new int[] {rData.getLength()});
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        callCount = 0;
    }
    
    /**
     * Takes the RDF for the current configuration of the given phase.
     */
    public void actionPerformed() {
        if (rData != xDataSource.getData() ||
            data.getLength() != rData.getLength() ||
            xDataSource.getXMax() != xMax) {
            reset();
        }
        
        double xMaxSquared = xMax*xMax;
        iterator.setPhase(phase);
        iterator.reset();
        while(iterator.hasNext()) {                 //iterate over all pairs
            AtomPair pair = (AtomPair)iterator.next();
            dr.Ev1Mv2(((AtomLeaf)pair.atom1).getPosition(),((AtomLeaf)pair.atom0).getPosition());
            nearestImageTransformer.nearestImage(dr);
            double r2 = dr.squared();       //compute pair separation
            if(r2 < xMaxSquared) {
                int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
                gSum[index]++;                        //add once for each atom
            }
        }
        callCount++;
    }
    
	/**
	 * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
	 */
	public Data getData() {
        if (rData != xDataSource.getData() ||
            data.getLength() != rData.getLength() ||
            xDataSource.getXMax() != xMax) {
            reset();
            //that zeroed everything.  just return the zeros.
            return data;
        }
        
        final double[] y = data.getData();
	    double norm = 0.5 * callCount * phase.getDensity()*(phase.getSpeciesMaster().getLeafList().size()-1);
	    double[] r = rData.getData();
	    double dx2 = 0.5*(xMax - xDataSource.getXMin())/r.length;
	    for(int i=0;i<r.length; i++) {
	        double vShell = space.sphereVolume(r[i]+dx2)-space.sphereVolume(r[i]-dx2);
	        y[i] = gSum[i] / (norm*vShell);
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
    protected Phase phase;
    protected final Space space;
    protected long[] gSum;
    protected DataFunction data;
    private IDataInfo dataInfo;
    protected DataDoubleArray rData;
    protected AtomsetIteratorPhaseDependent iterator;
    private final IVector dr;
    private NearestImageTransformer nearestImageTransformer;
    protected final DataSourceUniform xDataSource;
    protected double xMax;
    private String name;
    protected final DataTag tag;
    protected int callCount;
}
