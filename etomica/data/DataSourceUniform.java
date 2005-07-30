package etomica.data;

import etomica.Constants;
import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.data.types.DataDoubleArray;
import etomica.units.Dimension;

/**
 * A DataSource object that provides a set of uniformly spaced values between
 * two limits. Number of points and range can be configured. Range can be set
 * inclusive or exclusive of end points (independently). A physical dimension
 * (e.g., length, time) can be associated with the data at construction (only).
 * <br>
 * Used principally to provide a set of "x" values for a function or plot.
 * 
 * @author David Kofke
 */
/*
 * History
 * ?? 		DAK				Created
 * 8/4/04	DAK,AJS,NRC	Add Index method
 */

public class DataSourceUniform implements DataSource, java.io.Serializable {
    
    /**
     * Default constructor. Chooses 100 points between 0 and 1, inclusive,
     * and sets dataInfo label to "Uniform" and dimension to Dimension.NULL 
     * (which cannot be changed after construction).
     */
    public DataSourceUniform() {
        this("Uniform", Dimension.NULL);
    }
    
    public DataSourceUniform(String label, Dimension dimension) {
        this(label, dimension, 100, 0.0, 1.0);
    }
    
    /**
     * Constructs a DataSourceUniform using the indicated number of values
     * and x range, and INCLUSIVE limits.
     */
    public DataSourceUniform(String label, Dimension dimension, int nValues, double xMin, double xMax) {
    	this(label, dimension, nValues, xMin, xMax, INCLUSIVE, INCLUSIVE);
    }
    
    /**
     * Constructs a DataSourceUniform using the indicated number of values,
     * x range, and limit types, with dataInfo indicating that the data have
     * the given dimension.  Dimension cannot be changed after construction.
     */
    public DataSourceUniform(String label, Dimension dimension, int nValues, double xMin, double xMax, 
                   LimitType typeMin, LimitType typeMax) {
        data = new DataDoubleArray(label, dimension, nValues);
        if(nValues < 2) nValues = 2;
        calculateX(nValues, xMin, xMax, typeMin, typeMax);
    }
    
    /**
     * Method that performs calculation of values.
     */
    private void calculateX(int nValues, double xMin, double xMax, 
                   LimitType typeMin, LimitType typeMax) {
        
        //determine number of intervals between left and right end points
        double intervalCount = (nValues - 1);

        
        /*
         * See inner class LimitType (below) for more information.
         */
        if(typeMin == HALF_STEP) intervalCount += 0.5;
        else if(typeMin == EXCLUSIVE) intervalCount += 1.0;
        
        if(typeMax == HALF_STEP) intervalCount += 0.5;
        else if(typeMax == EXCLUSIVE) intervalCount += 1.0;
        
        //determine spacing between values
        dx = (xMax - xMin)/intervalCount;
        
        //determine first value
        double x0 = xMin;
        if(typeMin == HALF_STEP) x0 += 0.5*dx;
        else if(typeMin == EXCLUSIVE) x0 += dx;
        
        //calculate values

        if (nValues != data.getLength()) {
            data = new DataDoubleArray(data.getDataInfo(), new int[]{nValues});
        }
        x = data.getData();
        for(int i=0; i<nValues; i++) x[i] = x0 + i*dx;
        
        //save parameters for later use
        this.xMin = xMin;
        this.xMax = xMax;
        this.typeMin = typeMin;
        this.typeMax = typeMax;
    }//end of calculateX
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    /**
     * Returns the index of the x value closest to the argument.
     * @throws an IllegalArgumentException if the argument is outside the defined range.
     * @param the double you are testing
     */
    public int getIndex(double value) {
    	// Test whether or not the desired value is in the defined range.
    	if(value < xMin || value > xMax) throw new IllegalArgumentException("Value outside of defined range in DataSourceUniform.getIndex. Value: "+value+" (Min, Max): ("+xMin+", "+xMax+")");
    	int i =(int)((value - x[0]) / dx + 0.5);
    	// Comparing to indices; integers are faster.
    	if(i < 0) return 0;
    	if(i > x.length-1) return x.length - 1;
    	return i;
    }
    
    /**
     * Sets the number of uniformly spaced values.
     * Minimum allowable value is 2 (smaller values cause return with no action).
     */
    public void setNValues(int nValues) {
        if(nValues < 2) return;
        calculateX(nValues, xMin, xMax, typeMin, typeMax);
    }
    
    /**
     * Accessor method for number of values.
     */
     public int getNValues() {
         return x.length;
     }
    
    /**
     * Mutator method for the left limit of the range.
     */
     public void setXMin(double xMin) {
        calculateX(x.length, xMin, xMax, typeMin, typeMax);
     }
     /**
      * Accessor method for the left limit of the range.
      */
      public double getXMax() {return xMax;}
      
    /**
     * Mutator method for the right limit of the range.
     */
     public void setXMax(double xMax) {
        calculateX(x.length, xMin, xMax, typeMin, typeMax);
     }
     /**
      * Accessor method for the right limit of the range.
      */
      public double getXMin() {return xMin;}
     
     /**
      * Mutator method for the type of the left limit.
      */
      public void setTypeMin(LimitType typeMin) {
        calculateX(x.length, xMin, xMax, typeMin, typeMax);
      }
      /**
       * Accessor method for the type of the left limit.
       */
      public LimitType getTypeMin() {return typeMin;}
        
     /**
      * Mutator method for the type of the right limit.
      */
      public void setTypeMax(LimitType typeMax) {
        calculateX(x.length, xMin, xMax, typeMin, typeMax);
      }
      
      /**
       * Accessor method for the type of the right limit.
       */
      public LimitType getTypeMax() {return typeMax;}
       
      /**
       * Returns the uniformly spaced values.
        */
      public Data getData() {return data;}
            
	/**
	 * Typed constant that indicates the way limits of the range are interpreted.
	 * Choices for the left and right limits may be made independently.<br>
	 * <ul>
	 *   <li>INCLUSIVE indicates that the end point of the series equals the limit.
	 *   <li>HALF_STEP indicates that the end point lies one half interval-step from the limit.
	 *   <li>EXCLUSIVE indicates that the end point lies one full interval-step from the limit.
	 * </ul>
	 * 
	 * For example, with n = 3 points and limits of 0.0 and 1.0, and with the left point (0.0)
	 * taken to have limitType = INCLUSIVE, the values are as follows for the different
	 * choices of the limitType of the right point (1.0):
	 * <ul>
	 *   <li>INCLUSIVE: 0.0, 0.5, 1.0
	 *   <li>HALF_STEP: 0.0, 0.4, 0.8
	 *   <li>EXCLUSIVE: 0.0, 0.333, 0.667
	 * </ul>
	 *
	 * To get values that lie at the centers of a set of bins (e.g., for histogramming),
	 * use HALF_STEP for both limits.
	 */
	public static class LimitType extends Constants.TypedConstant {
        public LimitType(String label) {super(label);}       
        public Constants.TypedConstant[] choices() {return CHOICES;}
        public static final LimitType[] CHOICES = 
            new LimitType[] {
                new LimitType("Inclusive"),
                new LimitType("Half step"),
                new LimitType("Exclusive"),};
    }//end of LimitType
    public static final LimitType INCLUSIVE = LimitType.CHOICES[0];
    public static final LimitType HALF_STEP = LimitType.CHOICES[1];
    public static final LimitType EXCLUSIVE = LimitType.CHOICES[2];
    
    private LimitType typeMin, typeMax;
    private double xMin, xMax;
    private double[] x;
    private DataDoubleArray data;
    private double dx;
    
    /**
     * Main method to demonstrate and check class.
     * Should output the values described in the comments for the LimitType class.
     */
    public static void main(String args[]) {
        
        DataSourceUniform source = new DataSourceUniform();
        
        int n = 3;
        source.setXMin(0.0);
        source.setXMax(1.0);
        source.setNValues(n);
        source.setTypeMin(DataSourceUniform.INCLUSIVE);
        source.setTypeMax(DataSourceUniform.INCLUSIVE);
        for(int i=0; i<n; i++) System.out.print(((DataDoubleArray)source.getData()).getData()[i] + "  ");
        System.out.println();
        source.setTypeMax(DataSourceUniform.HALF_STEP);
        for(int i=0; i<n; i++) System.out.print(((DataDoubleArray)source.getData()).getData()[i] + "  ");
        System.out.println();
        source.setTypeMax(DataSourceUniform.EXCLUSIVE);
        for(int i=0; i<n; i++) System.out.print(((DataDoubleArray)source.getData()).getData()[i] + "  ");
        System.out.println("");
        System.out.println("Done");
        
    }//end of main
    
 }//end of class
