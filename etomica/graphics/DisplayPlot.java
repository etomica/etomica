package etomica.graphics;
import java.io.Serializable;
import java.util.Iterator;
import java.util.LinkedList;

import ptolemy.plot.Plot;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.data.DataSet;
import etomica.data.DataSetListener;
import etomica.data.DataTag;
import etomica.data.DataSet.DataCasterJudge;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.units.Dimension;
import etomica.units.Null;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;
import etomica.util.Arrays;

/**
 * Class for creating a plot of simulation data.  Data for plot is taken
 * from a DataSinkTable instance, which can be obtained via the getDataTable
 * method.  Data sets are added to the plot by piping the data to a new DataSink
 * obtained from the makeColumn method of the data table.  
 *
 * @author David Kofke
 */

public class DisplayPlot extends Display implements DataSetListener, EtomicaElement {
    
    /**
     * Creates a plot with a new, empty, DataSinkTable.
     */
    public DisplayPlot() {
        this(new DataSet(new DataCasterJudgeFunction()));
    }
    
    /**
     * Creates a plot using data from the given table.
     */
    public DisplayPlot(DataSet dataSet) {
        super();
        this.dataSet = dataSet;
        if (!(dataSet.getDataCasterJudge() instanceof DataCasterJudgeFunction)) {
            // consider just calling the default constructor if you're hitting this exception
            throw new IllegalArgumentException("DisplayPlot a DataSet with a DataCasterJudgeFunction");
        }
        dataSet.addDataListener(this);
        plot = new Plot();
        panel.add(plot);
        setName("Data Plot");
        units = new Unit[0];
        labelList = new LinkedList();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("X-Y graphical plot of data");
        return info;
    }
    
    /**
     * Returns the DataSinkTable instance that specifies the plot data.
     * Data sets are added to the plot by piping them to a new data sink
     * given by the makeColumn method of the table.
     * @return
     */
    public DataSet getDataSet() {
        return dataSet;
    }
    
    /**
     * Causes the display of the plot to be updated.
     */
    public void dataChanged(DataSet dataSet) {
        doUpdate();
    }
    
    /**
     * Updates the units array for the new column, using the default units.
     */
    public void dataCountChanged(DataSet dataSet) {
        int oldColumnCount = units.length;
        int newDataCount = dataSet.getDataCount();
        if(newDataCount > oldColumnCount) {
            String xLabel = "";
            units = (Unit[])Arrays.resizeArray(units, newDataCount);
            Dimension xDimension = null;
            if (units.length > 0) {
                DataInfo xDataInfo = ((DataInfoFunction)dataSet.getDataInfo(0)).getIndependentInfo(0);
                xDimension = xDataInfo.getDimension();
                if (xUnit == null) {
                    xUnit = xDimension.getUnit(UnitSystem.SIM);
                }
                xLabel = xDataInfo.getLabel();
                plot.setXLabel(xLabel);
            }
            for(int i=oldColumnCount; i<units.length; i++) {
                DataFunction iData = (DataFunction)dataSet.getData(i);
                DataInfoFunction dataInfo = (DataInfoFunction)dataSet.getDataInfo(i);
                if (iData.getXDimension() != 1 || !dataInfo.getIndependentInfo(0).getDimension().equals(xDimension)) {
                    throw new IllegalArgumentException("All data functions must have the same X dimension");
                }
                units[i] = (defaultUnit != null) ? defaultUnit : dataInfo.getDimension().getUnit(UnitSystem.SIM);
            }
        } else {
            //TODO have DisplayTable adjust appropriately to removal of columns; this works only for newColumnCount = 0
            //may need to remember columns and reassign units to match
            units = (Unit[])Arrays.resizeArray(units, newDataCount);
        }
        setDoLegend(doLegend);
    }

    /**
     * Redraws the plot.
     */
    public void doUpdate() {
        if(!plot.isShowing()) return;
        int nSource = dataSet.getDataCount();
        plot.clear(false);
        for(int k=0; k<nSource; k++) {
            DataFunction dataFunction = (DataFunction)dataSet.getData(k);
            final double[] xValues = dataFunction.getXData(0).getData();
            final double[] data = dataFunction.getData();
            for(int i=0; i<data.length; i++) {
                double y = units[k].fromSim(data[i]);
                if (!Double.isNaN(y)) {
                    plot.addPoint(k, xUnit.fromSim(xValues[i]), y, true);
                }
            }
        }
        plot.repaint();
    }

    /**
     * Mutator for flag determining if a legend is to be shown.
     * Default is true.
     */
    public void setDoLegend(boolean b) {
        doLegend = b;
        for(int i=0; i<dataSet.getDataCount(); i++) {
            plot.removeLegend(i);
            if (b) {
                Iterator iterator = labelList.iterator();
                String dataLabel = dataSet.getDataInfo(i).getLabel();
                DataTag[] thisDataTags = dataSet.getDataInfo(i).getTags();
                while (iterator.hasNext()) {
                    DataTagLabel tagLabel = (DataTagLabel)iterator.next();
                    DataTag[] tags = tagLabel.tags;
                    boolean found = true;
                    for (int j=0; j<tags.length; j++) {
                        found = false;
                        for (int k=0; k<thisDataTags.length; k++) {
                            if (thisDataTags[k] == tags[j]) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            break;
                        }
                    }
                    if (found) {
                        dataLabel = tagLabel.label;
                    }
                }
                plot.addLegend(i, dataLabel);
            }
        }
    }
    
    protected Data getDataFromSet(DataTag[] tags) {
        for(int i=0; i<dataSet.getDataCount(); i++) {
            Data thisData = dataSet.getData(i);
            DataTag[] thisDataTags = dataSet.getDataInfo(i).getTags();
            boolean found = true;
            for (int j=0; j<tags.length; j++) {
                found = false;
                for (int k=0; k<thisDataTags.length; k++) {
                    if (thisDataTags[k] == tags[j]) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    break;
                }
            }
            if (found) {
                return thisData;
            }
        }
        return null;
    }
    
    /**
     * Accessor for flag determining if a legend is to be shown.
     * Default is true.
     */
    public boolean isDoLegend() {return doLegend;}
    
    public void setLegend(DataTag[] dataTags, String label) {
        labelList.add(new DataTagLabel(dataTags, label));
        
    }
    
    /**
     * Accessor method to plot class so that its properties can be edited.
     */
    public final Plot getPlot() {return plot;}
        
    /**
     * Overrides superclass method to return a JPanel holding the plot.
     */
    public java.awt.Component graphic(Object obj) {return panel;}


    /**
     * Extend superclass method to update label with change of unit.
     */
    public void setXUnit(Unit u) {
        xUnit = u;
        if (dataSet.getDataCount() > 0) {
            DataInfoFunction dataInfoFunction = (DataInfoFunction)dataSet.getDataInfo(0);
            plot.setXLabel(dataInfoFunction.getIndependentInfo(0).getLabel() + " ("+xUnit.symbol()+")");
        }
    }
    
    /**
     * Sets the display units of the given column to the given unit.
     */
    public void setUnit(Data data, Unit newUnit) {
        for(int i=0; i<dataSet.getDataCount(); i++) {
            if(dataSet.getData(i) == data) {
                units[i] = newUnit;
                break;
            }
        }
    }
    
    public void setUnit(Unit newUnit) {
    	defaultUnit = newUnit;
        for(int i=0; i<dataSet.getDataCount(); i++) {
            units[i] = newUnit;
        }
    }

    /**
     * Sets the drawn size of the plot
     * @param width new width, in pixels
     * @param height new height, in pixels
     */
    public void setSize(int width, int height) {
        plot.setSize(width, height);
        panel.setSize(width, height);
    }
    
    private final DataSet dataSet;
    private Plot plot;
    private javax.swing.JPanel panel = new javax.swing.JPanel();
    private boolean doLegend = true;
    private Unit xUnit = Null.UNIT;
    private Unit[] units;
    private Unit defaultUnit;
    private LinkedList labelList;
    
    protected static class DataTagLabel {
        public DataTagLabel(DataTag[] tags, String label) {
            this.tags = (DataTag[])tags.clone();
            this.label = label;
        }
        public final DataTag[] tags;
        public final String label;
    }
    
    protected static class DataCasterJudgeFunction implements DataCasterJudge, Serializable {

        public DataProcessor getDataCaster(DataInfo dataInfo) {
            if (dataInfo instanceof DataInfoFunction) {
                return null;
            } else if(dataInfo instanceof DataInfoGroup) {
                for (int i = 0; i<((DataInfoGroup)dataInfo).getNDataInfo(); i++) {
                    if (!(((DataInfoGroup)dataInfo).getSubDataInfo(i) instanceof DataInfoFunction)) {
                        throw new IllegalArgumentException("DisplayPlot can only handle homogeneous groups of DataFunctions");
                    }
                }
                return null;
            }
            throw new RuntimeException("DisplayPlot can only handle DataFunctions");
        }

    }

}