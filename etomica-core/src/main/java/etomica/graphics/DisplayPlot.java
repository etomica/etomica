/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.data.*;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.Unit;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.units.systems.UnitSystem;

import javax.swing.*;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;

/**
 * Class for creating a plot of simulation data.  Data for plot is taken
 * from a DataSet instance, which can be obtained via the getDataSet
 * method.  Data sets are added to the plot by piping the data to a new DataSink
 * obtained from the makeDataSink method in DataSet.  The incoming data must be
 * either a DataFunction or a DataGroup of DataFunctions.  A DataGroup of
 * DataFunctions will be plotted as separate sets of data in the plot.
 *
 * @author David Kofke
 */

public class DisplayPlot extends Display implements DataSetListener {
    
    /**
     * Creates a plot with a new, empty, DataSinkTable.
     */
    public DisplayPlot() {
        this(new DataSet());
    }
    
    /**
     * Creates a plot using data from the given table.
     */
    public DisplayPlot(final DataSet dataSet) {
        super();
        this.dataSet = dataSet;
        dataSet.addDataListener(this);
        plot = new EtomicaPlot(this);
        panel = new JPanel();
        panel.add(plot);
        units = new Unit[0];
        labelList = new LinkedList<DataTagBag>();
        unitList = new LinkedList<DataTagBag>();
        drawLineList = new LinkedList<DataTagBag>();

        final JPopupMenu popupMenu = new JPopupMenu();
        JMenuItem resetMenuItem = new JMenuItem("reset");
        resetMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                plot.resetAxes();
                doUpdate();
            }
        });
        popupMenu.add(resetMenuItem);
        final JMenu rawMenu = new JMenu("raw data");
        popupMenu.add(rawMenu);
        PopupMenuListener popupMenuListener = new PopupMenuListener(){
            public void popupMenuCanceled(PopupMenuEvent e) {}

            public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {}

            public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
                rawMenu.removeAll();
                int nSource = dataSet.getDataCount();
                // populate the submenu with one element for each data set.
                for(int k=0; k<nSource; k++) {
                    String dataLabel = dataSet.getDataInfo(k).getLabel();
                    DataTagBag tagLabel = DataTagBag.getDataTagBag(labelList, dataSet.getDataInfo(k).getTags());
                    if (tagLabel != null) {
                        dataLabel = (String)tagLabel.object;
                    }
                    JMenuItem rawDataMenuItem = new JMenuItem(dataLabel);
                    rawDataMenuItem.addActionListener(new RawDataWindowOpener(k, DisplayPlot.this));
                    rawMenu.add(rawDataMenuItem);
                }
            }
        };
        popupMenu.addPopupMenuListener(popupMenuListener);
        final JMenu writeMenu = new JMenu("write data");
        popupMenu.add(writeMenu);
        popupMenuListener = new PopupMenuListener(){
            public void popupMenuCanceled(PopupMenuEvent e) {}

            public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {}

            public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
                writeMenu.removeAll();
                int nSource = dataSet.getDataCount();
                // populate the submenu with one element for each data set.
                for(int k=0; k<nSource; k++) {
                    String dataLabel = dataSet.getDataInfo(k).getLabel();
                    DataTagBag tagLabel = DataTagBag.getDataTagBag(labelList, dataSet.getDataInfo(k).getTags());
                    if (tagLabel != null) {
                        dataLabel = (String)tagLabel.object;
                    }
                    JMenuItem writeDataMenuItem = new JMenuItem(dataLabel);
                    writeDataMenuItem.addActionListener(new DataWriter(k, DisplayPlot.this));
                    writeMenu.add(writeDataMenuItem);
                }
            }
        };
        popupMenu.addPopupMenuListener(popupMenuListener);

        JMenu logMenuItem = new JMenu("toggle log");
        JMenuItem xlogMenuItem = new JMenuItem("x");
        xlogMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                boolean log = plot.getXLog();
                plot.setXLog(!log);
                doUpdate();
            }
        });
        logMenuItem.add(xlogMenuItem);
        JMenuItem ylogMenuItem = new JMenuItem("y");
        ylogMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                boolean log = plot.getYLog();
                plot.setYLog(!log);
                doUpdate();
            }
        });
        logMenuItem.add(ylogMenuItem);
        popupMenu.add(logMenuItem);

        plot.addMouseListener(new PopupListener(popupMenu));
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
    public void dataChanged(DataSet changedDataSet) {
        doUpdate();
    }
    
    /**
     * Updates the units array for the new column, using the default units.
     */
    public void dataCountChanged(DataSet changedDataSet) {
        // changedDataSet is our dataSet.  or at least, it better be.
        int oldColumnCount = units.length;
        int newDataCount = dataSet.getDataCount();
        if(newDataCount > oldColumnCount) {
            units = (Unit[]) Arrays.copyOf(units, newDataCount);
            Dimension xDimension = null;
            if (units.length > 0) {
                IDataInfo xDataInfo = ((DataInfoFunction)dataSet.getDataInfo(0)).getXDataSource().getIndependentDataInfo(0);
                xDimension = xDataInfo.getDimension();
                if (xUnit == null) {
                    xUnit = xDimension.getUnit(UnitSystem.SIM);
                }
                if (xLabel == null) {
                    xLabel = xDataInfo.getLabel();
                }
                plot.setXLabel(xLabel);
            }
            for(int i=oldColumnCount; i<units.length; i++) {
            	Unit dataUnit = null;

            	if(dataSet.getDataInfo(i) instanceof DataInfoFunction) {

	                DataInfoFunction dataInfo = (DataInfoFunction)dataSet.getDataInfo(i);
	                if (dataInfo.getXDataSource().getIndependentArrayDimension() != 1 || !dataInfo.getXDataSource().getIndependentDataInfo(0).getDimension().equals(xDimension)) {
	                    throw new IllegalArgumentException("All data functions must have the same X dimension");
	                }
	
	                dataUnit = (defaultUnit != null) ? defaultUnit : dataInfo.getDimension().getUnit(UnitSystem.SIM);
            	}
                DataTagBag tagUnit = DataTagBag.getDataTagBag(unitList, dataSet.getDataInfo(i).getTags());
                if (tagUnit != null) {
                    dataUnit = (Unit)tagUnit.object;
                }

                units[i] = dataUnit;
            }
        } else {
            units = Arrays.copyOf(units, newDataCount);

            // whack them all!
            for(int i=0; i<units.length; i++) {
                DataInfoFunction dataInfo = (DataInfoFunction)dataSet.getDataInfo(i);

                Unit dataUnit = (defaultUnit != null) ? defaultUnit : dataInfo.getDimension().getUnit(UnitSystem.SIM);
                
                // if the user specified a unit for this data specifically, use it.
                DataTagBag tagUnit = DataTagBag.getDataTagBag(unitList, dataSet.getDataInfo(i).getTags());
                if (tagUnit != null) {
                    dataUnit = (Unit)tagUnit.object;
                }
                units[i] = dataUnit;
            }
        }
        setDoLegend(doLegend);
    }

    /**
     * Redraws the plot.
     */
    public void doUpdate() {
        int nSource = dataSet.getDataCount();
        if (!plot.isShowing() && doClear) {
            // notify the plot we have updated data
            plot.notifyNeedUpdate();
            return;
        }
        // if we don't clear the plot, we'll draw now even if the plot isn't
        // showing, since the data won't be remembered otherwise
        if (doClear) {
            plot.clear(false);
        }
        boolean xLog = plot.getXLog();
        boolean yLog = plot.getYLog();
        for(int k=0; k<nSource; k++) {
        	double[] xValues = null;
            double[] data = null;
        	if(dataSet.getDataInfo(k) instanceof DataInfoFunction) {
                xValues = ((DataInfoFunction)dataSet.getDataInfo(k)).getXDataSource().getIndependentData(0).getData();
                data = ((DataFunction)dataSet.getData(k)).getData();

                DataTagBag tagDrawLine = DataTagBag.getDataTagBag(drawLineList, dataSet.getDataInfo(k).getTags());
                if (tagDrawLine != null) {
                    // don't (or perhaps do) draw lines for this data set
                    plot.setConnected((Boolean) tagDrawLine.object, k);
                }
                boolean drawLine = false;
                for(int i=0; i<data.length && i<xValues.length; i++) {
                    double y = units[k].fromSim(data[i]);
                    if (!Double.isNaN(y)) {
                        if ((yLog && y<=0) || (xLog && xValues[i]<=0)) {
                            if (!logWarn) {
                                System.err.println("non-positive values dropped from log scale");
                                logWarn = true;
                            }
                            continue;
                        }
                        plot.addPoint(k, xUnit.fromSim(xValues[i]), y, drawLine);
                        drawLine = true;
                    }
                    else {
                        drawLine = false;
                    }
                }
        	}
        }
        if (plot.isShowing()) {
            plot.repaint();
        }
    }

    /**
     * Mutator for flag determining if a legend is to be shown.
     * Default is true.
     */
    public void setDoLegend(boolean b) {
        // we can get called when there's just new data, so go through the 
        // whole song and dance even if we're not changing the flag.
        doLegend = b;
        plot.clearLegends();
        if (!b) {
            return;
        }
        for(int i=0; i<dataSet.getDataCount(); i++) {
            String dataLabel = dataSet.getDataInfo(i).getLabel();
            DataTagBag tagLabel = DataTagBag.getDataTagBag(labelList, dataSet.getDataInfo(i).getTags());
            if (tagLabel != null) {
                dataLabel = (String)tagLabel.object;
            }
            plot.addLegend(i, dataLabel);
        }
    }
    
    protected IData getDataFromSet(DataTag[] tags) {
        for(int i=0; i<dataSet.getDataCount(); i++) {
            IData thisData = dataSet.getData(i);
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
        DataTagBag bag = DataTagBag.getDataTagBagExact(labelList, dataTags);
        if (bag != null) {
            labelList.remove(bag);
        }
        labelList.add(new DataTagBag(dataTags, label));
        if (doLegend) setDoLegend(true);
    }

    public void setDoClear(boolean newDoClear) {
        doClear = newDoClear;
    }

    public boolean isDoClear() {
        return doClear;
    }

    /**
     * Sets the data set corresponding to the given tags to be have the data
     * points connected with lines or not.
     */
    public void setDoDrawLines(DataTag[] dataTags, boolean doDrawLines) {
        DataTagBag bag = DataTagBag.getDataTagBagExact(drawLineList, dataTags);
        if (bag != null) {
            drawLineList.remove(bag);
        }
        
        drawLineList.add(new DataTagBag(dataTags, doDrawLines));
        // if a data set's points aren't connected, we need to draw points
        if (!doDrawLines) {
            plot.setMarksStyle("points");
        }
    }
    
    /**
     * Accessor method to plot class so that its properties can be edited.
     */
    public final EtomicaPlot getPlot() {return plot;}
        
    /**
     * Overrides superclass method to return a JPanel holding the plot.
     */
    public java.awt.Component graphic() {
        return panel;
    }


    /**
     * Extend superclass method to update label with change of unit.
     */
    public void setXUnit(Unit u) {
        xUnit = u;
        if (dataSet.getDataCount() > 0) {
            DataInfoFunction dataInfoFunction = (DataInfoFunction)dataSet.getDataInfo(0);
            plot.setXLabel(dataInfoFunction.getXDataSource().getIndependentDataInfo(0).getLabel() + " ("+xUnit.symbol()+")");
        }
    }
    
    public void setXLabel(String newXLabel) {
        xLabel = newXLabel;
        plot.setXLabel(xLabel);
    }

    public void setUnit(DataTag[] dataTags, Unit newUnit) {
        DataTagBag bag = DataTagBag.getDataTagBagExact(unitList, dataTags);
        if (bag != null) {
            unitList.remove(bag);
        }

        unitList.add(new DataTagBag(dataTags, newUnit));

        // now go through and look for a current Data with these tags
        for(int i=0; i<units.length; i++) {
            // if the user specified a unit for this data specifically, use it.
            DataTagBag tagUnit = DataTagBag.getDataTagBag(unitList, dataSet.getDataInfo(i).getTags());
            if (tagUnit != null) {
                units[i]= (Unit)tagUnit.object;
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
    
    protected final DataSet dataSet;
    protected final EtomicaPlot plot;
    protected final LinkedList<DataTagBag> labelList;
    protected final LinkedList<DataTagBag> unitList;
    protected final LinkedList<DataTagBag> drawLineList;
    protected final javax.swing.JPanel panel;
    private boolean doLegend = true;
    protected boolean doClear = true;
    private Unit xUnit = Null.UNIT;
    protected String xLabel;
    private Unit[] units;
    private Unit defaultUnit;
    protected boolean logWarn;

    /**
     * Opens a window with the raw data from a plot.  The data is re-retrieved
     * from the DataSet, so it won't necessarily be consistent with the plot.
     *
     * @author Andrew Schultz
     */
    public static final class RawDataWindowOpener implements ActionListener {
        private final int kk;
        private final DisplayPlot displayPlot;

        public RawDataWindowOpener(int kk, DisplayPlot displayPlot) {
            this.kk = kk;
            this.displayPlot = displayPlot;
        }

        public void actionPerformed(ActionEvent e) {
            JFrame f = new JFrame();
            TextArea textArea = new TextArea();
            textArea.setEditable(false);
            textArea.setBackground(Color.white);
            textArea.setForeground(Color.black);
            DataSet dataSet = displayPlot.getDataSet();
            if(dataSet.getDataInfo(kk) instanceof DataInfoFunction) {
                String dataLabel = dataSet.getDataInfo(kk).getLabel();
                DataTagBag tagLabel = DataTagBag.getDataTagBag(displayPlot.labelList, dataSet.getDataInfo(kk).getTags());
                if (tagLabel != null) {
                    dataLabel = (String)tagLabel.object;
                }
                textArea.append(displayPlot.getPlot().getXLabel()+"\t"+dataLabel+"\n");
                double[] xValues = ((DataInfoFunction)dataSet.getDataInfo(kk)).getXDataSource().getIndependentData(0).getData();
                double[] data = ((DataFunction)dataSet.getData(kk)).getData();
                for(int i=0; i<data.length; i++) {
                    double y = displayPlot.units[kk].fromSim(data[i]);
                    if (!Double.isNaN(y)) {
                        textArea.append(displayPlot.xUnit.fromSim(xValues[i])+"\t"+y+"\n");
                    }
                }
            }
            f.add(textArea);
            f.pack();
            f.setSize(400,600);
            f.setVisible(true);
        }
    }

    /**
     * Writes a file with the raw data from a plot.  The data is re-retrieved
     * from the DataSet, so it won't necessarily be consistent with the plot.
     *
     * @author Andrew Schultz
     */
    public static final class DataWriter implements ActionListener {
        private final int kk;
        private final DisplayPlot displayPlot;

        public DataWriter(int kk, DisplayPlot displayPlot) {
            this.kk = kk;
            this.displayPlot = displayPlot;
        }

        public void actionPerformed(ActionEvent e) {
            try {
                FileWriter fw = new FileWriter("data"+kk+".dat");
                DataSet dataSet = displayPlot.getDataSet();
                String dataLabel = dataSet.getDataInfo(kk).getLabel();
                DataTagBag tagLabel = DataTagBag.getDataTagBag(displayPlot.labelList, dataSet.getDataInfo(kk).getTags());
                if (tagLabel != null) {
                    dataLabel = (String)tagLabel.object;
                }
                fw.write(displayPlot.getPlot().getXLabel()+"\t"+dataLabel+"\n");
                double[] xValues = ((DataInfoFunction)dataSet.getDataInfo(kk)).getXDataSource().getIndependentData(0).getData();
                double[] data = ((DataFunction)dataSet.getData(kk)).getData();
                for(int i=0; i<data.length; i++) {
                    double y = displayPlot.units[kk].fromSim(data[i]);
                    if (!Double.isNaN(y)) {
                        fw.write(displayPlot.xUnit.fromSim(xValues[i])+"\t"+y+"\n");
                    }
                }
                fw.close();
            }
            catch (IOException ex) {
                // eat it
            }
        }
    }

    public static class PopupListener extends MouseAdapter {
        
        public PopupListener(JPopupMenu popupMenu) {
            this.popupMenu = popupMenu;
        }
        
        public void mousePressed(MouseEvent e) {
            maybeShowPopup(e);
        }

        public void mouseReleased(MouseEvent e) {
            maybeShowPopup(e);
        }

        private void maybeShowPopup(MouseEvent e) {
            if (e.isPopupTrigger()) {
                popupMenu.show(e.getComponent(),
                           e.getX(), e.getY());
            }
        }
        
        protected final JPopupMenu popupMenu;
    }

}
