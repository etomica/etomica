/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Component;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Formatter;
import java.util.LinkedList;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;

import etomica.data.*;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressureHard;
import etomica.data.types.DataDoubleArray;
import etomica.simulation.prototypes.HSMD2D;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Presents data in a tabular form.
 * 
 * @author David Kofke
 */
public class DisplayTable extends Display implements DataTableListener {

    public DisplayTable() {
        this(new DataSinkTable());
    }

    public DisplayTable(DataSinkTable dataTable) {
        this.dataTable = dataTable;
        stringBuilder = new StringBuilder();
        formatter = new Formatter(stringBuilder);

        units = new Unit[dataTable.getDataCount()];
        for (int i = 0; i < units.length; i++) {
            units[i] = dataTable.getDataInfo(i).getDimension().getUnit(UnitSystem.SIM);
        }

        dataTable.addDataListener(this);
        unitList = new LinkedList<DataTagBag>();
        columnHeaderList = new LinkedList<DataTagBag>();
        tableSource = new MyTable(); //inner class, defined below
        table = new JTable(tableSource);
        panel = new javax.swing.JPanel(new java.awt.FlowLayout());

        setLabel("Data");
        setPrecision(4);
        setTransposed(false);
        setShowingRowLabels(true);
        setShowingColumnHeaders(true);
        panel.add(new JScrollPane(table));
        InputEventHandler listener = new InputEventHandler();
        panel.addKeyListener(listener);
        panel.addMouseListener(listener);
        panel.setSize(100, 150);
        if (!fitToWindow)
            table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    }

    public DataSinkTable getDataTable() {
        return dataTable;
    }

    /**
     * Causes the display of the plot to be updated.
     */
    public void dataChanged(DataSet dummyTable) {
        tableSource.fireTableDataChanged();
        repaint();
    }

    /**
     * Updates the units array for any new or deleted columns, 
     * using the default units for new columns.
     */
    public void dataCountChanged(DataSet dummyTable) {
        recomputeUnits();
        recomputeColumnHeaders();
        tableSource.fireTableStructureChanged();
    }

    public void tableRowHeadersChanged(DataSinkTable dummyTable) {
        if (!autoRowLabels) return;
        for (int i=0; i<rowLabels.length; i++) {
            rowLabels[i] = dataTable.getRowHeader(i);
        }
    }

    /**
     * Part of the DataTableListener interface.  Updates the row headers.
     */
    public void tableRowCountChanged(DataSinkTable dummyTable) {
        if (autoRowLabels) {
            rowLabels = new String[dataTable.getRowCount()];
            tableRowHeadersChanged(dummyTable);
        }
        tableSource.fireTableStructureChanged();
    }

    /**
     * If true, each data set fills a row; if false, data are arranged to fill
     * columns.
     */
    public boolean isTransposed() {
        return transposed;
    }

    /**
     * If fillHorizontal is true, each data set fills a row; if false, data are
     * arranged to fill columns.
     */
    public void setTransposed(boolean transposed) {
        if (this.transposed != transposed) {
            this.transposed = transposed;
            tableSource.fireTableStructureChanged();
        }
    }

    /**
     * If true, columns will be squeezed so table fits to window; if false,
     * table will exceed window size, and to view full width must be placed in a
     * scroll pane. Default is true.
     */
    public void setFitToWindow(boolean b) {
        fitToWindow = b;
        if (fitToWindow)
            table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
        else
            table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    }

    /**
     * Accessor method for flag determining if table is squeezed to fit in
     * window.
     */
    public boolean isFitToWindow() {
        return fitToWindow;
    }

    /**
     * Mutator method for whether table should include a column of "x" values.
     */
    public void setShowingRowLabels(boolean b) {
        showingRowLabels = b;
        c0 = showingRowLabels ? 1 : 0;
        tableSource.fireTableStructureChanged();
    }

    /**
     * Accessor method for flag indicating if table should include a column of
     * "x" values.
     */
    public boolean isShowingRowLabels() {
        return showingRowLabels;
    }

    /**
     * Accessor method of the precision, which specifies the number of
     * significant figures to be displayed.
     */
    public int getPrecision() {
        return digits;
    }

    /**
     * Accessor method of the precision, which specifies the number of
     * significant figures to be displayed.
     */
    public void setPrecision(int n) {
        if (n < 1) {
            throw new RuntimeException("precision must be positive");
        }
        digits = n;
        if (digits+4 > nColumns) {
            setNumDigitColumns(digits+4);
        }
    }
    
    public int getNumDigitColumns() {
        return nColumns;
    }
    
    public void setNumDigitColumns(int n) {
        nColumns = n;
        minFloat = Math.pow(10, digits-nColumns);
        maxFloat = Math.pow(10, nColumns);
    }

    public void repaint() {
        table.repaint();
    }

    public Component graphic(Object obj) {
        return panel;
    }

    public AbstractTableModel tableModel() {
        return tableSource;
    }

    /**
     * Descriptive text for each row.
     */
    public String[] getRowLabels() {
        return rowLabels;
    }

    /**
     * Sets values for descriptive text for each row (as defined before
     * transpose). Display of labels is determined by value of showingRowLabels
     * (if not transposing) or showingColumnHeaders (if transposing).
     * 
     * Setting the row labels will persist.  To cause the row labels to be
     * taken from the IDataInfo again after calling setRowLabels, setRowLabels
     * must be called again with null argument.
     */
    public void setRowLabels(String[] newRowLabels) {
        rowLabels = newRowLabels;
        autoRowLabels = rowLabels == null;
        if (autoRowLabels) {
            rowLabels = new String[dataTable.getRowCount()];
            tableRowHeadersChanged(null);
        }
    }

    public void setUnit(DataTag[] dataTags, Unit newUnit) {
        unitList.add(new DataTagBag(dataTags, newUnit));

        // now go through and look for a current Data with these tags, but
        // don't use these tags if a previous set of tags also matches.
        for(int i=0; i<units.length; i++) {
            // if the user specified a unit for this data specifically, use it.
            DataTagBag tagUnit = DataTagBag.getDataTagBag(unitList, dataTable.getDataInfo(i).getTags());
            if (tagUnit != null) {
                units[i]= (Unit)tagUnit.object;
            }
        }
    }
    
    /**
     * Reconstruct the units array.  Manually-set units for columns are stored
     * in a HashMap (keyed to the column), so this method looks for units in 
     * the HashMap and if that fails, uses the column's dimension's default unit.
     * Old units (from old columns) are not removed from the hash, because the 
     * DataSinkTables doesn't know with certainty that it no longer contains a 
     * column.  If the column magically comes back, we need to use any 
     * previously-set unit.
     */
    protected void recomputeUnits() {
        units = new Unit[dataTable.getDataCount()];
        for (int i=0; i<units.length; i++) {
            IDataInfo columnInfo = dataTable.getDataInfo(i);
            Unit dataUnit = defaultUnit == null ? columnInfo.getDimension().getUnit(UnitSystem.SIM) : defaultUnit;

            DataTagBag tagUnit = DataTagBag.getDataTagBag(unitList, dataTable.getDataInfo(i).getTags());
            if (tagUnit != null) {
                dataUnit = (Unit)tagUnit.object;
            }
            units[i] = dataUnit;
        }
    }
    
    /**
     * Reconstruct the column headers array.  Manually-set headers for columns
     * are used if they were set and if not, uses the column's own header is used.
     */
    protected void recomputeColumnHeaders() {
        columnHeaders = new String[dataTable.getDataCount()];
        for (int i=0; i<columnHeaders.length; i++) {
            IDataInfo columnInfo = dataTable.getDataInfo(i);
            String header = columnInfo.getLabel();
            String suffix = "";
            if (showingUnits) {
                suffix = units[i].symbol();
                if (!suffix.equals(""))
                    suffix = " (" + suffix + ")";
            }
            header += suffix;

            DataTagBag tagHeader = DataTagBag.getDataTagBag(columnHeaderList, dataTable.getDataInfo(i).getTags());
            if (tagHeader != null) {
                header = (String)tagHeader.object;
            }
            columnHeaders[i] = header;
        }
    }
    
    /**
     * Sets the units of all columns to the given unit.
     */
    public void setAllUnits(Unit newUnit) {
        defaultUnit = newUnit;
    }

    /**
     * Returns the unit assign to the i-th column (as defined before any
     * transpose).
     * 
     * @param i
     *            column index, numbered from zero, not including the row-label
     *            column
     * @return the unit of the indicated column
     */
    public Unit getUnit(int i) {
        return units[i];
    }

    /**
     * @return Returns the showingColumnHeaders flag.
     */
    public boolean isShowingColumnHeaders() {
        return showingColumnHeaders;
    }

    /**
     * @param showingColumnHeaders
     *            The showingColumnHeaders flag to set.
     */
    public void setShowingColumnHeaders(boolean showingColumnHeaders) {
        this.showingColumnHeaders = showingColumnHeaders;
    }

    /**
     * @return Returns the showingUnits flag, which indicates whether a units
     *         suffix is applied to the column headers (or row labels, if
     *         transposing).
     */
    public boolean isShowingUnits() {
        return showingUnits;
    }

    /**
     * Flag indicating whether a units suffix applied to the column headers (or
     * row labels, if transposing). Default is true.
     * 
     * @param showingUnits
     *            The showingUnits flag to set.
     */
    public void setShowingUnits(boolean showingUnits) {
        this.showingUnits = showingUnits;
    }

    /**
     * Returns the header of the row-label column. Not used if showingRowLabels
     * is false. Default is empty string.
     */
    public String getRowLabelColumnHeader() {
        return rowLabelColumnHeader;
    }

    /**
     * Sets the header of the row-label column. Not used if showingRowLabels is
     * false.
     */
    public void setRowLabelColumnHeader(String rowLabelColumnHeader) {
        this.rowLabelColumnHeader = rowLabelColumnHeader;
    }
    
    public void setColumnHeader(DataTag[] tags, String newHeader) {
        columnHeaderList.add(new DataTagBag(tags, newHeader));

        // now go through and look for a current Data with these tags, but
        // don't use these tags if a previous set of tags also matches.
        for(int i=0; i<columnHeaders.length; i++) {
            // if the user specified a unit for this data specifically, use it.
            DataTagBag tagHeader = DataTagBag.getDataTagBag(columnHeaderList, dataTable.getDataInfo(i).getTags());
            if (tagHeader != null) {
                columnHeaders[i]= (String)tagHeader.object;
            }
        }
    }

    private final JTable table;
    protected final DataSinkTable dataTable;
    private final MyTable tableSource;
    protected final JPanel panel;

    protected boolean showingRowLabels;
    protected boolean showingColumnHeaders;
    protected boolean showingUnits = true;
    private boolean fitToWindow = true;
    protected boolean transposed;
    protected int c0 = 0;
    protected String[] rowLabels = new String[0];
    protected boolean autoRowLabels = true;
    protected String rowLabelColumnHeader = "";
    protected String[] columnHeaders = new String[0];
    protected Unit[] units = new Unit[0];
    private LinkedList<DataTagBag> unitList;
    private LinkedList<DataTagBag> columnHeaderList;
    private Unit defaultUnit;
    protected int digits;
    protected int nColumns;
    protected double maxFloat, minFloat;
    

    //structures used to adjust precision of displayed values
    //  private final java.text.NumberFormat formatter =
    // java.text.NumberFormat.getInstance();
    protected final Formatter formatter;
    protected final StringBuilder stringBuilder;
//    protected final java.text.NumberFormat formatter = new etomica.util.ScientificFormat();

    private class MyTable extends AbstractTableModel {

        public Object getValueAt(int row, int column) {
            if (showingRowLabels && column == 0) {
                if (transposed) {
                    return columnLabel(row);
                }
                return (row < rowLabels.length) ? rowLabels[row] : "";
            }
            //r and c are the row/column indices for the internal
            // representation
            //of the table
            int r = transposed ? column - c0 : row;
            int c = transposed ? row : column - c0;
            DataDoubleArray columnData = (DataDoubleArray)dataTable.getData(c);
            double value = Double.NaN;
            if (r < columnData.getLength()) {
                value = units[c].fromSim(dataTable.getValue(r, c));
            }
            if (value == Double.NEGATIVE_INFINITY) {
                return "-infinity";
            }
            if (value == Double.POSITIVE_INFINITY) {
                return "infinity";
            } 
            if (Double.isNaN(value)){
            	return "NaN";
            }
            double av = Math.abs(value);
            if ((av < minFloat && value != 0) || av >= maxFloat) {
                formatter.format("%"+(3+digits)+"."+(digits-1)+"e", value);
            }
            else if (av < 1) {
                formatter.format("%"+(digits+5)+"."+(digits+3)+"f", value);
            }
            else {
                int n = (int)(Math.log(av)/Math.log(10));
                formatter.format("%"+(digits+5)+"."+(digits+3-n)+"f", value);
            }
            String out = formatter.out().toString();
            stringBuilder.delete(0, stringBuilder.length());
            return out;
        }

        public int getRowCount() {
            int n = transposed ? dataTable.getDataCount() : dataTable
                    .getRowCount();
            return n;
        }

        public int getColumnCount() {
            int n = c0
                    + (transposed ? dataTable.getRowCount() : dataTable
                            .getDataCount());
            return n;
        }

        public String getColumnName(int i) {
            if (i == 0 && showingRowLabels)
                return rowLabelColumnHeader;
            int c = i - c0;
            if (transposed) {
                return (rowLabels.length > c) ? rowLabels[c] : "";
            }
            return columnLabel(c);
        }

        private String columnLabel(int i) {
            return columnHeaders[i];
        }
    }

    private class InputEventHandler extends MouseAdapter implements KeyListener {

        public void keyPressed(KeyEvent evt) {
//            char c = evt.getKeyChar();
//            if (Character.isDigit(c)) {
//                System.out.println("Setting precision "
//                        + Character.getNumericValue(c));
//                setPrecision(Character.getNumericValue(c));
//                panel.repaint();
//            }
        }

        public void keyReleased(KeyEvent evt) {
        }

        public void keyTyped(KeyEvent evt) {
        }

        public void mouseClicked(MouseEvent evt) {
            panel.requestFocus();
        }

        public void mouseEntered(MouseEvent evt) {
            panel.requestFocus();
        }

        public void mouseExited(MouseEvent evt) {
            panel.transferFocus();
        }
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
    	final String APP_NAME = "Display Table";

    	etomica.space.Space sp = etomica.space2d.Space2D.getInstance();
        final HSMD2D sim = new HSMD2D();
        final SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
        sim.integrator.setIsothermal(true);

        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        MeterNMolecules nMeter = new MeterNMolecules();
        nMeter.setBox(sim.box);
        DataTableAverages dataTable = new DataTableAverages(sim.integrator);
        dataTable.addDataSource(pMeter);
        dataTable.addDataSource(nMeter);
        DisplayTable table = new DisplayTable(dataTable);

        table.setRowLabels(new String[] { "Current", "Average", "Error" });
        table.setTransposed(false);
        table.setShowingRowLabels(true);
//        table.setPrecision(7);

        graphic.getController().getReinitButton().setPostAction(graphic.getPaintAction(sim.box));

        graphic.add(table);
        graphic.makeAndDisplayFrame(APP_NAME);
    }//end of main

}
