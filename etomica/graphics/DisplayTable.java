package etomica.graphics;

import java.awt.Component;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.LinkedList;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.AbstractTableModel;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.DataInfo;
import etomica.data.DataSet;
import etomica.data.DataSinkTable;
import etomica.data.DataTableAverages;
import etomica.data.DataTableListener;
import etomica.data.DataTag;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressureHard;
import etomica.data.types.DataDoubleArray;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Presents data in a tabular form.
 * 
 * @author David Kofke
 * @see DisplayTableFunction
 */

/*
 * History of changes 7/20/02 Added key listener to set precision of displayed
 * values
 */

public class DisplayTable extends Display implements DataTableListener,
        EtomicaElement {

    public DisplayTable() {
        this(new DataSinkTable());
    }

    public DisplayTable(DataSinkTable dataTable) {
        this.dataTable = dataTable;

        units = new Unit[dataTable.getDataCount()];
        for (int i = 0; i < units.length; i++) {
            units[i] = dataTable.getDataInfo(i).getDimension().getUnit(UnitSystem.SIM);
        }

        dataTable.addDataListener(this);
        tableSource = new MyTable(); //inner class, defined below
        table = new JTable(tableSource);
        panel = new javax.swing.JPanel(new java.awt.FlowLayout());

        setLabel("Data");
        numberRenderer.setHorizontalAlignment(SwingConstants.RIGHT);
        setPrecision(4);
        setTransposed(false);
        setShowingRowLabels(true);
        setShowingColumnHeaders(true);
        table.setDefaultRenderer(Number.class, numberRenderer);
        table.setDefaultRenderer(Double.class, numberRenderer);
        panel.add(new JScrollPane(table));
        InputEventHandler listener = new InputEventHandler();
        panel.addKeyListener(listener);
        panel.addMouseListener(listener);
        panel.setSize(100, 150);
        if (!fitToWindow)
            table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo(
                "Tabular display of data from several sources");
        return info;
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
        tableSource.fireTableStructureChanged();
    }

    /**
     * Part of the DataTableListener interface.  Updates the row headers.
     */
    public void tableRowCountChanged(DataSinkTable dummyTable) {
        rowLabels = new String[dataTable.getRowCount()];
        for (int i=0; i<rowLabels.length; i++) {
            rowLabels[i] = dataTable.getRowHeader(i);
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
        return formatter.getMaximumFractionDigits();
    }

    /**
     * Accessor method of the precision, which specifies the number of
     * significant figures to be displayed.
     */
    public void setPrecision(int n) {
        formatter.setMaximumFractionDigits(n);
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
     * @param rowLabels
     */
    public void setRowLabels(String[] rowLabels) {
        this.rowLabels = rowLabels;
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
            DataInfo columnInfo = dataTable.getDataInfo(i);
            Unit dataUnit = defaultUnit == null ? columnInfo.getDimension().getUnit(UnitSystem.SIM) : defaultUnit;

            DataTagBag tagUnit = DataTagBag.getDataTagBag(unitList, dataTable.getDataInfo(i).getTags());
            if (tagUnit != null) {
                dataUnit = (Unit)tagUnit.object;
            }
            units[i] = dataUnit;
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
    protected  String rowLabelColumnHeader = "";
    protected Unit[] units = new Unit[0];
    private LinkedList unitList;
    private Unit defaultUnit;
    

    //structures used to adjust precision of displayed values
    //  private final java.text.NumberFormat formatter =
    // java.text.NumberFormat.getInstance();
    protected final java.text.NumberFormat formatter = new etomica.util.ScientificFormat();
    private final javax.swing.table.DefaultTableCellRenderer numberRenderer = new javax.swing.table.DefaultTableCellRenderer() {

        public void setValue(Object value) {
            setText((value == null) ? "" : formatter.format(value));
        }
    };

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
            return new Double(value);
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
            String suffix = "";
            if (showingUnits) {
                suffix = units[i].symbol();
                if (!suffix.equals(""))
                    suffix = "(" + suffix + ")";
            }
            return dataTable.getDataInfo(i).getLabel()
                    + (showingUnits ? suffix : "");
        }

    }

    private class InputEventHandler extends MouseAdapter implements KeyListener {

        public void keyPressed(KeyEvent evt) {
            char c = evt.getKeyChar();
            if (Character.isDigit(c)) {
                System.out.println("Setting precision "
                        + Character.getNumericValue(c));
                setPrecision(Character.getNumericValue(c));
                panel.repaint();
            }
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
        etomica.simulation.prototypes.HSMD2D sim = new etomica.simulation.prototypes.HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.integrator.setIsothermal(true);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space);
        pMeter.setIntegrator(sim.integrator);
        MeterNMolecules nMeter = new MeterNMolecules();
        nMeter.setPhase(sim.phase);
        DataTableAverages dataTable = new DataTableAverages(sim,sim.integrator);
        dataTable.addDataSource(pMeter);
        dataTable.addDataSource(nMeter);
        DisplayTable table = new DisplayTable(dataTable);

        table.setRowLabels(new String[] { "Current", "Average", "Error" });
        table.setTransposed(false);
        table.setShowingRowLabels(true);
        table.setPrecision(7);

        graphic.add(table);
        graphic.makeAndDisplayFrame();
    }//end of main

}
