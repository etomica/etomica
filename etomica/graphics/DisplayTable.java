package etomica.graphics;

import java.awt.Component;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;

import etomica.Action;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.AccumulatorAverage;
import etomica.data.DataBin;
import etomica.data.DataTable;
import etomica.data.DataTableAverages;
import etomica.data.DataTableListener;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressureHard;
import etomica.units.Unit;
import etomica.utility.Arrays;

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
        this(new DataTable());
    }

    public DisplayTable(DataTable dataTable) {
        this.dataTable = dataTable;

        units = new Unit[dataTable.getColumnCount()];
        for (int i = 0; i < units.length; i++) {
            units[i] = dataTable.getColumn(i).getDimension().defaultIOUnit();
        }

        dataTable.addTableListener(this);
        tableSource = new MyTable(); //inner class, defined below
        table = new JTable(tableSource);
        panel = new javax.swing.JPanel(new java.awt.FlowLayout());

        setLabel("Data");
        numberRenderer.setHorizontalAlignment(javax.swing.JLabel.RIGHT);
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

    public DataTable getDataTable() {
        return dataTable;
    }

    /**
     * Causes the display of the plot to be updated.
     */
    public void tableDataChanged() {
        tableSource.fireTableDataChanged();
        repaint();
    }

    /**
     * Updates the units array for the new column, using the default units.
     */
    public void tableColumnAdded(DataBin newColumn) {
        units = (Unit[]) Arrays.addObject(units, newColumn.getDimension()
                .defaultIOUnit());
        tableSource.fireTableStructureChanged();
    }

    /**
     * Causes the corresponding units element to be removed.
     */
    public void tableColumnRemoved(int index, DataBin oldColumn) {
        units = (Unit[]) Arrays.removeObject(units, units[index]);
        tableSource.fireTableStructureChanged();
    }

    /**
     * Has no effect. Part of the DataTableListener interface.
     */
    public void tableRowCountChanged(int oldCount, int newCount) {
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

    //    public Action makeLogTableAction() {return new LogTableAction();}

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

    /**
     * Sets the units of all columns to the given unit.
     */
    public void setAllUnits(Unit newUnit) {
        for (int i = 0; i < units.length; i++) {
            units[i] = newUnit;
        }
    }

    /**
     * Sets the unit of the i-th column (as defined before any transpose).
     * 
     * @param i
     *            column index, numbered from zero, not including the row-label
     *            column
     * @param newUnit
     *            the unit to be assigned to the indicated column
     */
    public void setUnit(int i, Unit newUnit) {
        units[i] = newUnit;
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
    private final DataTable dataTable;
    private final MyTable tableSource;
    private final JPanel panel;

    private boolean showingRowLabels;
    private boolean showingColumnHeaders;
    private boolean showingUnits = true;
    private boolean fitToWindow = true;
    private boolean transposed;
    private int c0 = 0;
    private String[] rowLabels = new String[0];
    private String rowLabelColumnHeader = "";
    private Unit[] units = new Unit[0];

    //structures used to adjust precision of displayed values
    //  private final java.text.NumberFormat formatter =
    // java.text.NumberFormat.getInstance();
    private final java.text.NumberFormat formatter = new etomica.utility.ScientificFormat();
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
                } else {
                    return (row < rowLabels.length) ? rowLabels[row] : "";
                }
            }
            //r and c are the row/column indices for the internal
            // representation
            //of the table
            int r = transposed ? column - c0 : row;
            int c = transposed ? row : column - c0;
            DataBin col = dataTable.getColumn(c);
            double value = Double.NaN;
            if (r < col.getDataLength()) {
                value = units[c].fromSim(col.getData()[r]);
            }
            return new Double(value);
        }

        public int getRowCount() {
            int n = transposed ? dataTable.getColumnCount() : dataTable
                    .getRowCount();
            return n;
        }

        public int getColumnCount() {
            int n = c0
                    + (transposed ? dataTable.getRowCount() : dataTable
                            .getColumnCount());
            return n;
        }

        public String getColumnName(int i) {
            if (i == 0 && showingRowLabels)
                return rowLabelColumnHeader;
            int c = i - c0;
            if (transposed) {
                return (rowLabels.length > c) ? rowLabels[c] : "";
            } else {
                return columnLabel(c);
            }
        }

        private String columnLabel(int i) {
            String suffix = "";
            if (showingUnits) {
                suffix = units[i].symbol();
                if (!suffix.equals(""))
                    suffix = "(" + suffix + ")";
            }
            return dataTable.getColumn(i).getLabel()
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

    private class LogTableAction implements Action {

        public void actionPerformed() {
            etomica.log.LogTable log = new etomica.log.LogTable();
            log.writeTable(tableSource, "output.xls");
        }

        public String getLabel() {
            return "Write table to log file";
        }
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.integrator.setIsothermal(true);

        //part that is unique to this demonstration
        //        DisplayTable table = new DisplayTable();
        //        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        //        pMeter.setPhase(sim.phase);
        //        AccumulatorAverageSegment pSegment = new AccumulatorAverageSegment(
        //                pMeter, sim.integrator, new AccumulatorAverage.Type[] {
        //                        AccumulatorAverage.MOST_RECENT,
        //                        AccumulatorAverage.AVERAGE,
        //                        AccumulatorAverage.ERROR },
        //                table.getDataTable().makeColumn(pMeter.getDimension()));
        //        MeterNMolecules nMeter = new MeterNMolecules();
        //        nMeter.setPhase(sim.phase);
        //        AccumulatorAverageSegment nSegment = new AccumulatorAverageSegment(
        //                nMeter, sim.integrator, new AccumulatorAverage.Type[] {
        //                        AccumulatorAverage.MOST_RECENT,
        //                        AccumulatorAverage.AVERAGE,
        //                        AccumulatorAverage.ERROR },
        //                table.getDataTable().makeColumn(nMeter.getDimension()));
        //end of unique part

        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        pMeter.setPhase(sim.phase);
        MeterNMolecules nMeter = new MeterNMolecules();
        nMeter.setPhase(sim.phase);
        DataTableAverages dataTable = new DataTableAverages(sim.integrator,
                new AccumulatorAverage.Type[] { AccumulatorAverage.MOST_RECENT,
                        AccumulatorAverage.AVERAGE, AccumulatorAverage.ERROR });
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