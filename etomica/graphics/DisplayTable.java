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
import etomica.DataSink;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageSegment;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressureHard;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.utility.Arrays;

/**
 * Presents data in a tabular form.
 *
 * @author David Kofke
 * @see DisplayTableFunction
 */
 
 /* History of changes
  * 7/20/02 Added key listener to set precision of displayed values
  */
 
public class DisplayTable extends Display implements EtomicaElement {
    
    private JTable table;
    private JPanel panel;
    private boolean showLabels = true;
    private boolean fitToWindow = true;
    private MyTableData tableSource;
    private DataGroup[] data = new DataGroup[0];
    private boolean fillHorizontal;
    
        //structures used to adjust precision of displayed values
//	private final java.text.NumberFormat formatter = java.text.NumberFormat.getInstance();
	private final java.text.NumberFormat formatter = new etomica.utility.ScientificFormat();
    private final javax.swing.table.DefaultTableCellRenderer numberRenderer =
        new javax.swing.table.DefaultTableCellRenderer() { 
            public void setValue(Object value) { 
        		setText((value == null) ? "" : formatter.format(value)); 
	        }
        };
        
    public DisplayTable() {
        setupDisplay();
        setLabel("Data");
        numberRenderer.setHorizontalAlignment(javax.swing.JLabel.RIGHT);
        setPrecision(4);
        InputEventHandler listener = new InputEventHandler();
        panel.addKeyListener(listener);
        panel.addMouseListener(listener);
        setFillHorizontal(true);
        setupDisplay();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabular display of data from several sources");
        return info;
    }
    
    /**
     * Implementation of parent's abstract method to set up table whenever 
     * data sources are added or changed.
     */
    protected void setupDisplay() {
        if(panel == null) panel = new javax.swing.JPanel(new java.awt.FlowLayout());
        else panel.removeAll();
        panel.setSize(100,150);
        tableSource = new MyTableData();   //inner class, defined below
        table = new JTable(tableSource);
        if(!fitToWindow) table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        table.setDefaultRenderer(Number.class, numberRenderer);
        table.setDefaultRenderer(Double.class, numberRenderer);
        panel.add(new JScrollPane(table));
    }
        
    /**
     * If true, each data set fills a row; if false, data are arranged to fill columns.
     */
    public boolean isFillHorizontal() {
        return fillHorizontal;
    }
    /**
     * If fillHorizontal is true, each data set fills a row; if false, data are arranged to fill columns.
     */
    public void setFillHorizontal(boolean fillHorizontal) {
        this.fillHorizontal = fillHorizontal;
    }
    /**
     * If true, columns will be squeezed so table fits to window; if false, table
     * will exceed window size, and to view full width must be placed in a scroll pane.
     * Default is true.
     */
    public void setFitToWindow(boolean b) {
        fitToWindow = b;
        setupDisplay();
    }
    /**
     * Accessor method for flag determining if table is squeezed to fit in window.
     */
    public boolean isFitToWindow() {return fitToWindow;}
    
    /**
     * Mutator method for whether table should include a column of "x" values.
     */
    public void setShowLabels(boolean b) {
        showLabels = b;
        setupDisplay();
    }
    /**
     * Accessor method for flag indicating if table should include a column of "x" values.
     */
    public boolean isShowLabels() {return showLabels;}
            
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public int getPrecision() {return formatter.getMaximumFractionDigits();}
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public void setPrecision(int n) {
	    formatter.setMaximumFractionDigits(n);
    }
    
//    public Action makeLogTableAction() {return new LogTableAction();}
    
    public void repaint() {table.repaint();}

    public Component graphic(Object obj) {return panel;}
    
    public AbstractTableModel tableModel() {return tableSource;}

    private class MyTableData extends AbstractTableModel {
        
        String[] columnNames;
        Class[] columnClasses;
        int nColumns;
        int y0;  //index of first y column (1 if showing x, 0 otherwise)
        
        MyTableData() {
            y0 = showLabels ? 1 : 0;
            nColumns += y0;
            columnNames = new String[nColumns];
            columnClasses = new Class[nColumns];
            if(showLabels) {
                columnNames[0] = "Property";
                columnClasses[0] = String.class;
            }
            for(int i=y0; i<nColumns; i++) {
                columnNames[i] = "";
                columnClasses[i] = Double.class;
            }
        }
        
        public Object getValueAt(int row, int column) {
            int r = fillHorizontal ? row : column;
            int c = fillHorizontal ? column : row;
//            if(showLabels && column == 0) return labels[row];
            return new Double(data[r].unit.fromSim(data[r].y[c]));
        }
        
        public int getRowCount() {return fillHorizontal ? data.length : nColumns;}
        public int getColumnCount() {return fillHorizontal ? nColumns : data.length;}
                
        public String getColumnName(int column) {return columnNames[column];}
        public Class getColumnClass(int column) {return columnClasses[column];}
    }
    
    public DataSink makeDataSink(Unit u) {
        DataGroup newGroup = (DataGroup)makeDataSink();
        newGroup.unit = u;
        data = (DataGroup[])Arrays.addObject(data, newGroup);
        return newGroup;
    }
    
    public DataSink makeDataSink() {
        DataGroup newGroup = new DataGroup();
        data = (DataGroup[])Arrays.addObject(data, newGroup);
        return newGroup;
    }
    
    public void removeDataSink(DataSink sink) {
        data = (DataGroup[])Arrays.removeObject(data, sink);
    }
    
    public void removeAllSinks() {
        data = new DataGroup[0];
    }

    private class DataGroup implements DataSink {
        double[] y;
        Dimension dimension = Dimension.UNDEFINED;
        Unit unit = Unit.UNDEFINED;
        String label = "";
        
        DataGroup() {
            y = new double[0];
        }

        public void putData(double[] values) {
            if(y.length != values.length) y = (double[])values.clone();
            else System.arraycopy(values, 0, y, 0, values.length);
            repaint();
        }

        public void setDimension(Dimension dimension) {
            this.dimension = dimension;
        }

        public void setLabel(String label) {
            this.label = label;
        }
    }

    
    private class InputEventHandler extends MouseAdapter implements KeyListener {
        
        public void keyPressed(KeyEvent evt) {
            char c = evt.getKeyChar();
            if(Character.isDigit(c)) {
            	System.out.println("Setting precision "+Character.getNumericValue(c));
                setPrecision(Character.getNumericValue(c));
                panel.repaint();
            }
        }
        public void keyReleased(KeyEvent evt) {}
        public void keyTyped(KeyEvent evt) {}
        
        public void mouseClicked(MouseEvent evt) {
            panel.requestFocus();
        }
        public void mouseEntered(MouseEvent evt) {panel.requestFocus();}
        public void mouseExited(MouseEvent evt) {panel.transferFocus();}
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
        DisplayTable table = new DisplayTable();
        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        pMeter.setPhase(sim.phase);
        AccumulatorAverageSegment pSegment = new AccumulatorAverageSegment(pMeter, sim.integrator, 
                new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE, AccumulatorAverage.ERROR},
                table.makeDataSink());
        MeterNMolecules nMeter = new MeterNMolecules();
        nMeter.setPhase(sim.phase);
        AccumulatorAverageSegment nSegment = new AccumulatorAverageSegment(nMeter, sim.integrator, 
                new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE, AccumulatorAverage.ERROR},
                table.makeDataSink());
        //end of unique part
                                            
        table.setPrecision(7);
        
        graphic.add(table);
        graphic.makeAndDisplayFrame();
    }//end of main

}