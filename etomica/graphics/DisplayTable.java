package etomica.graphics;
import etomica.*;

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.Box;
import javax.swing.JScrollPane;
import etomica.units.Unit;
import java.awt.event.*;

/**
 * Presents data in a tabular form.  Data is obtained from DatumSource objects.
 * Sources may be identified to the Display by passing an array of DatumSource objects (in
 * the setDatumSources method, or one at a time via the addDatumSource method).
 *
 * @author David Kofke
 * @see DisplayTableFunction
 */
 
 /* History of changes
  * 7/20/02 Added key listener to set precision of displayed values
  */
 
public class DisplayTable extends DisplayDatumSources implements EtomicaElement
{
    public String getVersion() {return "DisplayTable:01.05.29/"+super.getVersion();}

    public JTable table;
    javax.swing.JPanel panel;
    private boolean showLabels = true;
    private boolean fitToWindow = true;
    private MyTableData tableSource;
    
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
        this(Simulation.instance);
    }
    public DisplayTable(Simulation sim)  {
        super(sim);
        setupDisplay();
        setLabel("Data");
        numberRenderer.setHorizontalAlignment(javax.swing.JLabel.RIGHT);
        setPrecision(4);
        InputEventHandler listener = new InputEventHandler();
        panel.addKeyListener(listener);
        panel.addMouseListener(listener);
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
   //     if(panel == null) panel = Box.createVerticalBox();
        if(panel == null) panel = new javax.swing.JPanel(new java.awt.FlowLayout());
        else panel.removeAll();
        panel.setSize(100,150);
        if(ySource == null || ySource.length == 0) return;
        tableSource = new MyTableData();   //inner class, defined below
        table = new JTable(tableSource);
        if(!fitToWindow) table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        table.setDefaultRenderer(Number.class, numberRenderer);
        panel.add(new JScrollPane(table));
        doUpdate();
        
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
    
    public Action makeLogTableAction() {return new LogTableAction();}
    
    public void repaint() {table.repaint();}

    public Component graphic(Object obj) {return panel;}

    private class MyTableData extends AbstractTableModel {
        
        String[] columnNames;
        Class[] columnClasses;
        int nColumns;
        int y0;  //index of first y column (1 if showing x, 0 otherwise)
        
        MyTableData() {
            nColumns = whichValues.length;
            y0 = showLabels ? 1 : 0;
            nColumns += y0;
            columnNames = new String[nColumns];
            columnClasses = new Class[nColumns];
            if(showLabels) {
                columnNames[0] = "Property";
                columnClasses[0] = String.class;
            }
            for(int i=y0; i<nColumns; i++) {
                columnNames[i] = (whichValues[i-y0]!=null) ? whichValues[i-y0].toString() : "";
                columnClasses[i] = Double.class;
            }
        }
        
            //x-y values are updated in update() method
            //because we don't want to call currentValue
            //or average for each function entry
        public Object getValueAt(int row, int column) {
            if(showLabels && column == 0) return labels[row];
            else return new Double(yUnit[row].fromSim(y[row][column-y0]));
        }
        
        public int getRowCount() {return ySource.length;}
        public int getColumnCount() {return y0 + whichValues.length;}
                
        public String getColumnName(int column) {return columnNames[column];}
        public Class getColumnClass(int column) {return columnClasses[column];}
    }
    
    private class InputEventHandler extends MouseAdapter implements KeyListener {
        
        public void keyPressed(KeyEvent evt) {
            char c = evt.getKeyChar();
            if(Character.isDigit(c)) {
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
        
    private class LogTableAction extends Action {
        
        public void actionPerformed() {
            etomica.log.LogTable log = new etomica.log.LogTable();
            log.writeTable(tableSource, "output.xls");
        }
    }

    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        
        //part that is unique to this demonstration
        Default.BLOCK_SIZE = 20;
        MeterPressureHard pMeter = new MeterPressureHard();
        MeterNMolecules nMeter = new MeterNMolecules();
        DisplayTable table = new DisplayTable();
        //end of unique part
                                            
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        table.setDatumSources(pMeter);
        table.addDatumSources(nMeter);
        table.setWhichValues(new MeterAbstract.ValueType[] {MeterAbstract.CURRENT, MeterAbstract.AVERAGE});
        table.setPrecision(7);
        
        Simulation.makeAndDisplayFrame(sim);
    }//end of main
*/
}