package etomica.graphics;
import etomica.*;

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.Box;
import javax.swing.JScrollPane;
import java.awt.event.*;

/**
 * Presents function data in a tabular form.  Data is obtained from DataSource objects.
 * Sources may be identified to the Display by passing an array of DataSource objects (in
 * the setDataSources method, or one at a time via the addDataSource method).
 *
 * @author David Kofke
 * @ee DisplayTable
 */
 
 /* History of changes
  * 7/20/02 Added key listener to set precision of displayed values
  */
public class DisplayTableFunction extends DisplayDataSources implements EtomicaElement {
    public String getVersion() {return "DisplayTableFunction:01.05.29/"+super.getVersion();}

    public JTable table;
    Box panel;
    private boolean transposed = false;
    private boolean showXColumn = true;
    private boolean fitToWindow = false;
    private double multiplier = 1.0;//alterntive way to adjust precision
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
    
        
    public DisplayTableFunction() {
        this(Simulation.instance);
    }
    public DisplayTableFunction(Simulation sim)  {
        super(sim);
        setupDisplay();
        setLabel("Function");
        setWhichValue(MeterAbstract.AVERAGE);
        table = new JTable();
        numberRenderer.setHorizontalAlignment(javax.swing.JLabel.RIGHT);
        setPrecision(4);
        InputEventHandler listener = new InputEventHandler();
        panel.addKeyListener(listener);
        panel.addMouseListener(listener);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Tabular display of function data from several sources");
        return info;
    }
    
//    public void setTransposed(boolean b) {transposed = b; setupDisplay();}
//    public boolean isTransposed() {return transposed;}
    
    /**
     * Implementation of parent's abstract method to set up table whenever 
     * data sources are added or changed.
     */
    protected void setupDisplay() {
        if(panel == null) panel = Box.createVerticalBox();
        else panel.removeAll();
        panel.setSize(100,150);
        if(ySource == null || ySource.length == 0) return;
        tableSource = new MyTableData();   //inner class, defined below
        table = new JTable(tableSource);
        if(!fitToWindow) table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        table.setDefaultRenderer(Number.class, numberRenderer);
        panel.add(new JScrollPane(table));
        //trying to get table to redisplay correctly when 
        //changing whichValue.  Still not working correctly.
        panel.invalidate();
        panel.validate();
        doUpdate();
        panel.repaint();        
    }
    
    /**
     * If true, columns will be squeezed so table fits to window; if false, table
     * will exceed window size, and to view full width must be placed in a scroll pane.
     * Default is false.
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
    public void setShowXColumn(boolean b) {
        showXColumn = b;
        setupDisplay();
    }
    /**
     * Accessor method for flag indicating if table should include a column of "x" values.
     */
    public boolean isShowXColumn() {return showXColumn;}
    
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public void setPrecision(int n) {
	    formatter.setMaximumFractionDigits(n);
//	    multiplier = Math.pow(10.0, n);
    }
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public int getPrecision() {return formatter.getMaximumFractionDigits();}
    
    public Action makeLogTableAction() {return new LogTableAction();}
    
    public void repaint() {table.repaint();}

    public Component graphic(Object obj) {return panel;}

    private class MyTableData extends AbstractTableModel {
        
        String[] columnNames;
        Class[] columnClasses;
        int nColumns;
        int y0;  //index of first y column (1 if showing x, 0 otherwise)
        
        MyTableData() {
 //           if(!transposed) {
                nColumns = ySource.length;
                y0 = showXColumn ? 1 : 0;
                nColumns += y0;
                columnNames = new String[nColumns];
                columnClasses = new Class[nColumns];
                if(showXColumn) {
                    if(xSource != null) {
                        if(xSource instanceof DataSource.X) columnNames[0] = ((DataSource.X)xSource).getXLabel();
                        else columnNames[0] = xSource.getLabel();
                    }
                    else columnNames[0] = "";
                    columnClasses[0] = Double.class;
                }
                for(int i=y0; i<nColumns; i++) {
                    columnNames[i] = (ySource[i-y0] != null) ? ySource[i-y0].getLabel() : "";
                    columnClasses[i] = Double.class;
                }
     //       } else {//transposed
     //           nColumns = ?
        }
        
            //x-y values are updated in update() method
            //because we don't want to call currentValue
            //or average for each function entry
        public Object getValueAt(int row, int column) {
            if(showXColumn && column == 0) return new Double(multiplier*x[row]);
            else return new Double(multiplier*y[column-y0][row]);
        }
        
        public int getRowCount() {return (y != null && y[0] != null) ? y[0].length : 0;}
        public int getColumnCount() {return y0 + ((ySource != null) ? ySource.length : 0);}
                
        public String getColumnName(int column) {return columnNames[column];}
        public Class getColumnClass(int column) {return columnClasses[column];}
    }
    
    private class InputEventHandler extends MouseAdapter implements KeyListener {
        
        public void keyPressed(KeyEvent evt) {
            char c = evt.getKeyChar();
            if(Character.isDigit(c)) {
                setPrecision(Character.getNumericValue(c));
                panel.repaint();
                System.out.println("Changing table precision to "+c);
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
        Default.ATOM_SIZE = 1.0;                   
	    IntegratorHard integratorHard = new IntegratorHard();
	    SpeciesSpheres speciesSpheres = new SpeciesSpheres();
	    speciesSpheres.setNMolecules(300);
	    Phase phase = new Phase();
	    Potential2 potential = new P2HardSphere();
	    Controller controller = new Controller();
	    DisplayPhase displayPhase = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard.new Timer(integratorHard.chronoMeter());
	    timer.setUpdateInterval(10);
        //part that is unique to this demonstration
        Default.BLOCK_SIZE = 20;
        MeterFunction rdf = new MeterRDF();
        rdf.setActive(true);
        DisplayTableFunction rdfTable = new DisplayTableFunction();
        //end of unique part
        
		Simulation.instance.elementCoordinator.go();
		                                    
        potential.setIterator(new AtomPairIterator(phase));
//        potential.set(species.getAgent(phase));
		
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);
        Simulation.makeAndDisplayFrame(Simulation.instance);
        		                                    
        rdfTable.setDataSources(rdf);
//        rdfTable.setX(rdf.X());
    }//end of main
*/
}