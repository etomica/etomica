package etomica;
import etomica.units.Dimension;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.Box;
import java.awt.event.*;
import java.awt.*;
import javax.swing.JScrollPane;
import etomica.utility.Histogram;


/**
 * Meter for recording and averaging tensor data.  Includes a DisplayTensorTable method to show results.
 * 
 * @author Rob Riggleman
 */

public abstract class MeterTensor extends MeterAbstract {
    
    MeterAbstract.Accumulator[][] accumulator;
    Space.Tensor values, average, error, y;
    private int D;
    
    
    public MeterTensor(Simulation sim) {
        super(sim);
        D = sim.space().D();
        accumulator = new MeterAbstract.Accumulator[D][D];
        setActive(true);
        for (int i=0; i<D; i++) {
            for (int j=0; j<D; j++) {
                accumulator[i][j] = new MeterAbstract.Accumulator();
            }
        }
        values = sim.space().makeTensor();
        average = sim.space().makeTensor();
        error = sim.space().makeTensor();
        y = sim.space().makeTensor();
    }    
    
    public abstract Space.Tensor currentValue();
    
    /**
     * Current value of the i,j component of the virial tensor
     */
    public double currentValue(int i, int j) {
        return currentValue().component(i, j);
    }
    
    public void updateSums() {
        values.E(currentValue());
        for (int i=0; i<accumulator.length; i++) {
            for (int j=0; j<accumulator.length; j++) {
                accumulator[i][j].add(values.component(i, j));
            }
        }
    }
        
    
    public Space.Tensor average() {
        for (int i=0; i<D; i++) {
            for (int j=0; j<D; j++) {
                average.setComponent(i, j, accumulator[i][j].average());
            }
        }
        return average;
    }
    
    public double average(int i, int j) {
        return average().component(i, j);        
    }
    
    public Space.Tensor error() {
        for (int i=0; i<D; i++) {
            for (int j=0; j<D; j++) {
                error.setComponent(i, j, accumulator[i][j].error());
            }
        }
        return error;
    }
    
    public double error(int i, int j) {
        return error().component(i, j);
    }
    
    public Space.Tensor mostRecent() {
        for (int i=0; i<accumulator.length; i++) {
            for (int j=0; j<accumulator.length; j++) {
                y.setComponent(i, j, accumulator[i][j].mostRecent());
            }
        }
        return y;
    }
    
	/**
	* Accessor method to indicate if the meter should keep a histogram of all measured values.
	* Histogramming is not enabled for this type of meter (MeterTensor).
	*/
	public boolean isHistogramming() {return false;}
    	 
	/**
	* Accessor method to indicate if the meter should keep a histogram of all measured values.
	* Histogramming is not enabled for this type of meter (MeterTensor).
	*/
	public void setHistogramming(boolean b) {
	 //perform no action
	}

    public void reset() {
        for (int i=0; i<accumulator.length; i++) {
            for (int j=0; j<accumulator.length; j++) {
                accumulator[i][j].reset();
            }
        }
    }
    
	public interface Atomic {
	    public Space.Tensor currentValue(Atom a);
	}
	
	public interface Collisional extends IntegratorHard.CollisionListener {
	    public Space.Tensor collisionValue(AtomPair pair, Potential.Hard p);
	}
    
    
   /**
    * Creates a table to display the values of the tensor.  
    */    
    
    public class DisplayTensorTable extends Display {
        public JTable table;
        MyTableData dataSource;
        MeterTensor meterT;
        Box panel = Box.createVerticalBox();
        Button resetButton = new Button("Reset averages");
        private boolean showAverages = true;
        
        public DisplayTensorTable() {this(true);
            SymFocus aSymFocus = new SymFocus();
            this.addFocusListener(aSymFocus);
        }
        
        public DisplayTensorTable(boolean showAvgs) {
            super(MeterTensor.this.parentSimulation());
            showAverages = showAvgs;
            panel.setSize(100,150);
            dataSource = new MyTableData();
            table = new JTable(dataSource);
            panel.add(new JScrollPane(table));
            if (showAverages) {
                resetButton.addActionListener(new ActionListener() {
                    public void actionPerformed(java.awt.event.ActionEvent event) {DisplayTensorTable.this.resetAverages();}
                });
                panel.add(resetButton);
            }
            add(panel);
        }
        
        public void resetAverages() {
            meterT.reset();
        }
        public void setResetVisible(boolean b) {resetButton.setVisible(b);}
        public boolean getResetVisible() {return resetButton.isVisible();}
        
        public void setMeter(MeterTensor mt) {
            meterT = mt;
        }
        
        public void doUpdate() {}
        public void repaint() {table.repaint();}
        
        private class MyTableData extends AbstractTableModel {
            String[] columnNames = new String[] {"Property", "Average", "Error"};
            Class[] columnClasses = new Class[] {String.class, Double.class, Double.class};
            
            MyTableData() {
                if(D==1) {
                    if (showAverages) {
                        columnNames = new String[] {"Property", "Average x"};
                        columnClasses = new Class[] {String.class, Double.class};
                    }
                    else {
                        columnNames = new String[] {"Property", "Current x"};
                        columnClasses = new Class[] {String.class, Double.class};
                    }
                }
                else if(D==2) {
                    if (showAverages) {
                        columnNames = new String[] {"Property", "Average x", "Average y"};
                        columnClasses = new Class[] {String.class, Double.class, Double.class};
                    }
                    else {
                        columnNames = new String[] {"Property", "Current x", "Current y"};
                        columnClasses = new Class[] {String.class, Double.class, Double.class};
                    }
                }
                else if (D == 3) {
                    if (showAverages) {
                        columnNames = new String[] {"Property", "Average x", "Average y", "Average z"};
                        columnClasses = new Class[] {String.class, Double.class, Double.class, Double.class};
                    }
                    else {
                        columnNames = new String[] {"Property", "Current x", "Current y", "Current z"};
                        columnClasses = new Class[] {String.class, Double.class, Double.class, Double.class};
                    }
                }
            }
            
            public Object getValueAt(int row, int column) {
                MeterTensor mt = meterT;
                switch(column) {
                    case 0: {
                        if (D == 1) {
                            if(showAverages) {return mt.getLabel() + "              Average x";}
                            else { return mt.getLabel() + "              Current x";}
                        }
                        else if (D == 2) {
                            if (showAverages) {
                                if (row == 1) {return mt.getLabel() + "               Average y";}
                                else {return mt.getLabel() + "               Average x";}
                            }
                            else {
                                if(row == 1) {return mt.getLabel() + "               Current x";}
                                else {return mt.getLabel() + "               Current y";}
                            }
                        }
                        else {
                            if (showAverages) {
                                if (row == 1) {return mt.getLabel() + "               Average x";}
                                else if (row == 2) {return mt.getLabel() + "               Average y";}
                                else {return mt.getLabel() + "               Average z";}
                            }
                            else {
                                if (row == 1) {return mt.getLabel() + "               Current x";}
                                else if (row == 2) {return mt.getLabel() + "               Current y";}
                                else {return mt.getLabel() + "               Current z";}
                            }
                        }
                    }
                    case 1: return new Double(showAverages ? average(row, column - 1) : currentValue(row, column - 1));
                    case 2: return new Double(showAverages ? average(row, column - 1) : currentValue(row, column - 1));
                    default: return null;
                }
            }
            
            public int getRowCount() {return D;}
            public int getColumnCount() {return D + 1;}
            
            public String getColumnName(int column) {return columnNames[column];}
            public Class getColumnClass(int column) {return columnClasses[column];}
            
        }
        
        class SymFocus extends java.awt.event.FocusAdapter {
            public void focusGained(java.awt.event.FocusEvent event) {
                Object object = event.getSource();
                if (object == DisplayTensorTable.this)
                    DisplayTable_FocusGained(event);
            }
        }
        
        void DisplayTable_FocusGained(java.awt.event.FocusEvent event) {
        }
    }
}
    
    



