package simulate;

import simulate.*;
import java.beans.*;
import java.awt.*;
import javax.swing.JTable;
import javax.swing.table.*;

public class DisplayTable extends simulate.Display
{
    public JTable table;
    MyTableData dataSource;
    Meter[] meter = new Meter[0];
    int nMeters = 0;
    private boolean showingPE = false;  //flag to indicate if PotentialEnergy meter should be displayed
    private boolean showingKE = false;  //same for KineticEnergy meter
    
    public DisplayTable()
    {
        super();
        dataSource = new MyTableData();   //inner class, defined below
        table = new JTable(dataSource);
//        setVisible(false);  //get out of way of painting of table
//        add(table);
    }
    
    public void setShowingKE(boolean b) {showingKE = b;}
    public boolean isShowingKE() {return showingKE;}
    public void setShowingPE(boolean b) {showingPE = b;}
    public boolean isShowingPE() {return showingPE;}
    
    public void createOffScreen () {
        if (offScreen == null) {
            System.out.println(pixels);
            offScreen = createImage(pixels, pixels);
            osg = offScreen.getGraphics();
//            getParent().add(table);
            table.setBounds(this.getBounds());
        }
    }
        public void addMeter(Meter m) {
            if(m instanceof MeterPotentialEnergy && !showingPE) {return;}
            if(m instanceof MeterKineticEnergy && !showingKE) {return;}
            nMeters++;
            Meter[] temp = new Meter[nMeters];
            for(int i=0; i<meter.length; i++) {temp[i] = meter[i];}
            temp[nMeters-1] = m;
            meter = temp;
        }
        
        public void setPhase(Phase p) {
            for(Meter m=p.firstMeter; m!=null; m=m.nextMeter()) {this.addMeter(m);}
        }

    public void doUpdate() {
//        TableDataEvent event = new TableDataEvent(dataSource,0,0,0,0,TableDataEvent.RESET);
//        dataSource.fireTableDataEvent(event);
    }
    
    public void doPaint(Graphics g) {
        g.setColor(Color.red);
        g.drawRect(0,0,getSize().width-1,getSize().height-1);
        g.drawRect(1,1,getSize().width-3,getSize().height-3);
        table.paint(g);
    }
    
//    public void paint(Graphics g) {
//      createOffScreen();
//      doPaint(g);
//    }

    class MyTableData extends AbstractTableModel {
        
        private final String[] columnNames = new String[] {"Property", "Current", "Average"};
        private final Class[] columnClasses = {String.class, Double.class, Double.class};
        
        MyTableData() {}
        
        public Object getValueAt(int row, int column) {
            Meter m = meter[row];
            switch(column) {
                case 0: return m.getLabel();
                case 1: return new Double(m.currentValue());
                case 2: return new Double(m.average());
                default: return null;
            }
        }
        
        public int getRowCount() {return meter.length;}
        public int getColumnCount() {return 3;}  //label, current, average
                
        public String getColumnName(int column) {return columnNames[column];}
        public Class getColumnClass(int column) {return columnClasses[column];}
    }
}