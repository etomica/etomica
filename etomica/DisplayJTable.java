package simulate;

import simulate.*;
import java.beans.*;
import jclass.table3.*;
import java.awt.*;

public class DisplayJTable extends simulate.Display
{
    public Table table;
    MyTableData dataSource;
    Meter[] meter = new Meter[0];
    int nMeters = 0;
    
    public DisplayJTable()
    {
        super();
        table = new jclass.table3.Table();
        dataSource = new MyTableData();   //inner class, defined below
        table.setDataSource(dataSource);
//        setVisible(false);  //get out of way of painting of table
        displayTool = table;
    }
    
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
            nMeters++;
            Meter[] temp = new Meter[nMeters];
            for(int i=0; i<meter.length; i++) {temp[i] = meter[i];}
            temp[nMeters-1] = m;
            meter = temp;
        }
        
        public void setPhase(Phase p) {
            for(Meter m=p.firstMeter; m!=null; m=m.getNextMeter()) {this.addMeter(m);}
        }

    public void doUpdate() {
        TableDataEvent event = new TableDataEvent(dataSource,0,0,0,0,TableDataEvent.RESET);
        dataSource.fireTableDataEvent(event);
    }
    
    public void doPaint(Graphics g) {
//        g.setColor(Color.red);
//        g.drawRect(0,0,getSize().width-1,getSize().height-1);
//        g.drawRect(1,1,getSize().width-3,getSize().height-3);
//        table.paint(g);
    }
    
    public void paint(Graphics g) {
      createOffScreen();
      doPaint(g);
    }

    class MyTableData extends jclass.table3.TableDataSupport {
        
        MyTableData() {
        }
        
        public Object getTableDataItem(int row, int column) {
            if(row<meter.length) {
                double value = (column==0) ? meter[row].currentValue() : meter[row].average();
                return new Double(value);
            }
            else {
                return null;
            }
        }
        
        public int getNumRows() {return meter.length;}
        
        public int getNumColumns() {return 2;}
        
        public Object getTableRowLabel(int row) {
            return meter[row].getLabel();
 //           return Integer.toString(row);
        }
        
        public Object getTableColumnLabel(int column) {
            return (column==0) ? "Current" : "Average";
        }
    }
}