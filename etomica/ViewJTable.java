package simulate;

import simulate.*;
import java.beans.*;
import jclass.table3.*;
import java.awt.*;

public class ViewJTable extends simulate.View
{
    Table table;
    MyTableData dataSource;
    Meter meter;
    
    public ViewJTable()
    {
        table = new jclass.table3.Table();
        dataSource = new MyTableData();   //inner class, defined below
        table.setDataSource(dataSource);
    }
    
        public void setMeter(Meter m) {
            meter = m;
        }
        
    public void doUpdate() {
        TableDataEvent event = new TableDataEvent(dataSource,0,0,0,0,TableDataEvent.CHANGE_VALUE);
        dataSource.fireTableDataEvent(event);
    }
    
	public void setParentDisplay(Display p) {
	    super.setParentDisplay(p);
	    p.add(table);
	}
	 
    class MyTableData extends jclass.table3.TableDataSupport {
        
        MyTableData() {
        }
        
        public Object getTableDataItem(int row, int column) {
            if(column==0) {
                double value = meter.currentValue();
                return new Double(value);
            }
            else {
                return null;
            }
        }
        
        public int getNumRows() {return 3;}
        
        public int getNumColumns() {return 1;}
        
        public Object getTableRowLabel(int row) {
            return Integer.toString(row);
        }
        
        public Object getTableColumnLabel(int column) {
            return "Some Data";
        }
    }
}