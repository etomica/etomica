package simulate;

import simulate.*;
import java.beans.*;
import jclass.table3.*;
import java.awt.*;

public class ViewJTable extends simulate.View
{
    Table table;
    MyTableData dataSource;
    
    public ViewJTable()
    {
        table = new jclass.table3.Table();
        dataSource = new MyTableData();   //inner class, defined below
        table.setDataSource(dataSource);
    }
    
    public void doUpdate() {;}
    
    public void paint(Graphics g) {table.paint(g);}

    class MyTableData extends jclass.table3.TableDataSupport {
        
        MyTableData() {
        }
        
        public Object getTableDataItem(int row, int column) {
            if(column==0) {
 //               double value = parentDisplay.getMeter(row).currentValue(parentDisplay.phase);
//                return new Double(value);
                    return null;
            }
            else {
                return null;
            }
        }
        
        public int getNumRows() {return 1;}
        
        public int getNumColumns() {return 1;}
        
        public Object getTableRowLabel(int row) {
            return Integer.toString(row);
        }
        
        public Object getTableColumnLabel(int column) {
            return "Some Data";
        }
        
        public void update(Graphics g) {;}
    }
}