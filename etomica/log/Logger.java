package etomica.log;

import javax.swing.table.AbstractTableModel;

/* History of changes
   8/26/02 (DAK) new
*/

public interface Logger {
    
    public void write(AbstractTableModel table);
    public void append(AbstractTableModel table);
    
}