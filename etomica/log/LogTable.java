package etomica.log;

import javax.swing.table.AbstractTableModel;

/**
 * Writes a table to an excel-formatted file.
 */
public class LogTable {
    
    //if true, will write values without overwriting old ones
    private boolean append = false;
    private short r0 = 0;
    HSSFWorkbook wb = null;
    HSSFSheet s = null;

    
    public LogTable() {
    }
    
    public void setAppend(boolean b) {
        append = b;
        if(!append) r0 = 0;
    }
    public boolean isAppend() {return append;}
    
    public void writeTable(AbstractTableModel table, String file) {
        
        short rownum;
        // create a new file
        java.io.FileOutputStream out = null;
        try {
            out = new java.io.FileOutputStream(file);
        } catch(java.io.FileNotFoundException ex) {ex.printStackTrace();}
        // create a new workbook
        if(wb == null || !append) wb = new HSSFWorkbook();
        // create a new sheet
        if(s == null || !append) s = wb.createSheet();
        // declare a row object reference
        HSSFRow r = null;
        // declare a cell object reference
        HSSFCell c = null;
        
        int rowCount = table.getRowCount();
        int colCount = table.getColumnCount();
        
        for (rownum = (short) 0; rownum < rowCount; rownum++) {
            // create a row
            r = s.createRow((short)(r0+rownum));
            for (short cellnum = (short) 0; cellnum < colCount; cellnum++) {
                // create a numeric cell
                c = r.createCell(cellnum);
                // do some goofy math to demonstrate decimals
                c.setCellValue(table.getValueAt(rownum, cellnum).toString());
            }//end cell loop
        }//end row loop
        if(append) r0 += rowCount + 1;
        
        try {
            wb.write(out);
            out.close();
        } catch(java.io.IOException ex) {ex.printStackTrace();}
    }//end writeTable

}


