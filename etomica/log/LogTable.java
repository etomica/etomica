package etomica.log;

import etomica.*;
import javax.swing.table.AbstractTableModel;
import org.apache.poi.hssf.usermodel.*;

/**
 * Writes a table to an excel-formatted file.
 */
public class LogTable {
    
    public LogTable() {
    }
    
    public void writeTable(AbstractTableModel table, String file) {
        
        short rownum;

        // create a new file
        java.io.FileOutputStream out = null;
        try {
            out = new java.io.FileOutputStream(file);
        } catch(java.io.FileNotFoundException ex) {ex.printStackTrace();}
        // create a new workbook
        HSSFWorkbook wb = new HSSFWorkbook();
        // create a new sheet
        HSSFSheet s = wb.createSheet();
        // declare a row object reference
        HSSFRow r = null;
        // declare a cell object reference
        HSSFCell c = null;
        
        int rowCount = table.getRowCount();
        int colCount = table.getColumnCount();
        
        for (rownum = (short) 0; rownum < rowCount; rownum++) {
            // create a row
            r = s.createRow(rownum);
            for (short cellnum = (short) 0; cellnum < colCount; cellnum++) {
                // create a numeric cell
                c = r.createCell(cellnum);
                // do some goofy math to demonstrate decimals
                c.setCellValue(table.getValueAt(rownum, cellnum).toString());
            }//end cell loop
        }//end row loop
        
        try {
            wb.write(out);
            out.close();
        } catch(java.io.IOException ex) {ex.printStackTrace();}
    }//end writeTable

}


