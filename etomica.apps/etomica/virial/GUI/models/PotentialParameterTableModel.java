package etomica.virial.GUI.models;

import java.util.EventObject;

import javax.swing.CellEditor;
import javax.swing.event.CellEditorListener;
import javax.swing.table.AbstractTableModel;



public class PotentialParameterTableModel extends AbstractTableModel  {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public String[] columnNames = {"Parameter",
            "Value"};
	
	private Object[][] data;
	
	
	public PotentialParameterTableModel(){
		data = new String[12][2];
		
	}
	
	public Object[][] getData() {
		return data;
	}

	public boolean isCellEditable(int row, int col) {
        //Note that the data/cell address is constant,
        //no matter where the cell appears onscreen.
        if (col == 1) {
            return true;
        } else {
            return false;
        }
		
    }
	
	public void removeData(){
		for (int i = 0; i<12;i++){
			for (int j = 0;j<2;j++){
				if(data[i][j] != null){
				data[i][j] = null;
				fireTableRowsDeleted(i,j);
				}
			}
		}
	}
	
	
	 public void setValueAt(Object value, int row, int col) {
		 data[row][col] = value;
	     fireTableCellUpdated(row, col);
	  }
	
	@Override
	public int getRowCount() {
		// TODO Auto-generated method stub
		return 12;
	}

	@Override
	public int getColumnCount() {
		// TODO Auto-generated method stub
		return 2;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		// TODO Auto-generated method stub
		return data[rowIndex][columnIndex];
	}

	
	public static void main(String[] args){
		PotentialParameterTableModel ljt = new PotentialParameterTableModel();
		CreateSpeciesDM_LJ_2CLJQ p2LJ = new CreateSpeciesDM_LJ_2CLJQ();
		//ljt.UpdateObjectData(p2LJ.getParametersArray());
		//System.out.println(ljt.getData1());
	}

	

	
}
