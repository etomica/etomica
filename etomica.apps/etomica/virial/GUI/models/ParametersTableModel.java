package etomica.virial.GUI.models;

import javax.swing.table.AbstractTableModel;


import etomica.virial.GUI.components.CreateP22CLJQ;
import etomica.virial.GUI.components.CreateP2LJ;
import etomica.virial.GUI.components.CreateP2LJQ;

public class ParametersTableModel extends AbstractTableModel {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public String[] columnNames = {"Parameter",
            "Value"};
	
	private Object[][] data;
	
	
	public ParametersTableModel(){
		data = new String[12][2];
		
	}
	
	public Object[][] getData() {
		return data;
	}

	public boolean isCellEditable(int row, int col) {
        //Note that the data/cell address is constant,
        //no matter where the cell appears onscreen.
        if (col == 0) {
            return false;
        } else {
            return true;
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
		ParametersTableModel ljt = new ParametersTableModel();
		CreateP22CLJQ p2LJ = new CreateP22CLJQ();
		//ljt.UpdateObjectData(p2LJ.getParametersArray());
		//System.out.println(ljt.getData1());
	}
}
