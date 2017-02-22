/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import java.util.EventObject;

import javax.swing.CellEditor;
import javax.swing.event.CellEditorListener;
import javax.swing.table.AbstractTableModel;



public class ModelMoleculeListDisplayTable extends AbstractTableModel  {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public String[] columnNames = {"Parameter",
            "Value"};
	
	private Object[][] data;
	
	
	public ModelMoleculeListDisplayTable(){
		data = new String[12][2];
		
	}
	
	public Object[][] getData() {
		return data;
	}

	public boolean isCellEditable(int row, int col) {
        //Note that the data/cell address is constant,
        //no matter where the cell appears onscreen.

            return false;
		
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
		ModelMoleculeListDisplayTable ljt = new ModelMoleculeListDisplayTable();
		MolecularModel2CLJQ_SpeciesLJ p2LJ = new MolecularModel2CLJQ_SpeciesLJ();
		//ljt.UpdateObjectData(p2LJ.getParametersArray());
		//System.out.println(ljt.getData1());
	}

	

	
}
