/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import com.jgoodies.looks.LookUtils;

public class TableView {

  public static JTable build() {

    TableModel model = new SampleTableModel(createSampleTableData(),
        new String[] { "ID", "Score      ", "Marked?" });
    JTable table = new JTable(model);
    table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    table.getColumnModel().getColumn(0).setPreferredWidth(100);
    table.getColumnModel().getColumn(1).setPreferredWidth(200);
    table.getColumnModel().getColumn(2).setPreferredWidth(100);
    table.setRowSelectionInterval(2, 2);
    int tableFontSize = table.getFont().getSize();
    int minimumRowHeight = tableFontSize + 6;
    int defaultRowHeight = LookUtils.IS_LOW_RESOLUTION ? 17 : 18;
    table.setRowHeight(Math.max(minimumRowHeight, defaultRowHeight));
    return table;
  }

  private static Object[][] createSampleTableData() {

    return new Object[][] {
        { "Albert Ayler", "Greenwich Village", Boolean.TRUE },
        { "Carla Bley", "Escalator Over the Hill", Boolean.TRUE },
        { "Frank Zappa", "Yo' Mama", Boolean.TRUE },
        { "John Coltrane", "Ascension", Boolean.TRUE },
        { "Miles Davis", "In a Silent Way", Boolean.TRUE },
        { "Pharoa Sanders", "Karma", Boolean.TRUE },
        { "Wayne Shorter", "Juju", Boolean.TRUE }, { "", "", Boolean.FALSE },
        { "", "", Boolean.FALSE }, { "", "", Boolean.FALSE },
        { "", "", Boolean.FALSE }, { "", "", Boolean.FALSE }, };
  }

  private static final class SampleTableModel extends AbstractTableModel {

    private static final long serialVersionUID = 6802574060087202171L;
    private final String[]    columnNames;
    private final Object[][]  rowData;

    SampleTableModel(Object[][] rowData, String[] columnNames) {

      this.columnNames = columnNames;
      this.rowData = rowData;
    }

    public String getColumnName(int column) {

      return columnNames[column].toString();
    }

    public int getRowCount() {

      return rowData.length;
    }

    public int getColumnCount() {

      return columnNames.length;
    }

    public Class getColumnClass(int column) {

      return column == 2 ? Boolean.class : super.getColumnClass(column);
    }

    public Object getValueAt(int row, int col) {

      return rowData[row][col];
    }

    public boolean isCellEditable(int row, int column) {

      return true;
    }

    public void setValueAt(Object value, int row, int col) {

      rowData[row][col] = value;
      fireTableCellUpdated(row, col);
    }
  }
}