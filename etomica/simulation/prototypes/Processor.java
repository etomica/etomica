package etomica.simulation.prototypes;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
/**
 * @author mimi
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class Processor extends JFrame {
	private File fileName;

	private Object input;

	private FileInputStream fileInput;

	public static void main(String[] args) {
		final Processor app = new Processor();
		
		app.addWindowListener(
			new WindowAdapter() {
					public void windowClosing(WindowEvent e) {
//							app.closeFile();
							System.exit(0);
					}
			}
			);
		app.openFile(true);
		File fName = app.fileName;
		System.out.println("You chose the file: "+fName.getName());
		BufferedReader in = null;
//		DataInputStream in = null;
		try {
			in = new BufferedReader(new FileReader( fName ));
//			in = new DataInputStream(new BufferedInputStream(new FileInputStream( fName )));
		} catch(IOException ex) {}
		
		int nTemp = 5,
		    nLines = 40000;
		int [][] phase = new int [ nTemp][nLines];
		double[] t = new double[nTemp];
	//read temperatures
		try {
			String line = in.readLine();
			System.out.println(line);
			StringTokenizer st = new StringTokenizer(line);
			int k = 0;
			while(st.hasMoreTokens()) {
				t[k++] = Double.parseDouble(st.nextToken());
				System.out.println(t[k-1]);
			}
			System.out.println();
		} catch(IOException ex) {}

		for ( int i=0; i< 10; i++) {
//			System.out.println(i);
		try {
//			int x = dataInput.readInt();
			int[] x = new int[nTemp];
			String line = in.readLine();
			System.out.println(line);
			StringTokenizer st = new StringTokenizer(line);
			int k = 0;
			while(st.hasMoreTokens()) {
				x[k++] = Integer.parseInt(st.nextToken());
				System.out.println(x[k-1]);
			}
			System.out.println();
			
		}
		catch(IOException ex) {System.out.println("EOF");break;}
	}
	}
	private void openFile( boolean firstTime ) 
	{
	
	if ( firstTime ) { 
	
		JFileChooser fileChooser = new JFileChooser();
		
		fileChooser.setFileSelectionMode (JFileChooser.FILES_ONLY);
		int result = fileChooser.showOpenDialog( this );
		
		// user clicked Cancel buttom on dialog
		if ( result == JFileChooser.CANCEL_OPTION )
		return;
		
		fileName = fileChooser.getSelectedFile();
	}
//	if ( fileName == null || fileName.getName().equals( "" ) )
//					JOptionPane.showMessageDialog( this,
//					"Invalid File Name",
//					"Invalid File Name",
//					JOptionPane.ERROR_MESSAGE );
//	else {
//		// Open the file
//		try {
//			// close file from previous operatiom
////			if ( input != null )
////				input.close();
//			
////			input = new ObjectInputStream( fileInput );
//
//		}
//		catch ( IOException e ) {
//			JOptionPane.showMessageDialog( this,
//			"File does not exist",
//			"Invalid File Name",
//			JOptionPane.ERROR_MESSAGE );
//		}
//	}
}
		
}	
	

