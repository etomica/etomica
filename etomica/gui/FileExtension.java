
package simulate.gui;

import java.io.*;

/**
 * Filter class
 */

public class FileExtension implements FilenameFilter {
    private String extension;

    public FileExtension(String ext) {
	this.extension = ext;
    }
    
    public boolean accept(File dir, String name) {
	return name.endsWith(extension);
    }
}
