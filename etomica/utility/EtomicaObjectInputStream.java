package etomica.utility;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.util.HashMap;
import java.util.LinkedList;
/*
 * Created on Jul 6, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class EtomicaObjectInputStream extends ObjectInputStream {

    public LinkedList atomLists = new LinkedList();
    public HashMap linkerLists = new HashMap();
    
    /**
     * @param out
     * @throws java.io.IOException
     */
    public EtomicaObjectInputStream(InputStream out) throws IOException {
        super(out);
    }

    /**
     * @throws java.io.IOException
     * @throws java.lang.SecurityException
     */
    public EtomicaObjectInputStream() throws IOException, SecurityException {
        super();
    }

}
