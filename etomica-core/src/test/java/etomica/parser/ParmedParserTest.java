package etomica.parser;

import org.junit.Test;

import java.io.IOException;

import static etomica.parser.ParmedParser.parseGromacsResourceFiles;
import static org.junit.Assert.*;

public class ParmedParserTest {

    @Test
    public void testParseGromacsResourceFiles() {
        try {
            ParmedStructure p = parseGromacsResourceFiles("test.top", "test.gro");
        } catch(IOException e) {
            fail("Error processing resource files: " + e.getMessage());
        }
    }

}
