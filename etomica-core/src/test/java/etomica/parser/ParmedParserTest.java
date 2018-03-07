/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.parser;

import org.junit.jupiter.api.Test;

import java.io.IOException;

import static etomica.parser.ParmedParser.parseGromacsResourceFiles;
import static org.junit.jupiter.api.Assertions.fail;
import static org.junit.jupiter.api.Assumptions.*;

class ParmedParserTest {

    @Test
    void testParseGromacsResourceFiles() {
        assumeTrue(ParmedParser.class.getClassLoader().getResource("virtualenv/bin/parmed_json") != null);
        try {
            ParmedStructure p = parseGromacsResourceFiles("test.top", "test.gro");
        } catch(IOException e) {
            fail("Error processing resource files: " + e.getMessage());
        }
    }

}
