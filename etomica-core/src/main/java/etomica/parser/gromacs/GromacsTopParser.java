package etomica.parser.gromacs;

import org.jparsec.Parser;
import org.jparsec.Parsers;
import org.jparsec.Scanners;
import org.jparsec.pattern.CharPredicate;
import org.jparsec.pattern.CharPredicates;
import org.jparsec.pattern.Pattern;
import org.jparsec.pattern.Patterns;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class GromacsTopParser {

    private static final Parser<Void> NEWLINE = Scanners.isChar('\n');
    private static final Parser<Void> COMMENT = Scanners.lineComment(";");
    private static final Pattern SPACE_PAT = Patterns.or(
            Patterns.among(" \t"),
            Patterns.string("\\\n") // escaped newline
    );
    private static final Parser<Void> SPACE = SPACE_PAT.many1().toScanner("space");
    private static final Parser<Void> SEPARATOR = Parsers.or(COMMENT, SPACE).skipMany();
    private static final Parser<Void> BLANK = Parsers.or(COMMENT, SPACE, NEWLINE).skipMany();
    private static final Parser<Void> LINE_SEP = NEWLINE.followedBy(BLANK);

    private static final CharPredicate ITEM_CHAR = c -> {
        return c != '[' && CharPredicates.range('!', '~').isChar(c);
    };

    private static final Parser<String> ITEM = Patterns.many1(ITEM_CHAR).toScanner("item").source();
    private static final Parser<List<String>> LINE = SEPARATOR.next(ITEM.sepEndBy1(SEPARATOR));

    private static final Parser<String> SECTION_HEADER = Parsers.between(
            Scanners.isChar('[').followedBy(SEPARATOR),
            Patterns.isChar(CharPredicates.IS_ALPHA_NUMERIC_).many1().toScanner("section name").source(),
            SEPARATOR.followedBy(Scanners.isChar(']'))
    );

    private static final Parser<GromacsTopSection> SECTION = Parsers.sequence(
            SECTION_HEADER.followedBy(LINE_SEP),
            LINE.sepBy(LINE_SEP),
            (header, lines) -> {
                List<String[]> sectionLines = lines.stream()
                        .map(list -> list.toArray(new String[0]))
                        .collect(Collectors.toList());
                return new GromacsTopSection(header, sectionLines);
            }
    );

    private static final Parser<List<GromacsTopSection>> FILE = BLANK.next(SECTION.sepEndBy1(BLANK));


    public static List<GromacsTopSection> parseContents(Readable contents) {
        try {
            return FILE.parse(contents);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    public static class GromacsTopSection {
        public final String name;
        public final List<String[]> lines;

        public GromacsTopSection(String name, List<String[]> lines) {
            this.name = name;
            this.lines = lines;
        }

        @Override
        public String toString() {
            return new StringJoiner(", ", GromacsTopSection.class.getSimpleName() + "{", "}")
                    .add("name='" + name + "'")
                    .add("lines=" + Arrays.deepToString(lines.toArray()))
                    .toString();
        }
    }

    public static void main(String[] args) throws IOException {
//        System.out.println(SECTION.parse("[ foo ]\n;blah\n  hi hi hi"));
        Path path = Paths.get("/home/alex/workspace/mosdef-test/ethane-box.top");
        System.out.println(FILE.parse(Files.newBufferedReader(path)));
    }

}
