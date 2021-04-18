package etomica.species;

import etomica.atom.AtomType;
import etomica.space.Vector;
import etomica.space3d.Space3D;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class SpeciesFileParser {
    private static final Pattern SECTION_HEADER = Pattern.compile("\\w+:");

    public static SpeciesGeneral parsePath(String path) throws IOException {
        return parse(Files.newBufferedReader(Paths.get(path)));
    }

    public static SpeciesGeneral parse(Readable contents) {
        Scanner scanner = new Scanner(contents)
                .skip("#.*\\R"); // ignore comments starting with '#' TODO this is totally wrong...

        List<AtomType> types = parseTypesSection(scanner);
        Map<String, AtomType> typesMap = types.stream().collect(Collectors.toMap(
                atomType -> atomType.getElement().getSymbol(),
                Function.identity()
        ));

        List<AtomDecl> atoms = parseAtomsSection(scanner, typesMap);

        // TODO: orientation, constraints, check for eof

        SpeciesBuilder builder = new SpeciesBuilder(Space3D.getInstance()).setDynamic(true);
        for (AtomDecl atom : atoms) {
            builder.addAtom(atom.type, atom.position, atom.name);
        }

        return builder.build();
    }

    private static List<AtomType> parseTypesSection(Scanner scanner) {
        String header = scanner.next();
        if (!header.equals("types:")) {
            // throw
        }
        List<AtomType> types = new ArrayList<>();
        while (!scanner.hasNext(SECTION_HEADER)) {
            AtomType type = parseAtomType(scanner);
            types.add(type);
        }
        return types;
    }

    private static AtomType parseAtomType(Scanner scanner) {
        String symbol = scanner.next();
        double mass = scanner.nextDouble();
        AtomType type = AtomType.simple(symbol, mass);
        String extra = scanner.nextLine();
        if (!extra.isEmpty()) {
            // throw
        }
        return type;
    }

    private static List<AtomDecl> parseAtomsSection(Scanner scanner, Map<String, AtomType> types) {
        String header = scanner.next();
        if (!header.equals("atoms:")) {
            // throw
        }

        List<AtomDecl> atoms = new ArrayList<>();
        int idx = 1;
        while(!scanner.hasNext(SECTION_HEADER)) {
            AtomDecl atom = parseAtom(scanner, types, idx);
            atoms.add(atom);
        }
        return atoms;
    }

    private static AtomDecl parseAtom(Scanner scanner, Map<String, AtomType> types, int expectedIdx) {
        int idx = scanner.nextInt();
        if (!(idx == expectedIdx)) {
            // throw
        }
        String typeStr = scanner.next();
        if (!types.containsKey(typeStr)) {
            // throw
        }
        AtomType type = types.get(typeStr);
        // TODO name
        double x = scanner.nextDouble();
        double y = scanner.nextDouble();
        double z = scanner.nextDouble();
        String extra = scanner.nextLine();
        if (!extra.isEmpty()) {
            // throw
        }
        AtomDecl atom = new AtomDecl();
        atom.type = type;
        atom.name = "";
        atom.position = Vector.of(x, y, z);
        return atom;
    }

    private static class AtomDecl {
        AtomType type;
        String name;
        Vector position;
    }

    public static void main(String[] args) throws IOException {
        SpeciesGeneral species = parsePath("water.species");
        System.out.println(species);
    }
}
