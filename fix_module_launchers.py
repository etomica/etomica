import sys

filename = sys.argv[1]

with open(filename) as f:
    contents = f.read()

contents = contents.replace('$0', '"$0"')
contents = contents.replace('$DIR', '"$DIR"')

with open(filename, 'w') as f:
    f.write(contents)