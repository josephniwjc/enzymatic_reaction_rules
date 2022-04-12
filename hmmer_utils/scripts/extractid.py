import sys

in_path = sys.argv[-2]
out_path = sys.argv[-1]

with open(in_path) as f_input:
    with open(out_path, 'w') as f_output:
        for line in f_input.readlines():
            if line.startswith('#'):
                continue
            f_output.write(line.split(' ')[0] + '\n')
