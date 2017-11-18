import os
import re

leftoverFile = "/Users/Nate/Desktop/pounds_left_after_machine_write"
leftoverLines = []
with open(leftoverFile, 'rb') as inputFile:
    for line in inputFile:
        leftoverLines.append(line)

print leftoverLines


originalPoundLines = []
for root, dirs, files in os.walk("/Users/Nate/code/RMG-database/input"):
    for file in files:
        if not file == '.DS_Store':
            # print os.path.join(root, file)
            filepath = os.path.join(root,file)
            lines = []
            with open(filepath) as inputFile:
                for line in inputFile:
                    lines.append(line)

            for index, line in enumerate(lines):
                if re.search('\#', line) and index > 1:
                    if filepath + ' ' + line not in leftoverLines:
                    # if True:
                    #     print filepath, line
                        print filepath, index, line

