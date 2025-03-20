import json
import glob

def readInput(directory: str) -> list:
    pattern = '*.txt'

    seqs = []
    
    files = glob.glob(f'{directory}/{pattern}')
    for file in files:
        with open(file, 'r', encoding="utf8") as f:
            firstLine = f.readline() # read the first line and throw it out
            data = f.read()
            filename = file[file.index('\\')+1: file.index('.')]
            seqs.append((filename, data))
    
    return seqs

def buildCommonVariance(sequences: list, groupBy: list) -> list:
    MAX_CHARS = 1000
    groupedSeqs = {}
    bestVariances = {}
    for group in groupBy:
        groupedSeqs[group] = []
        bestVariances[group] = []
    groupedSeqs["unique"] = []
    bestVariances["unique"] = []

    for seq in sequences:
        start = seq[1][:4]
        matched = False
        for group in groupBy:
            if seq[1].startswith(group):
                matched = True
                groupedSeqs[group].append(seq)
                break
        if not matched:
            groupedSeqs["unique"].append(seq)

    commonVariances = []

    for key in groupedSeqs:
        commonVariance = [{} for _ in range(MAX_CHARS)] # char: freq 
        for seq in groupedSeqs[key]:
            for i, c in enumerate(seq[1]):
                commonVariance[i][c] = commonVariance[i].get(c, 0) + 1
        commonVariances.append((key, commonVariance))
    
    for group, dics in commonVariances:
        for dic in dics:
            total = 0
            best = None
            for char in dic:
                total += dic[char]
                if best == None or best[1] < dic[char]:
                    best = (char, dic[char])
            if not best == None:
                bestVariances[group].append((best[0], (best[1]/total)*100))
            
    return groupedSeqs, bestVariances


def output(directory: str, groupedSequences: list, bestVariances: list) -> None:
    with open(directory, 'w', encoding="utf8") as o:
        for group in groupedSequences:
            o.write('<h1><strong>'+group+'</strong></h1>\n<hr>\n')
            for seq in groupedSequences[group]:
                o.write('<h1>'+seq[0]+'</h1>\n')
                o.write('<p>')
                for i in range(len(seq[1])):
                    color = "red"
                    if bestVariances[group][i][0] == seq[1][i]:
                        color = "green"
                    if seq[1][i] == "K" or seq[1][i] == "R":
                        color = "blue"
                    o.write('<span style="color: '+color+'">'+seq[1][i]+'</span>')
                o.write('</p>\n\n')

if __name__ == '__main__':
    dir = "serovars"
    out = "output\\variance.html"
    groups = ["RLSS", "MAQV"] # Change me as needed :)
    sequences = readInput(dir)
    groupedSeqs, commonVariance = buildCommonVariance(sequences, groups)
    output(out, groupedSeqs, commonVariance)
    