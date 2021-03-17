def export_psc(psclist):
    with open('results/statistics/pscresults.txt',"w") as f:
            for row in psclist:
                f.write('%s %i %i %i\n' % (row[0],row[1][0],row[1][1],row[1][2]))
    