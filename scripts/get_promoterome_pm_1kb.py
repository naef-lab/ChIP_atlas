import pandas as pd

infile='resources/mm10_promoters_v2.gff'
outfile='resources/mm10_promoters_v2_pm1kb.gff'

win = 1000
prom = pd.read_csv('resources/mm10_promoters_v2.gff',sep='\t',header=None,comment='#')

middle = ((prom.loc[:,4] + prom.loc[:,3])/2).astype(int)
# start
prom.loc[:,3] = middle - win
# clip start at 0
prom.loc[ prom.loc[:,3]<0,3] = 0
# end
prom.loc[:,4] = middle + win

prom.to_csv(outfile,sep='\t',header=False,index=False)
