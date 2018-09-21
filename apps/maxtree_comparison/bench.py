# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 21:48:35 2013

@author: edwin
"""
import csv, sqlite3
import numpy as np

class MyStdv:
    def __init__(self):
        self.count = 0
        self.s = 0.0
        self.s2 = 0.0
        
    def step(self, x):
        self.count += 1
        self.s += x
        self.s2 += x*x
        
    def finalize(self):
        return np.sqrt((self.s2 / self.count) - (self.s / self.count) ** 2)

con = sqlite3.connect(":memory:")
con.create_aggregate("stdv", 1, MyStdv)

cur = con.cursor()
cur.execute("CREATE TABLE t (algo TEXT, image TEXT, nbits INTEGER, nthread INTEGER, size INTEGER, time REAL);")
cur.execute("CREATE TABLE algos (serial TEXT, parallel TEXT);") # serial <=> parallel

with open('log.txt','rb') as fin:
    # csv.DictReader uses first line in file for column headings by default
    dr = csv.reader(fin)
    dr.next() # skip first line
    to_db = [i for i in dr]

cur.executemany("INSERT INTO t (algo, image, nbits, nthread, size, time) VALUES (?, ?, ?, ?, ?, ?);", to_db)
cur.executemany("INSERT INTO algos (serial, parallel) VALUES (?, ?);",
                [("maxtree_serial_hqueue", "maxtree_parallel_hqueue"),
                 ("maxtree_serial_nister", "maxtree_parallel_pqueue"),
                 ("maxtree_parallel_ufind_line", "maxtree_parallel_ufind_line"),
                 ("maxtree_serial_ufindrank", "maxtree_parallel_ufindrank"),
                 ("maxtree_serial_ufind_wlc", "maxtree_parallel_ufind")])
                 


con.commit()

column2labels = {
    "size" : "Number of pixels $(\\times 10^6)$",
    "mean" : "Time (ms)",
    "speedup" : "Speed Up",
    "nbits" : "Number of bits",
    "nthread" : "Number of threads"
}


algos = {
    "maxtree_serial_berger" : ("Berger et al. [11]", "b:x"),
    "maxtree_serial_nister" : (u"Nistér and Stewénius [24]", "k--o"),
    "maxtree_serial_ufindrank" : ("Berger + rank", "g-.+"),
    "maxtree_serial_hqueue" : ("Salembier et al. [3]", "r-s"),
    "maxtree_serial_ufind_wlc" : ("Berger + level compression", "y-d"),
    "maxtree_serial_najman" : ("Najman and Couprie [14]", "m:+"),
    "maxtree_serial_wilkinson" : ("Wilkinson [19]", "c-+"),
    "maxtree_parallel_hqueue" : ("Salembier et al. [3]", "r-s"),
    "maxtree_parallel_pqueue" : ("Salembier unrecursive", "c-.+"),
    "maxtree_parallel_ufind_line": ("Matas et al. [20]", "k:*"),
    "maxtree_parallel_ufindrank": ("Berger + rank", "g-.+"),
    "maxtree_parallel_ufind" : ("Berger + level compression", "y-d")
}

def plotdata(cur):
    data = {} # algo -> ([x], [y], [err])        
    for row in cur:
        if not data.has_key(row[0]):
            data[row[0]] = ([],[],[])
        d = data[row[0]]
        d[0].append(row[1])
        d[1].append(row[2])
        d[2].append(row[3])
    clf()
    for k, d in data.iteritems():
        lbl, mark = algos[k]
        plot(d[0], d[1], mark, label=lbl)
        errorbar(d[0], d[1], d[2], fmt=mark[0] + mark[-1])
    xlabel(column2labels[cur.description[1][0]])
    ylabel(column2labels[cur.description[2][0]])
    

from matplotlib.ticker import FormatStrFormatter
    
######################################
# Display by Seq Algorithm by size   #
######################################
rq = '''
SELECT algo, size / 1000000 as size, AVG(time*1000) as mean, stdv(time*1000) / 2 as std
FROM t
WHERE nthread = -1 AND nbits = 8
GROUP BY algo, size;
'''
figure()
plotdata(cur.execute(rq))
pl = gca()
pl.set_xscale("log", basex=2)
pl.set_yscale("log", basey=2)
pl.set_xticks(2**arange(10))
pl.set_yticks(  50*2**arange(14) )
pl.yaxis.set_major_formatter(FormatStrFormatter('%d'))
pl.xaxis.set_major_formatter(FormatStrFormatter('%d'))
pl.set_xlim(1,256)
l = legend(loc=2, prop = {"size" : "medium"})
l.get_frame().set_visible(False)
draw()
savefig("res-fig1.pdf")


######################################
# Display by Seq Algorithm by nbits  #
######################################
rq = '''
SELECT algo, nbits, AVG(time*1000) as mean, stdv(time*1000) / 2 as std
FROM t
WHERE nthread = -1 AND size = 8388608
GROUP BY algo, nbits;
'''

figure()
plotdata(cur.execute(rq))
pl = gca()
pl.set_xlim(8,32)
pl.set_xticks(range(8,34,2))
pl.set_ylim(0,3500)
pl.set_yticks(range(0,3501,500) )
l = legend(loc=4, prop = {"size" : "medium"})
l.get_frame().set_visible(False)
draw()
savefig("res-fig2.pdf")

######################################
# Display by Number of threads  #
######################################
rq = '''
SELECT algo, nthread, AVG(time*1000) as mean, stdv(time*1000) / 2 as std
FROM t
WHERE nthread != -1 AND size = 8388608 AND nbits = 8
GROUP BY algo, nthread;
'''

figure()
plotdata(cur.execute(rq))
pl = gca()
pl.set_xlim(0.9,16)
pl.set_xticks(range(2,17,2))
pl.set_ylim(100,750)
l = legend(loc=1, prop = {"size" : "medium"})
l.get_frame().set_visible(False)
draw()
savefig("res-fig3.pdf")

###########################################
# Display by Number of threads (Speed up) #
###########################################
rq = '''
SELECT t2.algo, t2.nthread, avg(t1.time / t2.time) as speedup, stdv(t1.time / t2.time) / 2 as std
FROM t as t1
INNER JOIN algos ON t1.algo = algos.serial
INNER JOIN t as t2 ON t2.algo = algos.parallel AND t1.image = t2.image AND t1.size = t2.size AND t1.nbits = t2.nbits
WHERE t1.nthread <= 1 and t1.nbits = 8 AND t1.size = 8388608
GROUP BY t2.algo, t2.nthread;
'''

#for row in cur.execute(rq):
#    print row

figure()
plotdata(cur.execute(rq))
pl = gca()
pl.set_xlim(0.95,16)
pl.set_xticks(range(2,17,2))
pl.set_ylim(0.5,4)
l = legend(loc=4, prop = {"size" : "medium"})
l.get_frame().set_visible(False)
draw()
savefig("res-fig4.pdf")

######################################
# Display Parallel by Nbits  #
######################################
rq = '''
SELECT algo, nbits, AVG(time*1000) as mean, stdv(time*1000) / 2 as std
FROM t
WHERE nthread != -1 AND size = 8388608 AND nthread = 8
GROUP BY algo, nbits;
'''

figure()
plotdata(cur.execute(rq))
pl = gca()
pl.set_xlim(8,18)
pl.set_ylim(0,2500)
l = legend(loc=1, prop = {"size" : "medium"})
l.get_frame().set_visible(False)
draw()
savefig("res-fig5.pdf")
