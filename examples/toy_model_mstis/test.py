import openpathsampling as p

st = p.AnalysisStorage('mstis.nc')
st2 = p.Storage('mst.nc', 'w')

st2.save(st.snapshots[0])

q = st.snapshots.all().as_proxies()

cvs = st.cvs[:]

map(st2.save, cvs)

for cv in cvs:
    cv(q)

map(st.snapshots.mention, st.snapshots)
