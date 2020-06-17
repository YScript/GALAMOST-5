import gamst
	
mst = gamst.snapshot.read("lj.mst")
app = gamst.application.dynamics(info=mst, dt=0.001)

fn = gamst.force.dpd(info=mst, rcut=1.0)
fn.setParams(type_i="a", type_j="a", alpha=30.0, sigma=3.0)
app.add(fn)

inn = gamst.integration.gwvv(info=mst, group=['a'])
app.add(inn)

di = gamst.dump.data(info=mst, group=['a'], file='data.log', period=100)
app.add(di)

app.run(1000)